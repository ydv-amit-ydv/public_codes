"""NetMHCpan 4.2 runner: local Flask GUI for MHC-I binding affinity / eluted ligand.

Paste a FASTA / protein / peptide list, pick alleles (IEDB 27, supertypes, mouse
H-2 presets or any allele NetMHCpan knows) and lengths, and get a sortable
table of EL / BA scores and binder calls.

NetMHCpan itself is DTU-licensed and is NOT in this repository. Point
models.netmhcpan.home in engine_config.json at your own install (see README.md).
"""
import os, re, shlex, shutil, subprocess, sys, tempfile, threading
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache

import pandas as pd
from flask import Flask, jsonify, request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, discover, predcache, alleles as allele_lib   # noqa: E402
from common.jobs import run_cmd                              # noqa: E402

AA = set("ACDEFGHIKLMNPQRSTVWY")
ALLELE_RE = re.compile(r"^[A-Za-z0-9:\-_]+$")

# Presets in NetMHCpan spelling; each is filtered against `netMHCpan -listMHC` at runtime.
PRESETS = {k: [allele_lib.netmhcpan_name(x) for x in v["alleles"]] for k, v in allele_lib.PRESETS.items()}


def netmhc_path():
    """The netMHCpan script to run; a patched copy if its NMHOME placeholder was never edited."""
    home = config.model("netmhcpan")["home"]
    if not home:
        return os.path.join(config.search_root(), "netMHCpan-not-found")
    return discover.runnable_netmhcpan(home, config.tmp_dir())


def nm_env():
    """Clean env: system tcsh, not conda's (avoids the libtinfo warning); settings re-read per call."""
    env = {k: v for k, v in os.environ.items() if k not in ("LD_LIBRARY_PATH", "PYTHONPATH", "PYTHONHOME")}
    env["PATH"] = "/usr/bin:/bin"
    env["TMPDIR"] = config.tmp_dir()
    return env


app = Flask(__name__)


@lru_cache(maxsize=4)
def _list_alleles(path):
    try:
        r = subprocess.run([path, "-listMHC"], capture_output=True, text=True, env=nm_env())
    except OSError:                       # not installed / wrong path in Settings
        return []
    return [l.split()[0] for l in r.stdout.splitlines() if l.strip() and not l.startswith("#")]


def list_alleles():
    return _list_alleles(netmhc_path())


def to_num(s):
    try:
        return pd.to_numeric(s)
    except (ValueError, TypeError):
        return s


def parse_netmhcpan(stdout):
    """Parse every result table using its own header line."""
    rows, cols = [], None
    for line in stdout.splitlines():
        t = line.split()
        if not t:
            continue
        if t[0] == "Pos":
            cols = t
        elif cols and t[0].isdigit() and len(t) >= len(cols) - 1:
            n = len(cols) - 1                       # BindLevel is optional and two tokens
            row = dict(zip(cols[:n], t[:n]))
            row["BindLevel"] = " ".join(t[n:]).replace("<= ", "")
            rows.append(row)
    return pd.DataFrame(rows)


def read_sequences(text):
    """FASTA or bare lines -> [(name, seq)]."""
    seqs, name, buf = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if buf:
                seqs.append((name or f"seq{len(seqs)+1}", "".join(buf)))
            parts = line[1:].split()
            name, buf = (parts[0] if parts else None), []
        elif line:
            buf.append(line.upper())
    if buf:
        seqs.append((name or f"seq{len(seqs)+1}", "".join(buf)))
    return seqs


def build_peptides(text, mode, lengths, bounds=(8, 14)):
    """DataFrame[Protein, Start, Peptide] plus a list of warnings. `bounds` = accepted length of loose peptides."""
    warns, rows = [], []
    if mode == "peptides":
        for i, tok in enumerate(text.replace(">", " ").split(), 1):
            p = tok.upper()
            if not (bounds[0] <= len(p) <= bounds[1]) or set(p) - AA:
                warns.append(f"skipped '{tok}' (needs {bounds[0]}-{bounds[1]} standard aa)")
                continue
            rows.append(("peptide", i, p))
    else:
        for name, seq in read_sequences(text):
            bad = set(seq) - AA
            if bad:
                warns.append(f"{name}: non-standard residues {''.join(sorted(bad))}; windows containing them skipped")
            for L in lengths:
                for s in range(len(seq) - L + 1):
                    p = seq[s:s + L]
                    if not set(p) - AA:
                        rows.append((name, s + 1, p))
    return pd.DataFrame(rows, columns=["Protein", "Start", "Peptide"]), warns


def run_netmhcpan(peptides, alleles, p, job=None, stage="netmhcpan"):
    """run_netmhcpan_raw behind the prediction cache: only peptide-allele predictions not computed before are run.
    Returns (DataFrame, [log]) exactly like the raw call."""
    logs = []

    def compute(peps, als):
        df, lg = run_netmhcpan_raw(peps, als, p, job, stage)
        logs.extend(lg)
        return df
    sig = ("netmhcpan", netmhc_path(), bool(p["ba"]), p["rankS"], p["rankW"], p.get("t"), p.get("raw", ""))
    df = predcache.fetch("netmhcpan", sig, peptides, list(alleles), compute, allele_lib.canonical, job, stage)
    return df, logs


def run_netmhcpan_raw(peptides, alleles, p, job=None, stage="netmhcpan"):
    """Run NetMHCpan -p. Work is split into tasks of (one peptide length, a few alleles)
    that run in parallel; `job` (optional) gets progress and can cancel. Returns (DataFrame, [log])."""
    cfg = config.model("netmhcpan")
    perf = config.performance()
    workers = cfg.get("workers") or perf["workers"]
    by_len = defaultdict(list)
    for pep in sorted(set(peptides)):
        by_len[len(pep)].append(pep)          # -p needs one length per run
    # enough tasks to keep every worker busy: fewer alleles per call when there are many cores and few (length, allele) combinations
    per_call = max(1, min(cfg.get("alleles_per_call") or 3, -(-len(by_len) * len(alleles) // (2 * workers)) or 1))
    env = nm_env()
    workdir = tempfile.mkdtemp(dir=env["TMPDIR"])
    tasks = []
    for L, peps in sorted(by_len.items()):
        fp = os.path.join(workdir, f"len{L}.pep")
        with open(fp, "w") as f:
            f.write("\n".join(peps) + "\n")
        for i in range(0, len(alleles), per_call):
            tasks.append((L, fp, alleles[i:i + per_call]))
    script = netmhc_path()
    done = [0]
    lock = threading.Lock()
    if job:
        job.update(stage, done=0, total=len(tasks), msg=f"{len(tasks)} runs, {workers} in parallel")

    def one(t):
        L, fp, als = t
        if job:
            job.check()
        cmd = [script, "-p", "-f", fp, "-a", ",".join(als),
               "-rankS", str(p["rankS"]), "-rankW", str(p["rankW"])]
        if p["ba"]:
            cmd.append("-BA")
        if p.get("t") not in (None, ""):
            cmd += ["-t", str(p["t"])]
        if p.get("raw", "").strip():
            cmd += shlex.split(p["raw"])
        r = run_cmd(cmd, env=env, job=job, task={"model": "NetMHCpan", "what": f"{L}-mers \u00b7 {', '.join(a.replace('HLA-', '') for a in als)}",
                                                 "n": len(by_len[L]), "sample": by_len[L][:3]})
        with lock:
            done[0] += 1
            if job:
                job.update(stage, done=done[0], msg=f"{L}-mers: {done[0]}/{len(tasks)} runs finished")
        return t, cmd, r

    frames, logs, failed = [], [], []
    try:
        with ThreadPoolExecutor(max_workers=workers) as ex:
            for t, cmd, r in ex.map(one, tasks):
                logs.append("$ " + " ".join(cmd) + "\n" + r.stdout[-3000:] + r.stderr[-1500:])
                df = parse_netmhcpan(r.stdout)
                if df.empty:
                    failed.append(f"{t[0]}-mers {','.join(t[2])}")
                else:
                    frames.append(df)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
    if job and failed:
        job.log(f"NetMHCpan returned nothing for: {'; '.join(failed[:5])}" + (" ..." if len(failed) > 5 else ""))
    if not frames:
        return pd.DataFrame(), logs
    df = pd.concat(frames, ignore_index=True).rename(columns={"MHC": "Allele"})
    df["Allele"] = df["Allele"].str.replace("*", "", regex=False)
    for c in df.columns:
        if c not in ("Allele", "Peptide", "Core", "Icore", "Identity", "BindLevel"):
            df[c] = df[c].map(to_num)
    keep = ["Allele", "Peptide", "Core", "Score_EL", "%Rank_EL", "Score_BA", "%Rank_BA",
            "Aff(nM)", "BindLevel"]
    df = df[[c for c in keep if c in df.columns]]
    return df.rename(columns={"Score_EL": "NMP_Score_EL", "%Rank_EL": "NMP_%Rank_EL",
                              "Score_BA": "NMP_Score_BA", "%Rank_BA": "NMP_%Rank_BA",
                              "Aff(nM)": "NMP_Aff_nM", "BindLevel": "NMP_BindLevel"}), logs


@app.get("/api/alleles")
def api_alleles():
    avail = set(list_alleles())
    presets, missing = {}, {}
    for k, v in PRESETS.items():
        presets[k] = [a for a in v if a in avail]
        missing[k] = [a for a in v if a not in avail]
    return jsonify(all=sorted(avail), presets=presets, missing=missing)


@app.post("/api/run")
def api_run():
    d = request.get_json(force=True)
    alleles = [a.replace("*", "") for a in d.get("alleles", [])]
    if not alleles or any(not ALLELE_RE.match(a) for a in alleles):
        return jsonify(error="Select at least one valid allele."), 400
    lengths = sorted({int(x) for x in d.get("lengths", [9]) if 8 <= int(x) <= 14})
    pep_df, warns = build_peptides(d.get("text", ""), d.get("mode", "scan"), lengths)
    if pep_df.empty:
        return jsonify(error="No valid peptides.", warnings=warns), 400
    params = dict(rankS=float(d.get("rankS", 0.5)), rankW=float(d.get("rankW", 2.0)),
                  ba=bool(d.get("ba", True)), t=d.get("t"), raw=d.get("raw", ""))
    try:
        res, logs = run_netmhcpan(pep_df["Peptide"].tolist(), alleles, params)
    except (FileNotFoundError, config.ConfigError) as e:
        return jsonify(error=f"NetMHCpan is not set up: {e}. Check Settings.", warnings=warns), 500
    if res.empty:
        return jsonify(error="NetMHCpan returned no results.", log="\n".join(logs), warnings=warns), 500
    out = pep_df.drop_duplicates().merge(res, on="Peptide")
    out = out.astype(object).where(out.notna(), None)
    return jsonify(rows=out.to_dict("records"), columns=list(out.columns),
                   warnings=warns, log="\n\n".join(logs),
                   n_peptides=int(pep_df["Peptide"].nunique()))


PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>NetMHCpan BA/EL</title>
<style>
body{font:14px system-ui,sans-serif;margin:0;display:flex;height:100vh}
#side{width:430px;padding:14px;overflow:auto;background:#f4f6f8;border-right:1px solid #ccd}
#main{flex:1;padding:14px;overflow:auto}
textarea{width:100%;height:150px;font-family:monospace}
h3{margin:14px 0 6px} button{cursor:pointer}
.chip{display:inline-block;background:#dde8ff;border-radius:10px;padding:1px 8px;margin:2px;font-size:12px}
.chip b{cursor:pointer;margin-left:4px}
#alist{height:140px;overflow:auto;background:#fff;border:1px solid #ccd;padding:4px}
.presets button{margin:2px}
table{border-collapse:collapse;font-size:13px} th,td{border:1px solid #ccd;padding:2px 6px}
th{background:#eef;cursor:pointer;position:sticky;top:0} tr.SB{background:#d8f5d8} tr.WB{background:#fff6d0}
label.i{display:block;margin:3px 0} input[type=number]{width:70px}
pre{background:#222;color:#ddd;padding:8px;overflow:auto;max-height:300px}
</style></head><body>
<div id="side">
<h2>NetMHCpan 4.2 - MHC-I BA / EL</h2>
<h3>1. Sequence / peptides</h3>
<label class="i"><input type="radio" name="mode" value="scan" checked> Scan FASTA/protein into k-mers</label>
<label class="i"><input type="radio" name="mode" value="peptides"> Use peptides as given (one per line)</label>
<textarea id="text">&gt;demo
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF</textarea>
<div>Lengths: <span id="lens"></span></div>

<h3>2. Alleles</h3>
<div class="presets" id="presets"></div>
<input id="asearch" placeholder="search all alleles..." style="width:100%;margin-top:4px">
<div id="alist"></div>
<div id="chips"></div>
<button onclick="setSel([])">clear</button>

<h3>3. NetMHCpan options</h3>
<label class="i"><input type="checkbox" id="ba" checked> Binding affinity (-BA)</label>
<label class="i">Strong rank (-rankS) <input type="number" id="rankS" value="0.5" step="0.1"></label>
<label class="i">Weak rank (-rankW) <input type="number" id="rankW" value="2" step="0.1"></label>
<label class="i">Output threshold %rank (-t) <input type="number" id="t" step="0.1" placeholder="all"></label>
<label class="i">Raw extra args <input id="raw" style="width:250px"></label>

<p><button id="go" onclick="run()" style="font-size:16px;padding:6px 18px">Run</button> <span id="status"></span></p>
</div>
<div id="main">
<div id="summary"></div>
<label><input type="checkbox" id="only" onchange="render()"> binders only</label>
<button onclick="csv()">Download CSV</button>
<div id="warn" style="color:#a60"></div>
<div id="tbl"></div>
<details><summary>Raw log</summary><pre id="log"></pre></details>
</div>
<script>
let ALL=[],PRE={},MISS={},SEL=new Set(),DATA=null,sortKey=null,sortDir=1;
const $=id=>document.getElementById(id);
$("lens").innerHTML=[8,9,10,11,12,13,14].map(n=>`<label><input type="checkbox" class="len" value="${n}" ${n==9?"checked":""}>${n}</label> `).join("");
fetch("api/alleles").then(r=>r.json()).then(d=>{ALL=d.all;PRE=d.presets;MISS=d.missing;
 $("presets").innerHTML=Object.keys(PRE).map(k=>`<button onclick='setSel(PRE[${JSON.stringify(k)}])' title="${MISS[k].length?"not in this NetMHCpan: "+MISS[k].join(", "):""}">${k} (${PRE[k].length})</button>`).join("");
 drawList();});
function setSel(a){SEL=new Set(a);drawList();drawChips();}
function drawList(){const q=$("asearch").value.toLowerCase();
 const m=ALL.filter(a=>a.toLowerCase().includes(q)).slice(0,300);
 $("alist").innerHTML=m.map(a=>`<label class="i"><input type="checkbox" ${SEL.has(a)?"checked":""} onchange="tog('${a}',this.checked)"> ${a}</label>`).join("")+(m.length==300?"<i>refine search for more</i>":"");}
function tog(a,on){on?SEL.add(a):SEL.delete(a);drawChips();}
function drawChips(){$("chips").innerHTML=`<div>${SEL.size} selected</div>`+[...SEL].map(a=>`<span class="chip">${a}<b onclick="tog('${a}',false);drawList()">×</b></span>`).join("");}
$("asearch").oninput=drawList;
async function run(){
 $("go").disabled=true;$("status").textContent="running...";
 const body={text:$("text").value,mode:document.querySelector("input[name=mode]:checked").value,
  lengths:[...document.querySelectorAll(".len:checked")].map(x=>+x.value),alleles:[...SEL],
  ba:$("ba").checked,rankS:$("rankS").value,rankW:$("rankW").value,t:$("t").value,raw:$("raw").value};
 const r=await fetch("api/run",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(body)});
 const d=await r.json();$("go").disabled=false;$("status").textContent="";
 $("log").textContent=d.log||"";$("warn").textContent=(d.warnings||[]).join(" | ");
 if(d.error){$("summary").innerHTML="<b style='color:#b00'>"+d.error+"</b>";return;}
 DATA=d;render();}
function render(){if(!DATA)return;let rows=DATA.rows.slice();
 if($("only").checked)rows=rows.filter(r=>r.NMP_BindLevel);
 if(sortKey)rows.sort((a,b)=>(a[sortKey]>b[sortKey]?1:a[sortKey]<b[sortKey]?-1:0)*sortDir);
 const sb=DATA.rows.filter(r=>(r.NMP_BindLevel||"").includes("SB")).length, wb=DATA.rows.filter(r=>(r.NMP_BindLevel||"").includes("WB")).length;
 $("summary").innerHTML=`<b>${DATA.n_peptides}</b> unique peptides, ${DATA.rows.length} peptide-allele rows, ${sb} SB / ${wb} WB`;
 $("tbl").innerHTML="<table><tr>"+DATA.columns.map(c=>`<th onclick="srt('${c}')">${c}</th>`).join("")+"</tr>"+
  rows.map(r=>`<tr class="${(r.NMP_BindLevel||"").trim()}">`+DATA.columns.map(c=>`<td>${r[c]??""}</td>`).join("")+"</tr>").join("")+"</table>";}
function srt(k){sortDir=(sortKey==k)?-sortDir:1;sortKey=k;render();}
function csv(){if(!DATA)return;const l=[DATA.columns.join(",")].concat(DATA.rows.map(r=>DATA.columns.map(c=>r[c]??"").join(",")));
 const a=document.createElement("a");a.href=URL.createObjectURL(new Blob([l.join("\n")]));a.download="netmhcpan_results.csv";a.click();}
</script></body></html>"""


@app.get("/")
def index():
    return PAGE


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(os.environ.get("PORT", 5000)), threaded=True)
