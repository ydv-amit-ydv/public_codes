"""Flask runner: paste a construct, watch the cycles, download the results.

    python -m circopt.app              # http://127.0.0.1:5000
    python -m circopt.app --port 8080 --host 0.0.0.0

One run at a time per process; each run gets its own directory under `runs/`.
"""
from __future__ import annotations

import argparse
import json
import os
import threading
import traceback
import uuid
from datetime import datetime

from flask import Flask, abort, jsonify, redirect, render_template_string, request, send_file, url_for

from .config import DEFAULT_WEIGHTS, Params
from .pipeline import run_design

app = Flask(__name__)
RUNS: dict[str, dict] = {}
RUNS_DIR = os.environ.get("CIRCOPT_RUNS", os.path.abspath("runs"))

PAGE = """<!doctype html><html><head><meta charset="utf-8"><title>circRNA CDS optimiser</title>
<style>
 body{font:14px/1.5 system-ui,sans-serif;margin:0;background:#f6f7f9;color:#17202a}
 .wrap{max-width:1050px;margin:0 auto;padding:24px}
 h1{font-size:20px;margin:0 0 4px} .sub{color:#5b6b7b;margin-bottom:18px}
 .card{background:#fff;border:1px solid #dde3ea;border-radius:10px;padding:16px;margin-bottom:16px}
 label{display:block;font-weight:600;margin:10px 0 4px;font-size:13px}
 textarea,input,select{width:100%;box-sizing:border-box;font:12px ui-monospace,monospace;
   padding:7px;border:1px solid #cfd8e3;border-radius:6px}
 textarea{min-height:84px}
 .grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(190px,1fr));gap:10px}
 button{background:#1f6feb;color:#fff;border:0;border-radius:7px;padding:10px 18px;
   font-weight:600;cursor:pointer;font-size:14px}
 pre{background:#0d1117;color:#c9d1d9;padding:12px;border-radius:8px;overflow:auto;max-height:420px;
   font-size:12px;white-space:pre-wrap}
 table{border-collapse:collapse;width:100%;font-size:12px}
 th,td{border-bottom:1px solid #e6ebf1;padding:5px 7px;text-align:right}
 th:first-child,td:first-child{text-align:left}
 .note{background:#fff8e1;border:1px solid #f0d78c;padding:10px;border-radius:7px;font-size:12.5px}
 .err{background:#fff0f0;border:1px solid #f0b0b0}
 img{max-width:100%;border:1px solid #dde3ea;border-radius:8px;margin:8px 0}
 a{color:#1f6feb}
</style></head><body><div class="wrap">{{ body|safe }}</div></body></html>"""

FORM = """
<h1>circRNA vaccine CDS optimiser</h1>
<div class="sub">Paste a protein, DNA or RNA sequence. Everything else has a sensible default.</div>

<form method="post" action="{{ url_for('start') }}">
<div class="card">
 <label>Sequence</label>
 <textarea name="sequence" style="min-height:130px" placeholder="Paste one of:
  - a protein (MVSKGEELFTGVVPILVELDGDVNGHKF...)  -> the CDS is built for you
  - the CDS payload as DNA or RNA (ATGGTGAGCAAGGGC... / AUGGUGAGCAAGGGC...)
  - the whole circular molecule, IRES and all  -> the longest ORF is optimised, the rest held fixed
FASTA headers, spaces and line numbers are ignored.">{{ p.sequence }}</textarea>
 <div class="grid" style="margin-top:8px">
  <div><label>It is</label><select name="sequence_type">
    <option value="auto" {{ 'selected' if p.sequence_type=='auto' }}>detect automatically</option>
    <option value="protein" {{ 'selected' if p.sequence_type=='protein' }}>protein</option>
    <option value="dna" {{ 'selected' if p.sequence_type=='dna' }}>DNA</option>
    <option value="rna" {{ 'selected' if p.sequence_type=='rna' }}>RNA</option></select></div>
  <div><label>Which part</label><select name="sequence_role">
    <option value="auto" {{ 'selected' if p.sequence_role=='auto' }}>detect automatically</option>
    <option value="payload" {{ 'selected' if p.sequence_role=='payload' }}>CDS payload only</option>
    <option value="molecule" {{ 'selected' if p.sequence_role=='molecule' }}>whole circular molecule</option>
    </select></div>
  <div><label>Name</label><input name="name" value="{{ p.name }}"></div>
 </div>
</div>

<div class="card">
 <div class="grid">
  <div><label>Translate best in</label><select name="primary_species">
    <option {{ 'selected' if p.primary_species=='mouse' }}>mouse</option>
    <option {{ 'selected' if p.primary_species=='human' }}>human</option></select></div>
  <div><label>...while still good in</label><select name="secondary_species">
    <option {{ 'selected' if p.secondary_species=='human' }}>human</option>
    <option {{ 'selected' if p.secondary_species=='mouse' }}>mouse</option></select></div>
  <div><label>Time budget (minutes)</label><input name="time_budget_min" value="{{ p.time_budget_min }}"></div>
 </div>
 <div style="margin-top:10px">
  <input type="checkbox" name="auto" {{ 'checked' if p.auto }}> <b>Auto-fill every other parameter</b>
  <div class="sub" style="font-size:12px">Sets the GC band to what your sequence can actually
  reach at the target CAI, picks the structure thresholds, and sizes the search to the budget
  above. Anything you change by hand in Advanced is respected.</div>
 </div>
</div>

<details class="card"><summary><b>Optional: IRES and flanking elements</b></summary>
 <div class="sub" style="font-size:12px;margin:6px 0">Only used when you give the CDS payload.
 These are placed around it and never edited. If you paste the whole molecule instead, leave
 them empty &mdash; everything outside the CDS is held fixed automatically.</div>
 <label>5' scaffold (IRES, UTR, spacer)</label><textarea name="scaffold_5">{{ p.scaffold_5 }}</textarea>
 <label>3' scaffold</label><textarea name="scaffold_3">{{ p.scaffold_3 }}</textarea>
 <div style="margin-top:8px"><input type="checkbox" name="add_stop" {{ 'checked' if p.add_stop }}>
  add a stop codon if the payload has none (uncheck for an infinite-ORF design)</div>
</details>

<details class="card"><summary><b>Advanced</b> &mdash; leave alone unless you know you need it</summary>
 <div class="sub" style="font-size:12px;margin:6px 0">With auto-fill on, any field left at its
 default here is chosen for you; edit one and your value wins.</div>

 <h4>Codon usage</h4>
 <div class="grid">
  <div><label>Primary weight</label><input name="primary_weight" value="{{ p.primary_weight }}"></div>
  <div><label>CAI floor (hard, weaker species)</label><input name="cai_floor" value="{{ p.cai_floor }}"></div>
  <div><label>CAI target (term saturates)</label><input name="cai_target" value="{{ p.cai_target }}"></div>
  <div><label>Metric</label><select name="codon_metric">
    <option {{ 'selected' if p.codon_metric=='cai' }}>cai</option>
    <option {{ 'selected' if p.codon_metric=='tai' }}>tai</option></select></div>
  <div><label>tRNA copy file (for tAI)</label><input name="trna_copy_file" value="{{ p.trna_copy_file or '' }}"></div>
  <div><label>Min codon frequency</label><input name="min_codon_freq" value="{{ p.min_codon_freq }}"></div>
 </div>

 <h4>Structure</h4>
 <div class="grid">
  <div><label>Goal</label><select name="structure_goal">
    <option value="minimize_pairing" {{ 'selected' if p.structure_goal=='minimize_pairing' }}>minimise pairing</option>
    <option value="minimize_mfe" {{ 'selected' if p.structure_goal=='minimize_mfe' }}>minimise MFE (more structure)</option>
    <option value="target" {{ 'selected' if p.structure_goal=='target' }}>hold near target</option></select></div>
  <div><label>MFE/nt target</label><input name="mfe_per_nt_target" value="{{ p.mfe_per_nt_target }}"></div>
  <div><label>dsRNA duplex threshold (bp)</label><input name="dsrna_min_bp" value="{{ p.dsrna_min_bp }}"></div>
  <div><label>Stem target (bp)</label><input name="stem_target" value="{{ p.stem_target }}"></div>
  <div><label>Stem hard cap (bp)</label><input name="stem_cap" value="{{ p.stem_cap }}"></div>
  <div><label>Min paired fraction</label><input name="min_paired_fraction" value="{{ p.min_paired_fraction }}"></div>
  <div><label>AUG min p(unpaired)</label><input name="start_unpaired_min_prob" value="{{ p.start_unpaired_min_prob }}"></div>
  <div><label>Start context mean p(unpaired)</label><input name="start_context_mean_min" value="{{ p.start_context_mean_min }}"></div>
 </div>

 <h4>Composition and forbidden motifs</h4>
 <div class="grid">
  <div><label>GC min</label><input name="gc_min" value="{{ p.gc_min }}"></div>
  <div><label>GC max</label><input name="gc_max" value="{{ p.gc_max }}"></div>
  <div><label>GC window (nt)</label><input name="gc_window" value="{{ p.gc_window }}"></div>
  <div><label>GC window min</label><input name="gc_window_min" value="{{ p.gc_window_min }}"></div>
  <div><label>GC window max</label><input name="gc_window_max" value="{{ p.gc_window_max }}"></div>
  <div><label>Max poly-U run</label><input name="max_poly_u" value="{{ p.max_poly_u }}"></div>
 </div>
 <div class="grid" style="margin-top:6px">
  <div><input type="checkbox" name="avoid_re_sites" {{ 'checked' if p.avoid_re_sites }}> restriction sites</div>
  <div><input type="checkbox" name="avoid_slippery" {{ 'checked' if p.avoid_slippery }}> ribosome slippery sites</div>
  <div><input type="checkbox" name="avoid_cryptic_splice" {{ 'checked' if p.avoid_cryptic_splice }}> cryptic splice sites</div>
  <div><input type="checkbox" name="avoid_g4" {{ 'checked' if p.avoid_g4 }}> G-quadruplexes</div>
 </div>
 <label>Extra forbidden motifs (comma separated, IUPAC ok)</label>
 <input name="extra_forbidden" value="{{ ','.join(p.extra_forbidden) }}">

 <h4>Search</h4>
 <div class="grid">
  <div><label>Cycles</label><input name="rounds" value="{{ p.rounds }}"></div>
  <div><label>Generations per cycle</label><input name="generations_per_round" value="{{ p.generations_per_round }}"></div>
  <div><label>Population</label><input name="population" value="{{ p.population }}"></div>
  <div><label>Top N reported</label><input name="top_n" value="{{ p.top_n }}"></div>
  <div><label>Polish steps</label><input name="polish_codons" value="{{ p.polish_codons }}"></div>
  <div><label>Workers (0 = all cores)</label><input name="workers" value="{{ p.workers }}"></div>
  <div><label>Seed</label><input name="seed" value="{{ p.seed }}"></div>
  <div><label>Fast-fold pair span</label><input name="inner_span" value="{{ p.inner_span }}"></div>
 </div>
 <div class="grid" style="margin-top:6px">
  <div><input type="checkbox" name="use_derna" {{ 'checked' if p.use_derna }}> seed with DERNA (slow)</div>
  <div><input type="checkbox" name="stress_test" {{ 'checked' if p.stress_test }}> noLP/noGU stress test</div>
  <div><input type="checkbox" name="linear_compare" {{ 'checked' if p.linear_compare }}> linear comparison</div>
 </div>

 <h4>Or give the molecule explicitly</h4>
 <div class="sub" style="font-size:12px">Overrides the sequence box above.</div>
 <label>Full circular RNA</label><textarea name="rna">{{ p.rna }}</textarea>
 <label>Protein (to locate the CDS)</label><textarea name="protein">{{ p.protein }}</textarea>
 <div class="grid">
  <div><label>CDS start</label><input name="cds_start" value="{{ '' if p.cds_start is none else p.cds_start }}"></div>
  <div><label>CDS end (exclusive)</label><input name="cds_end" value="{{ '' if p.cds_end is none else p.cds_end }}"></div>
 </div>
 <label>Protected regions inside the CDS (name,start,end per line)</label>
 <textarea name="protected">{{ protected_text }}</textarea>
</details>

<button type="submit">Optimise</button>
</form>
{% if runs %}<div class="card"><b>Recent runs</b><table>
<tr><th>run</th><th>name</th><th>state</th><th>started</th></tr>
{% for r in runs %}<tr><td><a href="{{ url_for('show', rid=r.id) }}">{{ r.id[:8] }}</a></td>
<td>{{ r.name }}</td><td>{{ r.state }}</td><td>{{ r.started }}</td></tr>{% endfor %}
</table></div>{% endif %}
"""

RESULT = """
<h1>{{ r.name }} <span class="sub">({{ r.state }})</span></h1>
<div class="sub">run {{ r.id[:8] }} &middot; started {{ r.started }} &middot;
<a href="{{ url_for('index') }}">new run</a></div>
{% if r.state == 'running' %}<div class="card note">Running — this page refreshes every 4 s.</div>
<meta http-equiv="refresh" content="4">{% endif %}
{% if r.error %}<div class="card err"><b>Failed</b><pre>{{ r.error }}</pre></div>{% endif %}
{% if summary %}<div class="card"><b>Original → best</b><table>
<tr><th>metric</th><th>original</th><th>best</th><th>delta</th></tr>
{% for row in summary %}<tr><td>{{ row[0] }}</td><td>{{ row[1] }}</td><td>{{ row[2] }}</td>
<td>{{ row[3] }}</td></tr>{% endfor %}</table></div>{% endif %}
{% if files %}<div class="card"><b>Downloads</b><ul>
{% for f in files %}<li><a href="{{ url_for('getfile', rid=r.id, name=f) }}">{{ f }}</a></li>{% endfor %}
</ul></div>{% endif %}
{% if plots %}<div class="card"><b>Plots</b>
{% for f in plots %}<div>{{ f }}<br><img src="{{ url_for('getfile', rid=r.id, name='plots/'+f) }}"></div>{% endfor %}
</div>{% endif %}
<div class="card"><b>Log</b><pre>{{ r.log }}</pre></div>
"""


def _parse_protected(text: str) -> list[dict]:
    out = []
    for line in (text or "").splitlines():
        line = line.strip()
        if not line:
            continue
        parts = [x.strip() for x in line.replace("\t", ",").split(",")]
        if len(parts) < 3:
            raise ValueError(f"protected region needs name,start,end -> got {line!r}")
        out.append({"name": parts[0], "start": int(parts[1]), "end": int(parts[2])})
    return out


def _form_to_params(form) -> Params:
    ints = {"cds_start", "cds_end", "gc_window", "max_poly_u", "max_homopolymer", "stem_target", "stem_cap", "dsrna_min_bp",
            "rounds", "generations_per_round", "population", "top_n", "polish_codons", "workers",
            "seed", "inner_span", "elite_k", "ir_kmer", "access_window", "access_span",
            "early_orf_nt", "start_upstream", "start_downstream"}
    floats = {"primary_weight", "cai_floor", "min_codon_freq", "gc_min", "gc_max", "gc_window_min",
              "gc_window_max", "min_paired_fraction", "start_unpaired_min_prob",
              "start_context_mean_min", "mutation_rate", "mfe_per_nt_target", "cai_target",
              "time_budget_min"}
    bools = {"avoid_re_sites", "avoid_slippery", "avoid_cryptic_splice", "avoid_g4", "use_derna",
             "stress_test", "linear_compare", "gquad", "no_lonely_pairs", "auto", "add_stop"}
    strs = {"name", "protein", "rna", "primary_species", "secondary_species", "codon_metric",
            "structure_goal", "sequence", "sequence_type", "sequence_role",
            "scaffold_5", "scaffold_3", "trna_copy_file"}
    d: dict = {}
    for k in strs:
        if form.get(k) is not None:
            d[k] = form.get(k).strip()
    for k in ints:
        if form.get(k):
            d[k] = int(float(form[k]))
    for k in floats:
        if form.get(k):
            d[k] = float(form[k])
    for k in bools:
        d[k] = form.get(k) is not None
    d["protected"] = _parse_protected(form.get("protected", ""))
    ef = (form.get("extra_forbidden") or "").strip()
    d["extra_forbidden"] = [x.strip().upper() for x in ef.split(",") if x.strip()]
    d["protein"] = "".join(d.get("protein", "").split()).upper()
    d["rna"] = "".join(d.get("rna", "").split()).upper().replace("T", "U")
    if not d.get("trna_copy_file"):
        d.pop("trna_copy_file", None)
    # a blank value means "not given", not "empty string"
    for k in ("sequence", "scaffold_5", "scaffold_3"):
        if k in d and not d[k].strip():
            d.pop(k)
    return Params.from_dict(d)


def _suggest_orf(rna: str) -> str:
    """Find the longest AUG..stop ORF in the pasted RNA, to put real numbers in the error page.

    This only ever *suggests*: silently guessing a vaccine CDS would be a bad idea, so the
    user still has to confirm it by pasting the protein or the coordinates.
    """
    from .codon import CODON_TO_AA, translate
    rna = "".join((rna or "").split()).upper().replace("T", "U")
    if len(rna) < 6 or set(rna) - set("ACGU"):
        return ""
    best = None
    for i in range(len(rna) - 2):
        if rna[i:i + 3] != "AUG":
            continue
        for j in range(i + 3, len(rna) - 2, 3):
            if CODON_TO_AA.get(rna[j:j + 3]) == "*":
                if best is None or (j + 3 - i) > (best[1] - best[0]):
                    best = (i, j + 3)
                break
    if not best or (best[1] - best[0]) < 30:
        return ""
    prot = translate(rna[best[0]:best[1]]).rstrip("*")
    return (f"The longest ORF in the RNA you pasted runs {best[0]}..{best[1]} "
            f"({len(prot)} aa). If that is your CDS, either paste its protein "
            f"<code>{prot[:40]}{'...' if len(prot) > 40 else ''}</code> or enter "
            f"cds_start={best[0]} and cds_end={best[1]}.")


def _worker(rid: str, p: Params) -> None:
    r = RUNS[rid]

    def log(msg=""):
        r["log"] += str(msg) + "\n"

    try:
        res = run_design(p, r["dir"], log=log)
        r["paths"] = res["paths"]
        r["state"] = "done"
    except Exception:
        r["error"] = traceback.format_exc()
        r["state"] = "failed"


@app.route("/")
def index():
    p = Params()
    protected_text = "\n".join(f"{r.get('name','region')},{r['start']},{r['end']}"
                               for r in p.protected)
    runs = [{"id": k, "name": v["name"], "state": v["state"], "started": v["started"]}
            for k, v in sorted(RUNS.items(), key=lambda kv: kv[1]["started"], reverse=True)][:10]
    body = render_template_string(FORM, p=p, protected_text=protected_text, runs=runs,
                                  weights=DEFAULT_WEIGHTS)
    return render_template_string(PAGE, body=body)


@app.post("/run")
def start():
    try:
        p = _form_to_params(request.form)
        # Build the Design here, not in the worker: it is what checks that the protein matches
        # the CDS, that the CDS is in frame and that the protected regions are sane. Doing it
        # up front turns those into a readable page instead of a traceback in a run log.
        import copy

        from .design import Design
        from .pipeline import resolve_inputs
        Design(resolve_inputs(copy.deepcopy(p)))   # copy: the worker re-resolves and logs it
    except Exception as e:
        hint = ""
        if "locate the CDS" in str(e) or "could not locate" in str(e):
            hint = _suggest_orf(request.form.get("rna", ""))
        body = (f"<div class='card err'><b>Bad input</b><p>{e}</p>"
                + (f"<p class='note'>{hint}</p>" if hint else "")
                + f"<p><a href='{url_for('index')}'>back to the form</a></p></div>")
        return render_template_string(PAGE, body=body), 400
    rid = uuid.uuid4().hex
    rdir = os.path.join(RUNS_DIR, rid)
    os.makedirs(rdir, exist_ok=True)
    RUNS[rid] = {"id": rid, "name": p.name, "state": "running", "log": "", "dir": rdir,
                 "started": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "error": "",
                 "paths": {}}
    threading.Thread(target=_worker, args=(rid, p), daemon=True).start()
    return redirect(url_for("show", rid=rid))


@app.route("/run/<rid>")
def show(rid):
    r = RUNS.get(rid) or abort(404)
    files, plots, summary = [], [], []
    if r["state"] == "done":
        for f in sorted(os.listdir(r["dir"])):
            if os.path.isfile(os.path.join(r["dir"], f)):
                files.append(f)
        pd = os.path.join(r["dir"], "plots")
        if os.path.isdir(pd):
            plots = sorted(x for x in os.listdir(pd) if x.endswith(".png"))
        comp = os.path.join(r["dir"], "comparison_original_vs_best.csv")
        if os.path.exists(comp):
            import csv as _csv
            with open(comp) as fh:
                summary = list(_csv.reader(fh))[1:]
    body = render_template_string(RESULT, r=r, files=files, plots=plots, summary=summary)
    return render_template_string(PAGE, body=body)


@app.route("/run/<rid>/log.json")
def logjson(rid):
    r = RUNS.get(rid) or abort(404)
    return jsonify({"state": r["state"], "log": r["log"], "error": r["error"]})


@app.route("/run/<rid>/file/<path:name>")
def getfile(rid, name):
    r = RUNS.get(rid) or abort(404)
    base = os.path.abspath(r["dir"])
    path = os.path.abspath(os.path.join(base, name))
    if not path.startswith(base + os.sep) or not os.path.isfile(path):
        abort(404)
    return send_file(path)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="circopt.app")
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=5000)
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args(argv)
    os.makedirs(RUNS_DIR, exist_ok=True)
    print(f"circopt runner on http://{a.host}:{a.port}  (runs -> {RUNS_DIR})")
    app.run(host=a.host, port=a.port, debug=a.debug, threaded=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
