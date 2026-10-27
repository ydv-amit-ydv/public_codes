"""Tests for the engine server, config handling and both model wrappers.

    python3 test_engine.py                         # stubs only
    BIGMHC_PYTHON=... BIGMHC_HOME=... python3 test_engine.py   # also run the real BigMHC

NetMHCpan and BigMHC are replaced by stub executables whose output copies the
real formats (NetMHCpan 4.2d table; BigMHC .prd with mhc,pep,len,BigMHC_* columns).
"""
import json, os, stat, sys, tempfile, textwrap

REAL_PY, REAL_HOME = os.environ.get("BIGMHC_PYTHON"), os.environ.get("BIGMHC_HOME")   # optional real BigMHC
TMP = tempfile.mkdtemp()
CFG = os.path.join(TMP, "cfg.json")
os.environ["ENGINE_CONFIG"] = CFG
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from common import alleles, config           # noqa: E402

# ---- alleles ------------------------------------------------------------------
assert alleles.canonical("A0201") == "A0201"                    # no colon: left alone, never guessed
assert alleles.canonical("HLA-A02:01") == "HLA-A*02:01"
assert alleles.canonical("a*02:01:01") == "HLA-A*02:01"
assert alleles.netmhcpan_name("HLA-A*02:01") == "HLA-A02:01"
assert alleles.netmhcpan_name("H-2-Kb") == "H-2-Kb"
assert len(alleles.PRESETS["IEDB human reference 27"]["alleles"]) == 27

# ---- config ---------------------------------------------------------------------
assert config.model("bigmhc")["device"] == "cpu"                # defaults with no file
for bad in ('{"models": {"bigmhc": {"jobs": "two"}}}', '{"models": []}', "not json", '{"models":{"x":{"enabled":1}}}'):
    try:
        config.save_text(bad); raise SystemExit(f"accepted bad config: {bad}")
    except config.ConfigError:
        pass
assert not os.path.exists(CFG), "a rejected edit must not create/overwrite the file"

# ---- stubs ------------------------------------------------------------------------
NMP = os.path.join(TMP, "nmp"); os.makedirs(NMP)
open(f"{NMP}/netMHCpan", "w").write(textwrap.dedent('''\
    #!/bin/bash
    setenv() { NMHOME_SEEN="$2"; }
    setenv NMHOME /wrong/placeholder
    if [ -n "$EXPECTED_NMHOME" ] && [ "$NMHOME_SEEN" != "$EXPECTED_NMHOME" ]; then echo "bad NMHOME $NMHOME_SEEN" >&2; exit 3; fi
    [ -n "$STUB_SLEEP" ] && [ "$1" != "-listMHC" ] && sleep "$STUB_SLEEP"
    if [ "$1" = "-listMHC" ]; then printf 'HLA-A01:01\\nHLA-A02:01\\nHLA-B07:02\\nH-2-Kb\\n'; exit 0; fi
    cat <<'EOT'
     Pos            MHC           Peptide      Core Of Gp Gl Ip Il             Icore        Identity Score_EL %Rank_EL Score_BA %Rank_BA  Aff(nM) BindLevel
       1    HLA-A*02:01          SIINFEKL SII-NFEKL  0  0  0  3  1          SIINFEKL         PEPLIST 0.016383    5.510 0.106826   17.634 15739.83
       2    HLA-A*02:01         GILGFVFTL GILGFVFTL  0  0  0  0  0         GILGFVFTL         PEPLIST 0.711032    0.197 0.781118    0.084    10.68 <= SB
    EOT
    '''))
os.chmod(f"{NMP}/netMHCpan", 0o755)

BM = os.path.join(TMP, "bm"); os.makedirs(f"{BM}/src")
open(f"{BM}/src/predict.py", "w").write(textwrap.dedent('''\
    import sys, pandas as pd
    a = dict(x.lstrip("-").split("=", 1) for x in sys.argv[1:])
    df = pd.read_csv(a["i"], header=None, names=["mhc", "pep"])
    df["len"] = df["pep"].str.len()
    df["BigMHC_" + a["m"].upper()] = df["pep"].str.len() / 20.0
    for k in range(7): print(f"loading model /x/bat{k}...done")
    print("starting prediction devices: []")
    for i in (1, 2): print(f"batch {i}/2")
    df.to_csv(a["o"], index=False)
    '''))


# stand-ins for MixMHCpred 3.0 and PRIME 2.1 that write the real file formats
# (GILGFVFTL is a strong hit for everything, any other peptide is rank 50)
MIX = os.path.join(TMP, "mix"); os.makedirs(f"{MIX}/code"); os.makedirs(f"{MIX}/lib")
open(f"{MIX}/code/main.py", "w").write(""); open(f"{MIX}/lib/alleles_list.txt", "w").write("Allele\t8\t9\nA0201\t1\t1\nB0702\t1\t1\n")
STUBPY = textwrap.dedent('''\
    import sys
    a = sys.argv[1:]; g = lambda k: a[a.index(k) + 1]
    peps = [l.strip() for l in open(g("-i")) if l.strip()]; als = g("-a").split(",")
    prime = "PRIME" in sys.argv[0]
    rank = lambda p, k: (0.1 if k == 0 else 0.2) if p == "GILGFVFTL" else 50.0
    with open(g("-o"), "w") as f:
        f.write("####################\\n# Output (stub)\\n####################\\n")
        if prime:
            f.write("Peptide\\t%Rank_bestAllele\\t" + "\\t".join(f"%Rank_{x}\\tScore_{x}\\t%RankBinding_{x}" for x in als) + "\\n")
            for p in peps: f.write(p + "\\t1\\t" + "\\t".join(f"{rank(p, 1)}\\t0.5\\t{rank(p, 0)}" for x in als) + "\\n")
        else:
            f.write("Peptide\\tScore_bestAllele\\t" + "\\t".join(f"Score_{x}\\t%Rank_{x}" for x in als) + "\\n")
            for p in peps: f.write(p + "\\t1\\t" + "\\t".join(f"0.3\\t{rank(p, 0)}" for x in als) + "\\n")
    ''')
open(f"{MIX}/MixMHCpred", "w").write("#!/bin/bash\nexec python3 - \"$@\" <<'EOF'\n" + STUBPY.replace("sys.argv[0]", "'MixMHCpred'") + "EOF\n")
PRM = os.path.join(TMP, "prm"); os.makedirs(f"{PRM}/lib")
open(f"{PRM}/lib/PRIME.cc", "w").write(""); open(f"{PRM}/lib/PRIME.x", "wb").write(b"\x7fELF stub")
open(f"{PRM}/lib/alleles.txt", "w").write("A0201\nB0702\nA0101\n")
open(f"{PRM}/PRIME", "w").write("#!/bin/bash\nexec python3 - \"$@\" <<'EOF'\n" + STUBPY.replace("sys.argv[0]", "'PRIME'") + "EOF\n")


# stand-in for MHCflurry's command: it only gives a good processing score when flanks are supplied
MFD = os.path.join(TMP, "mfd"); os.makedirs(MFD)
open(f"{MFD}/mhcflurry-predict", "w").write("#!%s\n" % sys.executable + textwrap.dedent('''\
    import sys, csv
    a = sys.argv[1:]; inp = a[0]; out = a[a.index("--out") + 1]
    rows = list(csv.DictReader(open(inp)))
    with open(out, "w") as f:
        f.write("allele,peptide,mhcflurry_affinity,mhcflurry_affinity_percentile,mhcflurry_processing_score,mhcflurry_presentation_score,mhcflurry_presentation_percentile\\n")
        for r in rows:
            hit = r["peptide"] == "GILGFVFTL"
            if r["allele"].startswith("H-2"): f.write(f"{r['allele']},{r['peptide']},,,,,\\n"); continue      # --no-throw: NaN
            flank = bool(r.get("n_flank")) or bool(r.get("c_flank"))
            f.write(f"{r['allele']},{r['peptide']},{25.0 if hit else 5000.0},{0.1 if hit else 9.0},{0.9 if flank else 0.1},{0.8 if hit else 0.01},{0.2 if hit else 40.0}\\n")
    '''))
os.chmod(f"{MFD}/mhcflurry-predict", 0o755)

# stand-in for the IEDB T Cell Class I bundle: refuses what the real one refuses
BUN = os.path.join(TMP, "bundle", "ng_tc1-v0.2-beta"); os.makedirs(f"{BUN}/src"); os.makedirs(f"{BUN}/.venv/bin")
open(f"{BUN}/pyproject.toml", "w").write("[project]\nname='ng_tc1'\n")
os.symlink(sys.executable, f"{BUN}/.venv/bin/python")
open(f"{BUN}/src/tcell_mhci.py", "w").write(textwrap.dedent('''\
    import json, sys
    a = sys.argv[1:]; j = json.load(open(a[a.index("-j") + 1])); pred = j["predictors"][0]
    seqs = [l for l in j["input_sequence_text"].splitlines()]
    if any(c not in ">0123456789ACDEFGHIKLMNPQRSTVWYX" for l in seqs for c in l): sys.exit("stub: stray character reached the bundle")
    recs = list(zip(seqs[0::2], seqs[1::2]))
    if pred["method"] == "netchop":
        print("sequence_number\\tposition\\tamino_acid\\tprediction_score")
        for i, (_, sq) in enumerate(recs, 1):
            for k, c in enumerate(sq, 1): print(f"{i}\\t{k}\\t{c}\\t{0.9 if c == 'L' else 0.1}")
    else:
        if "," in j["alleles"]: sys.exit("stub: several alleles in one job (the real bundle mislabels them)")
        if j["alleles"] == "HLA-A*99:99": sys.stderr.write("Allele(s) incompatible with netctlpan: HLA-A*99:99. try quotes.\\n"); sys.exit(1)
        lo, hi = j["peptide_length_range"]
        print("sequence_number\\tstart\\tpeptide\\tallele\\tmhc_prediction\\ttap_prediction_score\\tcleavage_prediction_score\\tcombined_prediction_score\\tpercentile_rank")
        for i, (_, sq) in enumerate(recs, 1):
            for L in range(lo, hi + 1):
                for st in range(len(sq) - L + 1):
                    pep = sq[st:st + L]; hit = pep == "GILGFVFTL"
                    print(f"{i}\\t{st + 1}\\t{pep}\\t{j['alleles']}\\t0.5\\t1.0\\t0.9\\t{0.96 if hit else 0.1}\\t{0.4 if hit else 50.0}")
    '''))

settings = {"tmp_dir": os.path.join(TMP, "work"),
            "models": {"netmhcpan": {"home": NMP}, "bigmhc": {"python": sys.executable, "home": BM},
                       "mixmhcpred": {"home": MIX, "python": sys.executable}, "prime": {"home": PRM, "python": sys.executable},
                       "mhcflurry": {"python": sys.executable, "command": f"{MFD}/mhcflurry-predict"},
                       "iedb_bundle": {"home": BUN, "python": f"{BUN}/.venv/bin/python"}}}

# ---- discovery: copy-and-run layout (no explicit paths anywhere) ----------------------
import shutil                                  # noqa: E402
LE = os.path.join(TMP, "local_engines")
inner = os.path.join(LE, "netmhc", "netMHCpan-4.2", "netMHCpan-4.2")           # the tarball nests like this
os.makedirs(os.path.join(inner, "Linux_x86_64")); os.makedirs(os.path.join(inner, "data"))
shutil.copy(f"{NMP}/netMHCpan", f"{inner}/netMHCpan"); os.chmod(f"{inner}/netMHCpan", 0o755)
bm_home = os.path.join(LE, "bigmhc", "bigmhc"); os.makedirs(f"{bm_home}/models"); os.makedirs(f"{bm_home}/src")
shutil.copy(f"{BM}/src/predict.py", f"{bm_home}/src/predict.py")
os.makedirs(os.path.join(LE, "netmhc", "decoy", "netMHCpan"))                     # a *directory* named netMHCpan is not a script
shutil.copytree(MIX, os.path.join(LE, "mixmhcpred", "MixMHCpred")); shutil.copytree(PRM, os.path.join(LE, "PRIME"))
os.environ["ENGINE_MODELS_DIR"] = LE
os.environ["BIGMHC_PYTHON"] = sys.executable
os.environ["EXPECTED_NMHOME"] = inner
from common import discover                    # noqa: E402
assert discover.find_netmhcpan(LE) == inner and discover.find_bigmhc(LE) == bm_home
assert discover.find_netmhcpan(os.path.join(TMP, "empty-nowhere")) is None
assert discover.find_mixmhcpred(LE) == os.path.join(LE, "mixmhcpred", "MixMHCpred") and discover.find_prime(LE) == os.path.join(LE, "PRIME")
assert config.model("mixmhcpred")["python"] and config.model("prime")["home"].endswith("PRIME")
m = config.model("netmhcpan"); assert m["home"] == inner and "home" in m["auto"], m
run = discover.runnable_netmhcpan(inner, TMP)
assert run != f"{inner}/netMHCpan" and f"setenv NMHOME {inner}" in open(run).read()   # patched copy, original untouched
assert "/wrong/placeholder" in open(f"{inner}/netMHCpan").read()
fixed = os.path.join(TMP, "fixed"); os.makedirs(fixed); shutil.copy(f"{inner}/netMHCpan", fixed)
t = open(f"{fixed}/netMHCpan").read().replace("/wrong/placeholder", fixed); open(f"{fixed}/netMHCpan", "w").write(t)
assert discover.runnable_netmhcpan(fixed, TMP) == f"{fixed}/netMHCpan"             # already correct: used as is
assert config.model("bigmhc")["home"] == bm_home and config.model("bigmhc")["python"] == sys.executable
assert config.tmp_dir().endswith("epitope_engine")

import engine_server                            # noqa: E402
from werkzeug.test import Client                # noqa: E402
auto = Client(engine_server.application)
assert not os.path.exists(CFG)                  # nothing configured at all
r = auto.get("/netmhcpan/api/alleles").get_json(); assert "HLA-A02:01" in r["all"], r
r = auto.post("/netmhcpan/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"]}).get_json()
assert r["rows"][0]["NMP_BindLevel"] == "SB", r                                      # ran through the patched copy
r = auto.post("/bigmhc/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"], "models": ["el"]}).get_json()
assert r["rows"][0]["BigMHC_EL"] > 0, r
eff = auto.get("/api/config").get_json()["effective"]
assert eff["models"]["netmhcpan"]["home"] == inner and eff["models"]["netmhcpan"]["auto"] == ["home"]
# nothing to find -> clear messages, not crashes
os.environ["ENGINE_MODELS_DIR"] = os.path.join(TMP, "nothing-here"); os.makedirs(os.environ["ENGINE_MODELS_DIR"])
del os.environ["BIGMHC_PYTHON"]
engine_server.netmhcpan_app._list_alleles.cache_clear()
r = auto.post("/netmhcpan/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"]})
assert r.status_code == 500 and "Settings" in r.get_json()["error"]
h = auto.post("/api/health/netmhcpan").get_json(); assert not h["ok"] and "not found under" in h["checks"][0]["detail"]
r = auto.post("/bigmhc/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"], "models": ["el"]})
assert r.status_code == 500 and "not found" in r.get_json()["error"] or "No Python" in r.get_json()["error"]
for k in ("ENGINE_MODELS_DIR", "EXPECTED_NMHOME"): os.environ.pop(k, None)
engine_server.netmhcpan_app._list_alleles.cache_clear()

# ---- server -----------------------------------------------------------------------
c = Client(engine_server.application)

r = c.get("/api/config").get_json()
assert r["path"] == CFG and "bigmhc" in r["text"]                # defaults shown before any save

bad = c.put("/api/config", json={"text": '{"models": 5}'})
assert bad.status_code == 400 and not os.path.exists(CFG)
ok = c.put("/api/config", json={"text": json.dumps(settings)})
assert ok.status_code == 200 and os.path.exists(CFG)

# NetMHCpan through the mount point
a = c.get("/netmhcpan/api/alleles").get_json()
assert "HLA-A02:01" in a["all"]
p = a["presets"]["IEDB human reference 27"]
assert p == ["HLA-A01:01", "HLA-A02:01", "HLA-B07:02"] and "HLA-B58:01" in a["missing"]["IEDB human reference 27"]
r = c.post("/netmhcpan/api/run", json={"text": "GILGFVFTL SIINFEKL", "mode": "peptides",
                                       "alleles": ["HLA-A02:01"], "rankS": 0.5, "rankW": 2, "ba": True}).get_json()
by = {x["Peptide"]: x for x in r["rows"]}
assert by["GILGFVFTL"]["NMP_BindLevel"] == "SB" and by["SIINFEKL"]["NMP_BindLevel"] == ""
assert by["GILGFVFTL"]["NMP_Aff_nM"] == 10.68

# a bad NetMHCpan path in Settings must degrade gracefully, not 500
nope = json.loads(json.dumps(settings)); nope["models"]["netmhcpan"]["home"] = "/nonexistent"
c.put("/api/config", json={"text": json.dumps(nope)})
r = c.get("/netmhcpan/api/alleles"); assert r.status_code == 200 and r.get_json()["all"] == []
r = c.post("/netmhcpan/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"]})
assert r.status_code == 500 and "Settings" in r.get_json()["error"]
c.put("/api/config", json={"text": json.dumps(settings)})
netmhcpan_app = engine_server.netmhcpan_app; netmhcpan_app._list_alleles.cache_clear()

# BigMHC through the mount point, using a different interpreter path than the server's own
assert c.get("/bigmhc/api/presets").get_json()["IEDB human reference 27"][0] == "HLA-A*01:01"
assert "Mouse H-2 (all)" not in c.get("/bigmhc/api/presets").get_json()
r = c.post("/bigmhc/api/run", json={"text": ">p\nGILGFVFTLN", "mode": "scan", "lengths": [9],
                                    "alleles": ["HLA-A02:01"], "models": ["el", "im"]}).get_json()
assert len(r["rows"]) == 2 and {"BigMHC_EL", "BigMHC_IM"} <= set(r["columns"]), r
assert r["rows"][0]["allele"] == "HLA-A*02:01"

# health checks
h = c.post("/api/health/netmhcpan").get_json(); assert h["ok"], h
h = c.post("/api/health/netmhcpan?deep=1").get_json(); assert h["ok"], h
h = c.post("/api/health/bigmhc").get_json()
assert not h["ok"] and any("imports" in k["name"] for k in h["checks"])   # stub env has no torch: must say so

# edit on the fly: break the path, the very next request reports it; fix it, it works again
broken = json.loads(json.dumps(settings)); broken["models"]["bigmhc"]["home"] = "/nonexistent"
c.put("/api/config", json={"text": json.dumps(broken)})
r = c.post("/bigmhc/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"], "models": ["el"]})
assert r.status_code == 500 and "not found" in r.get_json()["error"]
c.put("/api/config", json={"text": json.dumps(settings)})
r = c.post("/bigmhc/api/run", json={"text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A02:01"], "models": ["el"]})
assert r.status_code == 200

# disabled model
off = json.loads(json.dumps(settings)); off["models"]["bigmhc"]["enabled"] = False
c.put("/api/config", json={"text": json.dumps(off)})
assert c.post("/api/health/bigmhc").get_json()["ok"] is False
c.put("/api/config", json={"text": json.dumps(settings)})

# hand-corrupted file is reported, not crashed on
open(CFG, "w").write("{oops")
g = c.get("/api/config").get_json(); assert g["error"] and "oops" in g["text"]
open(CFG, "w").write(json.dumps(settings))


# ---- jobs / pipeline ---------------------------------------------------------------------
import time                                      # noqa: E402
from common.jobs import Job                      # noqa: E402
import bigmhc_runner                             # noqa: E402
import pandas as pd                              # noqa: E402

# BigMHC progress lines drive the stage (loading 7 networks = 0-40%, batches = 40-100%)
j = Job(); j.stage("bigmhc_el", "BigMHC EL")
seen = []
orig = j.update
j.update = lambda *a, **k: (seen.append((k.get("done"), k.get("total"))), orig(*a, **k))[1]
bigmhc_runner.predict_pairs(pd.DataFrame({"peptide": ["GILGFVFTL"], "allele": ["HLA-A02:01"]}), ["el"], job=j)
dones = [d_ for d_, _ in seen if d_ is not None]
assert dones == sorted(dones) and 100 in [t for _, t in seen], seen          # one 0-100 scale, never backwards
assert j.stages[0]["state"] == "done" and j.stages[0]["total"] == 100, j.stages
assert 40 in dones, dones                                                    # 7 networks loaded = 40%


# the pipeline preflights BigMHC's imports, so give the stub env fake torch/psutil modules
FAKE = os.path.join(TMP, "fake"); os.makedirs(FAKE)
open(f"{FAKE}/torch.py", "w").write("__version__ = 'stub'\n"); open(f"{FAKE}/psutil.py", "w").write("")
os.makedirs(f"{FAKE}/mhcflurry"); open(f"{FAKE}/mhcflurry/__init__.py", "w").write("__version__ = 'stub'\n")
open(f"{FAKE}/mhcflurry/downloads.py", "w").write(f"def get_path(name, test_exists=True):\n    return {MFD!r}\n")
os.environ["PYTHONPATH"] = FAKE


def wait_job(jid, timeout=60):
    t0 = time.time()
    while time.time() - t0 < timeout:
        d = c.get(f"/api/jobs/{jid}").get_json()
        if d["state"] in ("done", "failed", "cancelled"):
            return d
        time.sleep(0.2)
    raise SystemExit("job timed out")


params = {"text": ">p\nGILGFVFTLN", "mode": "scan", "lengths": [9], "alleles": ["HLA-A*02:01", "H-2-Kb"],
          "models": {"netmhcpan": True, "bigmhc_el": True, "bigmhc_im": True}, "chain": True,
          "filters": {"rank_el_max": 2, "bm_im_min": None}}
pv = c.post("/api/pipeline/preview", json=params).get_json()
assert pv["n_peptides"] == 2 and pv["chained"] and any("BigMHC is human" in w for w in pv["warnings"]), pv
assert c.post("/api/pipeline/run", json={**params, "models": {}}).status_code == 400          # nothing selected
assert c.post("/api/pipeline/run", json={**params, "alleles": []}).status_code == 400

r = c.post("/api/pipeline/run", json=params); assert r.status_code == 200, r.get_json()
d = wait_job(r.get_json()["id"])
assert d["state"] == "done" and d["progress"] == 1.0, d
assert [x["name"] for x in d["stages"]] == ["netmhcpan", "bigmhc_el", "bigmhc_im", "merge"]
assert all(x["state"] == "done" for x in d["stages"]), d["stages"]
sm = d["summary"]
assert sm["chained"] and sm["pairs"] == 1 and sm["pairs_passed"] == 1 and sm["peptides_passed"] == 1, sm
res = c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()
row = res["rows"][0]
assert row["Peptide"] == "GILGFVFTL" and row["Allele"] == "HLA-A*02:01" and row["Pass_All"] is True
assert row["NMP_%Rank_EL"] == 0.197 and row["BigMHC_EL"] == 0.45 and row["BigMHC_IM"] == 0.45, row
summ = c.get(f"/api/jobs/{d['id']}/result?table=summary").get_json()["rows"][0]
assert summ["Alleles_passed"] == 1 and summ["Alleles_passed_list"] == "HLA-A*02:01", summ
assert c.get(f"/api/jobs/{d['id']}/csv?table=pairs").data.startswith(b"Protein,Start,Peptide")

# a filter nothing passes: BigMHC is told there is nothing to score, the run still succeeds
r = c.post("/api/pipeline/run", json={**params, "filters": {"rank_el_max": 0.01}})
d = wait_job(r.get_json()["id"]); assert d["state"] == "done" and d["summary"]["pairs_passed"] == 0, d
assert any("nothing to score" in x["msg"] for x in d["stages"]), d["stages"]

# unchained: BigMHC scores every human pair
r = c.post("/api/pipeline/run", json={**params, "chain": False, "alleles": ["HLA-A*02:01", "HLA-B*07:02"]})
d = wait_job(r.get_json()["id"]); assert d["state"] == "done" and not d["summary"]["chained"], d

# cancel kills the running NetMHCpan process promptly
os.environ["STUB_SLEEP"] = "40"
r = c.post("/api/pipeline/run", json={**params, "models": {"netmhcpan": True}})
jid = r.get_json()["id"]; time.sleep(1.5)
t0 = time.time(); c.post(f"/api/jobs/{jid}/cancel")
d = wait_job(jid, 15); assert d["state"] == "cancelled" and time.time() - t0 < 10, d
del os.environ["STUB_SLEEP"]

# an unready model is refused up front with the reason, not after a long wait
broken = json.loads(json.dumps(settings)); broken["models"]["bigmhc"]["home"] = "/nonexistent"
c.put("/api/config", json={"text": json.dumps(broken)})
r = c.post("/api/pipeline/run", json=params)
assert r.status_code == 400 and "bigmhc is not ready" in r.get_json()["error"], r.get_json()
c.put("/api/config", json={"text": json.dumps(settings)})
assert c.get("/api/jobs/doesnotexist").status_code == 404
page = c.get("/").get_data(as_text=True)
assert "Epitope Filter Engine" in page and 'fetch("api/alleles")' in page
al = c.get("/api/alleles").get_json()
assert "HLA-A*02:01" in al["all"] and al["presets"]["Mouse H-2 (all)"]["human"] is False


# ---- MixMHCpred + PRIME in the pipeline -----------------------------------------------------
full = {"text": ">p\nGILGFVFTLN", "mode": "scan", "lengths": [9], "alleles": ["HLA-A*02:01", "H-2-Kb"],
        "models": {"netmhcpan": True, "mixmhcpred": True, "prime": True, "bigmhc_el": True, "bigmhc_im": True},
        "chain": True, "filters": {"rank_el_max": 2, "mix_rank_max": 1, "prime_rank_max": 1}}
pv = c.post("/api/pipeline/preview", json=full).get_json()
assert any("MixMHCpred does not support: H-2-Kb" in w for w in pv["warnings"]) and pv["n_pairs_fast"] == 4, pv
r = c.post("/api/pipeline/run", json=full); assert r.status_code == 200, r.get_json()
d = wait_job(r.get_json()["id"])
assert d["state"] == "done", d
assert [x["name"] for x in d["stages"]] == ["netmhcpan", "mixmhcpred", "prime", "bigmhc_el", "bigmhc_im", "merge"], d["stages"]
row = c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"][0]
assert row["MixMHC_%Rank"] == 0.1 and row["PRIME_%Rank"] == 0.2 and row["PRIME_%RankBinding"] == 0.1, row
assert row["Pass_NMP"] and row["Pass_Mix"] and row["Pass_PRIME"] and row["Pass_BigMHC"] and row["Pass_All"], row
assert d["summary"]["pairs_passed"] == 1 and d["summary"]["chained"], d["summary"]
srow = c.get(f"/api/jobs/{d['id']}/result?table=summary").get_json()["rows"][0]
assert srow["Best_MixMHC_%Rank"] == 0.1 and srow["Best_PRIME_%Rank"] == 0.2, srow

# the mix filter gates BigMHC: a peptide failing it is never sent to BigMHC
strict = {**full, "text": "GILGFVFTL\nILGFVFTLN", "mode": "peptides", "alleles": ["HLA-A*02:01"], "filters": {"mix_rank_max": 1}}
d = wait_job(c.post("/api/pipeline/run", json={**strict, "models": {"mixmhcpred": True, "bigmhc_el": True}}).get_json()["id"])
rows = {r["Peptide"]: r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
assert rows["GILGFVFTL"]["BigMHC_EL"] is not None and rows["ILGFVFTLN"]["BigMHC_EL"] is None, rows
assert rows["GILGFVFTL"]["Pass_All"] and not rows["ILGFVFTLN"]["Pass_All"] and d["summary"]["chained"]

# PRIME alone: all selected alleles it knows get scored; one it does not know cannot pass
d = wait_job(c.post("/api/pipeline/run", json={**full, "models": {"prime": True}, "alleles": ["HLA-A*02:01", "HLA-B*07:02"],
                                              "filters": {"prime_rank_max": 1}}).get_json()["id"])
assert d["state"] == "done" and d["summary"]["pairs"] == 4 and d["summary"]["pairs_passed"] == 2, d["summary"]
# every selected allele unsupported by the tool: nothing scored, nothing passes, and it says why
d = wait_job(c.post("/api/pipeline/run", json={**full, "models": {"mixmhcpred": True}, "alleles": ["H-2-Kb"]}).get_json()["id"])
assert d["state"] == "done" and d["summary"]["pairs_passed"] == 0 and any("does not support" in w for w in d["summary"]["warnings"]), d
# a filter for a model that was not run is reported, not silently dropped
d = wait_job(c.post("/api/pipeline/run", json={**full, "models": {"prime": True}, "filters": {"mix_rank_max": 1}}).get_json()["id"])
assert any("Filters ignored" in w and "mix_rank_max" in w for w in d["summary"]["warnings"]), d["summary"]
# PRIME without MixMHCpred installed is refused up front
nomix = json.loads(json.dumps(settings)); nomix["models"]["mixmhcpred"]["home"] = "/nonexistent"
c.put("/api/config", json={"text": json.dumps(nomix)})
r = c.post("/api/pipeline/run", json={**full, "models": {"prime": True}})
assert r.status_code == 400 and "mixmhcpred is not ready" in r.get_json()["error"], r.get_json()
c.put("/api/config", json={"text": json.dumps(settings)})
for name in ("mixmhcpred", "prime"):
    h = c.post(f"/api/health/{name}?deep=1").get_json(); assert h["ok"], h


# ---- IEDB Class I immunogenicity (Calis 2013), native implementation --------------------------------
import calis_immunogenicity as ci                 # noqa: E402
# scores IEDB's website produced for HEV ORF1 peptides (default mask), lengths 8-14, and scores from
# IEDB's own immunogenicity_predictor.py (next-generation T Cell Class I v0.2-beta) for the other masks
CALIS = {"web_default": {"AARWFEAN": 0.48096, "LAARWFEA": 0.47518, "LAARWFEAN": 0.49804, "FGPWFRAIE": 0.49279, "FGPWFRAIEK": 0.59134, "ELAARWFEAN": 0.53573, "LFGPWFRAIEK": 0.58759, "TGNFWFHPEGL": 0.57076, "ALFGPWFRAIEK": 0.64905, "CALFGPWFRAIE": 0.6256, "CALFGPWFRAIEK": 0.72415, "RNAWRRRGFAAFT": 0.67204, "RNAWRRRGFAAFTP": 0.7302, "FCALFGPWFRAIEK": 0.72549}, "peps": ["AARWFEAN", "LAARWFEA", "LAARWFEAN", "FGPWFRAIE", "FGPWFRAIEK", "ELAARWFEAN"], "custom_2_5_9": [0.36318, 0.28234, 0.28234, 0.37879, 0.41884, 0.46247], "by_allele_A2402": [0.44416, 0.41354, 0.41354, 0.45977, 0.49601, 0.40267], "by_allele_B0801": [0.36318, 0.28234, 0.28234, 0.37879, 0.41884, 0.46247]}
for pep, want in CALIS["web_default"].items():
    assert abs(ci.score_peptides([pep])[0] - want) < 1e-9, (pep, want)
P6 = CALIS["peps"]
got, _ = ci.predict(P6, ["HLA-A*02:01"], "custom", "2,5,9"); assert list(got["IEDB_Immunogenicity"]) == [x for _, x in sorted(zip(P6, CALIS["custom_2_5_9"]))] or \
    dict(zip(got.Peptide, got.IEDB_Immunogenicity)) == dict(zip(P6, CALIS["custom_2_5_9"]))
for a, key in (("HLA-A*24:02", "by_allele_A2402"), ("HLA-B*08:01", "by_allele_B0801")):
    got, notes = ci.predict(P6, [a], "by_allele"); assert not notes
    assert dict(zip(got.Peptide, got.IEDB_Immunogenicity)) == dict(zip(P6, CALIS[key])), a
got, notes = ci.predict(P6, ["HLA-A*99:99"], "by_allele")
assert dict(zip(got.Peptide, got.IEDB_Immunogenicity)) == CALIS["web_default"] | {} or notes   # falls back to default, says so
assert notes and "default mask" in notes[0]
for bad in ("", "0", "a,b", "-2"):
    try:
        ci.parse_positions(bad); raise SystemExit(f"accepted bad positions {bad!r}")
    except ValueError:
        pass

iedb = {"text": "RNAWRRRGFAAFTP\nFGPWFRAIEK", "mode": "peptides", "alleles": ["HLA-A*02:01", "HLA-B*07:02"],
        "models": {"iedb_imm": True}, "filters": {"iedb_imm_min": 0.6}, "iedb": {"mask_choice": "default"}}
d = wait_job(c.post("/api/pipeline/run", json=iedb).get_json()["id"])
rows = c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]
by = {(r["Peptide"], r["Allele"]): r for r in rows}
assert by[("RNAWRRRGFAAFTP", "HLA-A*02:01")]["IEDB_Immunogenicity"] == 0.7302 and by[("RNAWRRRGFAAFTP", "HLA-A*02:01")]["Pass_All"]
assert by[("FGPWFRAIEK", "HLA-B*07:02")]["IEDB_Immunogenicity"] == 0.59134 and not by[("FGPWFRAIEK", "HLA-B*07:02")]["Pass_All"]   # 0.59 < 0.6
assert d["summary"]["pairs"] == 4 and d["summary"]["pairs_passed"] == 2, d["summary"]
srow = c.get(f"/api/jobs/{d['id']}/result?table=summary").get_json()["rows"][0]
assert srow["Max_IEDB_Immunogenicity"] == 0.7302, srow
# allele-specific masks differ between alleles
d = wait_job(c.post("/api/pipeline/run", json={**iedb, "filters": {}, "iedb": {"mask_choice": "by_allele"},
                                              "alleles": ["HLA-A*02:01", "HLA-A*24:02"]}).get_json()["id"])
sc = {(r["Peptide"], r["Allele"]): r["IEDB_Immunogenicity"] for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
assert sc[("FGPWFRAIEK", "HLA-A*02:01")] != sc[("FGPWFRAIEK", "HLA-A*24:02")], sc
# bad custom positions are refused with a readable message
r = c.post("/api/pipeline/run", json={**iedb, "iedb": {"mask_choice": "custom", "positions": "x,y"}})
assert r.status_code == 400 and "Positions to mask" in r.get_json()["error"], r.get_json()
r = c.post("/api/pipeline/run", json={**iedb, "iedb": {"mask_choice": "bogus"}}); assert r.status_code == 400


# ---- MHCflurry, NetCTLpan, NetChop ---------------------------------------------------------------------
import shutil as _sh                                                       # noqa: E402
SEQ = "MSLLTEVETYVLSIIPSGPLKAEIAQRLEDVFAGKNTDLEVLMEWLKTRPILSPLTKGILGFVFTLTVPSERG"
proc = {"text": ">flu desc\n" + SEQ, "mode": "scan", "lengths": [9], "alleles": ["HLA-A*02:01", "HLA-B*07:02", "H-2-Kb"],
        "chain": False, "filters": {}}
# MHCflurry: flanks are passed in scan mode (processing 0.9) and absent for loose peptides (0.1); unsupported allele reported
d = wait_job(c.post("/api/pipeline/run", json={**proc, "models": {"mhcflurry": True}}).get_json()["id"])
assert d["state"] == "done", d
rows = {(r["Peptide"], r["Allele"]): r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"] if r["Allele"] != "H-2-Kb"}
g = rows[("GILGFVFTL", "HLA-A*02:01")]
assert g["MHCF_Affinity_nM"] == 25.0 and g["MHCF_Processing"] == 0.9 and g["MHCF_Presentation_%Rank"] == 0.2 and g["Pass_MHCF"], g
assert any("MHCflurry does not support: H-2-Kb" in w for w in d["summary"]["warnings"]), d["summary"]["warnings"]
d = wait_job(c.post("/api/pipeline/run", json={**proc, "mode": "peptides", "text": "GILGFVFTL", "alleles": ["HLA-A*02:01"],
                                              "models": {"mhcflurry": True}}).get_json()["id"])
assert c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"][0]["MHCF_Processing"] == 0.1
# filters and ranking
d = wait_job(c.post("/api/pipeline/run", json={**proc, "alleles": ["HLA-A*02:01"], "models": {"mhcflurry": True},
                                              "filters": {"mhcf_pres_rank_max": 1, "mhcf_aff_max": 100}}).get_json()["id"])
assert d["summary"]["pairs_passed"] == 1, d["summary"]
# MHCflurry gates BigMHC when chained
d = wait_job(c.post("/api/pipeline/run", json={**proc, "alleles": ["HLA-A*02:01"], "chain": True,
                                              "models": {"mhcflurry": True, "bigmhc_el": True},
                                              "filters": {"mhcf_pres_rank_max": 1}}).get_json()["id"])
rows = {r["Peptide"]: r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
assert rows["GILGFVFTL"]["BigMHC_EL"] is not None and rows["MSLLTEVET"]["BigMHC_EL"] is None and d["summary"]["chained"], d["summary"]

if _sh.which("tcsh", path="/usr/bin:/bin") and _sh.which("gawk", path="/usr/bin:/bin"):
    seqq = SEQ[:30] + '"' + SEQ[31:]                                # a stray quote must not reach the bundle
    d = wait_job(c.post("/api/pipeline/run", json={**proc, "text": ">flu\n" + seqq, "lengths": [9, 12], "alleles": ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*99:99"],
                                                  "models": {"netctlpan": True, "netchop": True}}).get_json()["id"])
    assert d["state"] == "done", d
    rows = {(r["Peptide"], r["Allele"]): r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
    a2, b7 = rows[("GILGFVFTL", "HLA-A*02:01")], rows[("GILGFVFTL", "HLA-B*07:02")]
    assert a2["CTLpan_%Rank"] == 0.4 and a2["CTLpan_Combined"] == 0.96 and b7["CTLpan_%Rank"] == 0.4, (a2, b7)     # one job per allele, labelled right
    assert a2["NetChop_Cterm"] == 0.9, a2                      # cleavage after the peptide's last residue (L)
    w = d["summary"]["warnings"]
    assert any("NetCTLpan does not support: HLA-A*99:99" in x for x in w) and any("8-11-mers only" in x for x in w) and any("replaced by X" in x for x in w), w
    d = wait_job(c.post("/api/pipeline/run", json={**proc, "alleles": ["HLA-A*02:01"], "models": {"netctlpan": True},
                                                  "filters": {"ctlpan_rank_max": 1}}).get_json()["id"])
    assert d["summary"]["pairs_passed"] == 1, d["summary"]
    d = wait_job(c.post("/api/pipeline/run", json={**proc, "alleles": ["HLA-A*02:01"], "models": {"netchop": True},
                                                  "filters": {"netchop_min": 0.5}}).get_json()["id"])
    assert 0 < d["summary"]["pairs_passed"] < d["summary"]["pairs"], d["summary"]
    for name in ("mhcflurry", "iedb_bundle"):
        h = c.post(f"/api/health/{name}?deep=1").get_json(); assert h["ok"], h
else:
    print("skipped NetCTLpan/NetChop tests: tcsh and gawk are not installed")
# they need the protein: loose peptides are refused with the reason
r = c.post("/api/pipeline/run", json={**proc, "mode": "peptides", "text": "GILGFVFTL", "models": {"netctlpan": True}})
assert r.status_code == 400 and "need the protein" in r.get_json()["error"], r.get_json()
# an unready bundle is refused up front
nob = json.loads(json.dumps(settings)); nob["models"]["iedb_bundle"]["home"] = "/nonexistent"
c.put("/api/config", json={"text": json.dumps(nob)})
r = c.post("/api/pipeline/run", json={**proc, "models": {"netchop": True}})
assert r.status_code == 400 and "iedb_bundle is not ready" in r.get_json()["error"], r.get_json()
c.put("/api/config", json={"text": json.dumps(settings)})

os.environ.pop("PYTHONPATH", None)

# ---- optional: real BigMHC ----------------------------------------------------------
if REAL_PY and REAL_HOME:
    real = json.loads(json.dumps(settings))
    real["models"]["bigmhc"] = {"python": REAL_PY, "home": REAL_HOME}
    c.put("/api/config", json={"text": json.dumps(real)})
    h = c.post("/api/health/bigmhc?deep=1").get_json(); assert h["ok"], h
    r = c.post("/bigmhc/api/run", json={"text": "GILGFVFTL", "mode": "peptides",
                                        "alleles": ["HLA-A02:01"], "models": ["el", "im"]}).get_json()
    row = r["rows"][0]
    assert 0.2 < row["BigMHC_EL"] < 1 and 0 < row["BigMHC_IM"] < 1, row
    print("real BigMHC ok:", row["BigMHC_EL"], row["BigMHC_IM"])
print("all tests passed")
