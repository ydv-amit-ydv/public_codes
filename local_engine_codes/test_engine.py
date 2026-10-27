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

settings = {"tmp_dir": os.path.join(TMP, "work"),
            "models": {"netmhcpan": {"home": NMP}, "bigmhc": {"python": sys.executable, "home": BM}}}

# ---- discovery: copy-and-run layout (no explicit paths anywhere) ----------------------
import shutil                                  # noqa: E402
LE = os.path.join(TMP, "local_engines")
inner = os.path.join(LE, "netmhc", "netMHCpan-4.2", "netMHCpan-4.2")           # the tarball nests like this
os.makedirs(os.path.join(inner, "Linux_x86_64")); os.makedirs(os.path.join(inner, "data"))
shutil.copy(f"{NMP}/netMHCpan", f"{inner}/netMHCpan"); os.chmod(f"{inner}/netMHCpan", 0o755)
bm_home = os.path.join(LE, "bigmhc", "bigmhc"); os.makedirs(f"{bm_home}/models"); os.makedirs(f"{bm_home}/src")
shutil.copy(f"{BM}/src/predict.py", f"{bm_home}/src/predict.py")
os.makedirs(os.path.join(LE, "netmhc", "decoy", "netMHCpan"))                     # a *directory* named netMHCpan is not a script
os.environ["ENGINE_MODELS_DIR"] = LE
os.environ["BIGMHC_PYTHON"] = sys.executable
os.environ["EXPECTED_NMHOME"] = inner
from common import discover                    # noqa: E402
assert discover.find_netmhcpan(LE) == inner and discover.find_bigmhc(LE) == bm_home
assert discover.find_netmhcpan(os.path.join(TMP, "empty-nowhere")) is None
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
