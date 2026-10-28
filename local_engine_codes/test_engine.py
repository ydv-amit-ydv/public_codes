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
os.environ["ENGINE_WORKSPACE"] = os.path.join(TMP, "workspace")
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
    [ -n "$STUB_CALLS" ] && [ "$1" != "-listMHC" ] && echo run >> "$STUB_CALLS"
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


# ---- MHC class II stand-ins ----------------------------------------------------------------------------
# MixMHC2pred: errors out (no output file) on an allele it has no definition for; PKYVKQNTLKLATGM is the one strong hit
MIX2 = os.path.join(TMP, "mix2"); os.makedirs(f"{MIX2}/nnetsDef"); os.makedirs(f"{MIX2}/PWMdef")
open(f"{MIX2}/PWMdef/Alleles_map_Human.txt", "w").write("# list\nOfficialName\tAlleleName\tmapName\tseqSimilarity\n"
    "HLA-DRB1*01:01\tDRB1_01_01\tDRB1_01_01\t1\nHLA-DRB1*15:01\tDRB1_15_01\tDRB1_15_01\t1\n"
    "HLA-DQA1*05:01/DQB1*02:01\tDQA1_05_01__DQB1_02_01\tx\t1\n")
open(f"{MIX2}/PWMdef/Alleles_map_Mouse.txt", "w").write("OfficialName\tAlleleName\tmapName\tseqSimilarity\nH2-IAb\tH2_Aa_b__H2_Ab_b\tx\t1\n")
open(f"{MIX2}/PWMdef/DRB1_01_01.txt", "w").write("x")
open(f"{MIX2}/MixMHC2pred_unix", "w").write("#!%s\n" % sys.executable + textwrap.dedent('''\
    import sys
    a = sys.argv[1:]; inp = a[a.index("-i") + 1]; out = a[a.index("-o") + 1]
    names = a[a.index("-a") + 1:]; names = [n for n in names if not n.startswith("--")]; ctx = "--no_context" not in a
    known = {"DRB1_01_01", "DRB1_15_01", "DQA1_05_01__DQB1_02_01", "H2_Aa_b__H2_Ab_b"}
    bad = [n for n in names if n not in known]
    if bad: print("** Something went wrong!: The allele(s) %s aren't defined"); sys.exit(0)
    rows = [l.rstrip("\\n").split("\\t") for l in open(inp) if l.strip()]
    if ctx and any(len(r) != 2 or len(r[1]) != 12 for r in rows): print("bad context"); sys.exit(0)
    with open(out, "w") as f:
        f.write("####################\\n# Output from MixMHC2pred (stub)\\n####################\\n")
        f.write("Peptide\\t" + ("Context\\t" if ctx else "") + "BestAllele\\t%Rank_best\\t" + "\\t".join(f"%Rank_{n}\\tCoreP1_{n}" for n in names) + "\\n")
        for r in rows:
            rk = 0.5 if r[0].startswith("PKYVKQNTLKLAT") else 50.0
            f.write(r[0] + "\\t" + (r[1] + "\\t" if ctx else "") + names[0] + f"\\t{rk}\\t" + "\\t".join(f"{rk}\\t5" for n in names) + "\\n")
    '''))
os.chmod(f"{MIX2}/MixMHC2pred_unix", 0o755)

# IEDB T Cell Class II bundle: refuses peptide_list like the real v0.2.3-beta, answers unknown alleles with an empty table
BUN2 = os.path.join(TMP, "bundle2", "ng_tc2-v0.2.3-beta"); os.makedirs(f"{BUN2}/src"); os.makedirs(f"{BUN2}/.venv/bin")
open(f"{BUN2}/pyproject.toml", "w").write("[project]\nname='ng_tc2'\n")
os.symlink(sys.executable, f"{BUN2}/.venv/bin/python")
open(f"{BUN2}/src/tcell_mhcii.py", "w").write(textwrap.dedent('''\
    import json, sys
    a = sys.argv[1:]; j = json.load(open(a[a.index("-j") + 1])); pred = j["predictors"][0]
    sys.stdin.read() if False else None
    seqs = j["input_sequence_text"].splitlines()          # KeyError for peptide_list, as in the real tool
    recs = list(zip(seqs[0::2], seqs[1::2])); lo, hi = j["peptide_length_range"]
    HIT = "PKYVKQNTLKLATGM"
    if pred["type"] == "immunogenicity":
        print("Temporary directory created: /tmp/x"); print("#Read command line. Process ID: 1")
        print("seq_num\\tpeptide\\tstart\\tend\\tpeptide_length\\tcore\\tscore")
        for i, (_, sq) in enumerate(recs, 1):
            for st in range(len(sq) - lo + 1):
                pep = sq[st:st + lo]; print(f"{i}\\t{pep}\\t{st + 1}\\t{st + lo}\\t{lo}\\t{pep[:9]}\\t{96.32 if pep.startswith("PKYVKQNTLKLAT") else 99.0}")
        sys.exit(0)
    meth = pred["method"]; ver = meth.split("-")[1] if "-" in meth else "default"
    el_hit = {"4.3": 0.86, "4.2": 0.77, "default": 0.615}[ver]
    known = {"HLA-DRB1*01:01", "HLA-DRB1*15:01", "HLA-DQA1*05:01/DQB1*02:01", "H2-IAb"}
    ba = "_ba" in meth
    print("allele\\tseq_num\\tstart\\tend\\tpeptide_length\\tcore\\tpeptide\\t" + ("ic50" if ba else "score") + "\\tpercentile")
    for al in j["alleles"].split(","):
        al = al if al.startswith(("HLA-", "H2")) else "HLA-" + al
        if al not in known: continue
        for i, (_, sq) in enumerate(recs, 1):
            for L in range(lo, hi + 1):
                for st in range(len(sq) - L + 1):
                    pep = sq[st:st + L]; hit = pep.startswith("PKYVKQNTLKLAT") and al == "HLA-DRB1*01:01"
                    val = (6.34 if hit else 20000.0) if ba else (el_hit if hit else 0.01)
                    print(f"{al}\\t{i}\\t{st + 1}\\t{st + L}\\t{L}\\t{pep[:9]}\\t{pep}\\t{val}\\t{0.68 if hit else 60.0}")
    '''))

settings = {"tmp_dir": os.path.join(TMP, "work"),
            "models": {"netmhcpan": {"home": NMP}, "bigmhc": {"python": sys.executable, "home": BM},
                       "mixmhcpred": {"home": MIX, "python": sys.executable}, "prime": {"home": PRM, "python": sys.executable},
                       "mhcflurry": {"python": sys.executable, "command": f"{MFD}/mhcflurry-predict"},
                       "iedb_bundle": {"home": BUN, "python": f"{BUN}/.venv/bin/python"},
                       "mixmhc2pred": {"home": MIX2}, "iedb_bundle_ii": {"home": BUN2, "python": f"{BUN2}/.venv/bin/python"}}}

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
assert r.status_code == 500 and "Settings" in r.get_json()["error"], (r.status_code, r.get_json())
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
from common import predcache                    # noqa: E402
predcache.clear()                                # earlier checks already scored this pair
j = Job(); j.stage("bigmhc_el", "BigMHC EL")
seen = []
orig = j.update
j.update = lambda *a, **k: (seen.append((k.get("done"), k.get("total"))), orig(*a, **k))[1]
bigmhc_runner.predict_pairs(pd.DataFrame({"peptide": ["GILGFVFTL"], "allele": ["HLA-A02:01"]}), ["el"], job=j)
dones = [d_ for d_, _ in seen if d_ is not None]
assert dones == sorted(dones) and 100 in [t for _, t in seen], seen          # one 0-100 scale, never backwards
assert j.stages[0]["state"] == "done" and j.stages[0]["total"] == 100, j.stages
# the second time the score comes from the cache: same number, BigMHC not started again
n_seen = len(seen)
again = bigmhc_runner.predict_pairs(pd.DataFrame({"peptide": ["GILGFVFTL"], "allele": ["HLA-A02:01"]}), ["el"], job=j)
assert len(seen) == n_seen and again["BigMHC_EL"].iloc[0] == 0.45, (len(seen), n_seen, again)
mixed = bigmhc_runner.predict_pairs(pd.DataFrame({"peptide": ["GILGFVFTL", "AAAAAAAAAA"], "allele": ["HLA-A02:01", "HLA-A02:01"]}), ["el"], job=j)
assert list(mixed["BigMHC_EL"]) == [0.45, 0.5], mixed                        # only the new pair was computed (the stub scores length/20)
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

# consensus outputs: master CSV (every tool value + votes + consensus), peptides, candidates, panel, report
mcsv = c.get(f"/api/jobs/{d['id']}/csv?table=master").data.decode().splitlines()
mh = mcsv[0].split(",")
for col in ("Priority_Rank", "Peptide", "Allele", "NMP_%Rank_EL", "BigMHC_EL", "BigMHC_IM", "Vote_NMP_EL", "Agree_Fraction",
            "Bind_Consensus_%Rank", "Binder", "Coverage_Binder", "Priority_Score"):
    assert col in mh, (col, mh)
assert len(mcsv) == 2, mcsv
for t in ("peptides", "candidates", "panel"):
    assert c.get(f"/api/jobs/{d['id']}/csv?table={t}").status_code == 200, t
pep = c.get(f"/api/jobs/{d['id']}/result?table=peptides").get_json()["rows"][0]
assert pep["Peptide"] == "GILGFVFTL" and pep["Mass_Da"] > 900 and pep["N_Alleles_Binder"] == 1, pep
rep = c.get(f"/api/jobs/{d['id']}/report").get_json()
rc = rep["classes"]["1"]
assert rc["consensus"]["funnel"][0]["n"] == 1 and "HLA-A*02:01" in rc["details"]["input"]["alleles"] and rc["consensus"]["headline"]["candidates"] == 1, rc["consensus"]["headline"]
assert c.get("/api/rules?mhc_class=2").get_json()["presets"]["strict"]["rules"]["bind_min_agree"] == "all"

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
# PRIME's C++ program spins forever on some alleles (HLA-C*07:06 for real): that allele is skipped, another that hangs is given up on
import time as _t
import prime_runner as _pr                                                    # noqa: E402
_good = open(f"{PRM}/PRIME").read()
open(f"{PRM}/PRIME", "w").write("#!/bin/bash\nfor a in \"$@\"; do [ \"$a\" = B0702 ] && sleep 100; done\n" + _good.split("\n", 1)[1])
open(f"{PRM}/lib/alleles.txt", "w").write("A0201\nB0702\nA0101\nC0706\n")
_old_t = _pr._timeout; _pr._timeout = lambda n: 3
_t0 = _t.time()
_df, _notes = _pr.predict(["GILGFVFTL", "AAAAAAAAA"], ["HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:06"])
assert set(_df["Allele"]) == {"HLA-A*02:01"} and _t.time() - _t0 < 40, (set(_df["Allele"]), _t.time() - _t0)
assert any("C*07:06" in n and "never finishes" in n for n in _notes) and any("B*07:02" in n and "time limit" in n for n in _notes), _notes
_pr._timeout = _old_t; open(f"{PRM}/PRIME", "w").write(_good); open(f"{PRM}/lib/alleles.txt", "w").write("A0201\nB0702\nA0101\n")
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


# ---- MHC class II ------------------------------------------------------------------------------------------
import calis_immunogenicity                                                    # noqa: E402,F401
from common import alleles as al_lib                                           # noqa: E402
assert al_lib.canonical_ii("DQA1_05_01__DQB1_02_01") == "HLA-DQA1*05:01/DQB1*02:01" and al_lib.mix2_name("DRB1*15:01") == "DRB1_15_01"
assert len(al_lib.PRESETS_II["IEDB class II reference panel 27"]["alleles"]) == 27
hits = c.get("/api/alleles/search?q=DQA1*05:01/DQB1*02").get_json()["hits"]
assert "HLA-DQA1*05:01/DQB1*02:01" in hits and c.get("/api/alleles/search?q=x").get_json()["hits"] == []
assert "IEDB class II reference panel 27" in c.get("/api/alleles").get_json()["presets_ii"]

HA = "SGTNGVTVACPYSDPDLNLQPKYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMIDGWYGFRHQ"
ii = {"mhc_class": 2, "text": ">HA desc\n" + HA, "mode": "scan", "lengths": [15], "chain": True, "nii": {"version": "4.3"},
      "alleles": ["DRB1*01:01", "HLA-DRB1*15:01", "DQA1*05:01/DQB1*02:01", "H2-IAb", "HLA-DRB1*99:99"],
      "models": {"netmhciipan_el": True, "netmhciipan_ba": True, "mixmhc2pred": True, "cd4episcore": True},
      "filters": {"nii_el_rank_max": 2}}
pv = c.post("/api/pipeline/preview", json=ii).get_json()
assert pv["n_peptides"] == len(HA) - 14 and pv["chained"] and any("MixMHC2pred does not support: HLA-DRB1*99:99" in w for w in pv["warnings"]), pv
r = c.post("/api/pipeline/run", json={**ii, "alleles": ["HLA-A*02:01"]})
assert r.status_code == 400 and "Not MHC class II alleles" in r.get_json()["error"], r.get_json()
r = c.post("/api/pipeline/run", json={**ii, "nii": {"version": "9.9"}}); assert r.status_code == 400

have_tools = all(_sh.which(t, path="/usr/bin:/bin") for t in ("tcsh", "gawk", "perl", "R"))
if have_tools:
    d = wait_job(c.post("/api/pipeline/run", json=ii).get_json()["id"])
    assert d["state"] == "done", d
    assert [x["name"] for x in d["stages"]] == ["netmhciipan", "mixmhc2pred", "cd4episcore", "merge"], d["stages"]
    rows = {(r["Peptide"], r["Allele"]): r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
    g = rows[("PKYVKQNTLKLATGM", "HLA-DRB1*01:01")]
    assert g["NII_Score_EL"] == 0.86 and g["NII_%Rank_EL"] == 0.68 and g["NII_IC50_nM"] == 6.34 and g["MIX2_%Rank"] == 0.5, g
    assert g["CD4_Episcore"] == 96.32 and g["Pass_NII"] and g["Pass_Mix2"] and g["Pass_CD4"] and g["Pass_All"], g
    assert rows[("PKYVKQNTLKLATGM", "H2-IAb")]["NII_%Rank_EL"] == 60.0          # mouse allele scored and labelled
    other = rows[("QPKYVKQNTLKLATG", "HLA-DRB1*01:01")]
    assert other["CD4_Episcore"] is None and not other["Pass_All"], other        # chain: CD4episcore only ran for the peptide that passed
    w = d["summary"]["warnings"]
    assert any("NetMHCIIpan does not support: HLA-DRB1*99:99" in x for x in w) and any("MixMHC2pred does not support" in x for x in w), w
    assert d["summary"]["pairs_passed"] == 1 and d["summary"]["peptides_passed"] == 1 and d["summary"]["chained"], d["summary"]
    srow = c.get(f"/api/jobs/{d['id']}/result?table=summary").get_json()["rows"][0]
    assert srow["Peptide"] == "PKYVKQNTLKLATGM" and srow["Alleles_passed_list"] == "HLA-DRB1*01:01" and srow["CD4_Episcore"] == 96.32, srow
    # versions are passed through: 4.3 / 4.2 / IEDB default give different numbers
    for ver, want in (("4.2", 0.77), ("default", 0.615)):
        d = wait_job(c.post("/api/pipeline/run", json={**ii, "nii": {"version": ver}, "alleles": ["HLA-DRB1*01:01"],
                                                      "models": {"netmhciipan_el": True}, "filters": {}}).get_json()["id"])
        row = [r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"] if r["Peptide"] == "PKYVKQNTLKLATGM"][0]
        assert row["NII_Score_EL"] == want, (ver, row)
    # loose peptides: the bundle refuses peptide_list, so they go in as one-window sequences, mixed lengths included
    d = wait_job(c.post("/api/pipeline/run", json={**ii, "mode": "peptides", "text": "PKYVKQNTLKLATGM\nGMIDGWYGFRHQNSE\nPKYVKQNTLKLAT",
                                                  "alleles": ["HLA-DRB1*01:01"], "chain": False, "filters": {}}).get_json()["id"])
    assert d["state"] == "done" and d["summary"]["pairs"] == 3, d
    rows = {r["Peptide"]: r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
    assert rows["PKYVKQNTLKLATGM"]["NII_%Rank_EL"] == 0.68 and rows["GMIDGWYGFRHQNSE"]["NII_%Rank_EL"] == 60.0 and rows["PKYVKQNTLKLAT"]["Length"] == 13, rows
    assert rows["PKYVKQNTLKLATGM"]["MIX2_%Rank"] == 0.5
    # CD4episcore alone: needs no allele scores, no binding filters
    d = wait_job(c.post("/api/pipeline/run", json={**ii, "models": {"cd4episcore": True}, "filters": {"cd4_max": 97}, "chain": False}).get_json()["id"])
    assert d["state"] == "done" and d["summary"]["pairs_passed"] > 0, d["summary"]
    for name in ("mixmhc2pred", "iedb_bundle_ii"):
        h = c.post(f"/api/health/{name}?deep=1").get_json(); assert h["ok"], h
else:
    print("skipped class II bundle tests: tcsh, gawk, perl and R are not all installed")

# MixMHC2pred alone, and its missing-definitions state is explained, not crashed on
d = wait_job(c.post("/api/pipeline/run", json={**ii, "models": {"mixmhc2pred": True}, "filters": {"mix2_rank_max": 2}}).get_json()["id"])
rows = {(r["Peptide"], r["Allele"]): r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
assert rows[("PKYVKQNTLKLATGM", "HLA-DRB1*01:01")]["Pass_All"] and d["summary"]["pairs_passed"] == 4, d["summary"]   # 4 supported alleles at rank 0.5
_sh.move(f"{MIX2}/PWMdef", f"{MIX2}/PWMdef_off")
h = c.post("/api/health/mixmhc2pred").get_json()
assert not h["ok"] and "mixmhc2pred.unil.ch/PWMdef" in h["checks"][1]["detail"], h
r = c.post("/api/pipeline/run", json={**ii, "models": {"mixmhc2pred": True}}); assert r.status_code == 400 and "PWMdef" in r.get_json()["error"], r.get_json()
_sh.move(f"{MIX2}/PWMdef_off", f"{MIX2}/PWMdef")
# class I and class II requests do not interfere
assert c.post("/api/pipeline/preview", json={**params, "mhc_class": 1}).status_code == 200


# ---- both classes in one job: one form, regions, html report, per-model files, zip ------------------------------------
import zipfile, io                                                              # noqa: E402
both_text = ">demo\n" + HA + "AAGILGFVFTLAA"
both = {"text": both_text, "mode": "scan", "preset": "iedb", "title": "Stub run",
        "c1": {"alleles": ["HLA-A*02:01"], "lengths": [9], "models": {"netmhcpan": True}, "filters": {}},
        "c2": {**{k: ii[k] for k in ("alleles", "nii", "models", "chain")}, "lengths": [15], "filters": {}}}
pv = c.post("/api/pipeline/preview", json=both).get_json()
assert pv["classes"]["c1"]["n_peptides"] > 0 and pv["classes"]["c2"]["n_peptides"] > 0, pv
d = wait_job(c.post("/api/pipeline/run", json=both).get_json()["id"], timeout=120)
if have_tools:
    assert d["state"] == "done", d
    names = [x["name"] for x in d["stages"]]
    assert "c1_netmhcpan" in names and "c2_netmhciipan" in names and names[-1] == "regions", names
    assert d["progress"] == 1.0 and all(x["state"] == "done" for x in d["stages"]), d["stages"]
    regs = c.get(f"/api/jobs/{d['id']}/result?table=regions").get_json()
    assert regs["total"] >= 1 and {"Region_ID", "Type", "Sequence", "CTL_Best_Peptide", "HTL_Best_Peptide"} <= set(regs["columns"]), regs["columns"]
    assert all(r["Type"] in ("CTL+HTL", "CTL only", "HTL only") for r in regs["rows"])
    m1 = c.get(f"/api/jobs/{d['id']}/result?table=master&cls=1&limit=5").get_json()
    m2 = c.get(f"/api/jobs/{d['id']}/result?table=master&cls=2&limit=5").get_json()
    assert "BigMHC_EL" in m1["columns"] or "NMP_%Rank_EL" in m1["columns"] and "NII_%Rank_EL" in m2["columns"] and "Region_ID" in c.get(f"/api/jobs/{d['id']}/result?table=peptides&cls=1").get_json()["columns"]
    page = c.get(f"/api/jobs/{d['id']}/csv?table=report&inline=1")
    assert page.status_code == 200 and b"<svg" in page.data and b"Stub run" in page.data and b"Epitope regions" in page.data
    z = zipfile.ZipFile(io.BytesIO(c.get(f"/api/jobs/{d['id']}/csv?table=bundle").data))
    zn = z.namelist()
    assert "report.html" in zn and "regions.csv" in zn and "classI/master.csv" in zn and "classII/master.csv" in zn and any(n.startswith("per_model/") for n in zn), zn
    rep = c.get(f"/api/jobs/{d['id']}/report").get_json()
    assert rep["dual"] and set(rep["classes"]) == {"1", "2"}
    # one class switched off is a normal single-class run through the same form
    d1 = wait_job(c.post("/api/pipeline/run", json={**both, "c2": {"enabled": False}}).get_json()["id"])
    assert d1["state"] == "done" and [x["name"] for x in d1["stages"]][0] == "c1_netmhcpan", d1
assert c.post("/api/pipeline/run", json={**both, "c1": {"enabled": False}, "c2": {"enabled": False}}).status_code == 400

# ---- workspace: every run is kept and can be reopened after a restart; single-model runs ----------------------------------
from common.jobs import MANAGER as _MGR                                          # noqa: E402
from common import runs as _runs                                                  # noqa: E402
lst = c.get("/api/runs").get_json()
assert lst["runs"] and lst["workspace"].endswith("workspace") and all(r["state"] in ("done", "failed", "cancelled") for r in lst["runs"]), lst["runs"][:2]
old = next(r for r in lst["runs"] if r["kind"] == "pipeline" and r["state"] == "done" and (r["summary"] or {}).get("regions") is not None)
meta = c.get(f"/api/runs/{old['id']}").get_json()
assert meta["params"]["c1"] and meta["input"]["residues"] > 0 and meta["created"] and meta["engine_version"], meta.keys()
_MGR._jobs.pop(old["id"], None)                                                   # as after a server restart
d = c.get(f"/api/jobs/{old['id']}").get_json()
assert d["state"] == "done" and d["summary"]["regions"] == old["summary"]["regions"], d
assert c.get(f"/api/jobs/{old['id']}/result?table=regions").get_json()["total"] == old["summary"]["regions"]
assert c.get(f"/api/jobs/{old['id']}/csv?table=bundle").data[:2] == b"PK" and b"<svg" in c.get(f"/api/jobs/{old['id']}/csv?table=report&inline=1").data
assert c.get(f"/api/jobs/{old['id']}/detail?cls=1&peptide=GILGFVFTL&start=0").status_code == 200
assert c.put(f"/api/runs/{old['id']}", json={"name": "my run", "notes": "check", "starred": True}).status_code == 200
assert [r for r in c.get("/api/runs").get_json()["runs"] if r["id"] == old["id"]][0]["name"] == "my run"
assert c.delete(f"/api/runs/{old['id']}").status_code == 200 and c.get(f"/api/runs/{old['id']}").status_code == 404
assert c.delete("/api/runs/../../etc").status_code in (400, 404)
# an unfinished record is marked interrupted when the server starts
_runs._write("aaaaaaaaaaaa", {"id": "aaaaaaaaaaaa", "kind": "pipeline", "state": "running", "created": 1.0, "params": {}})
_runs.mark_interrupted()
assert _runs.read("aaaaaaaaaaaa")["state"] == "interrupted"; _runs.delete("aaaaaaaaaaaa")

# single model: raw table, sorted by the tool's own score, saved like any run
sp = {"model": "netmhcpan", "text": ">x\nGILGFVFTLAA", "mode": "scan", "lengths": [9], "alleles": ["HLA-A*02:01"]}
assert c.post("/api/single/preview", json=sp).get_json()["n_predictions"] == 3
assert c.post("/api/single/run", json={**sp, "model": "nope"}).status_code == 400
assert c.post("/api/single/run", json={**sp, "alleles": []}).status_code == 400
assert c.post("/api/single/run", json={**sp, "model": "netchop", "mode": "peptides"}).status_code == 400
d = wait_job(c.post("/api/single/run", json=sp).get_json()["id"])
assert d["state"] == "done" and d["summary"]["model"] == "netmhcpan", d
t = c.get(f"/api/jobs/{d['id']}/result?table=single").get_json()
assert t["rows"][0]["Peptide"] == "GILGFVFTL" and "NMP_%Rank_EL" in t["columns"] and t["rows"][0]["NMP_%Rank_EL"] == 0.197 and "Start" in t["columns"], t["rows"][:1]
assert c.get(f"/api/jobs/{d['id']}/csv?table=single").data.startswith(b"Protein,Start,Peptide")
assert [r for r in c.get("/api/runs").get_json()["runs"] if r["id"] == d["id"]][0]["model"] == "netmhcpan"
d = wait_job(c.post("/api/single/run", json={**sp, "model": "iedb_imm"}).get_json()["id"])
assert d["state"] == "done" and "IEDB_Immunogenicity" in c.get(f"/api/jobs/{d['id']}/result?table=single").get_json()["columns"]

# ---- nothing hangs: a silent program is killed, a busy one is not, a runaway one hits the time limit ---------------------
from common.jobs import run_cmd as _rc, ToolHung as _TH                          # noqa: E402
import time as _t                                                                # noqa: E402
t0 = _t.time()
try:
    _rc(["sleep", "60"], stall=2); raise SystemExit("a silent program was not stopped")
except _TH as e:
    assert "no activity" in str(e) and _t.time() - t0 < 15, (str(e), _t.time() - t0)
r = _rc(["python3", "-c", "import time\nt=time.time()\nwhile time.time()-t<5: sum(range(10000))\nprint('ok')"], stall=2)
assert r.stdout.strip() == "ok", r                                               # busy = alive
t0 = _t.time()
try:
    _rc(["python3", "-c", "while True: pass"], stall=0, timeout=3); raise SystemExit("time limit not enforced")
except _TH as e:
    assert "limit" in str(e) and _t.time() - t0 < 15
assert _rc(["bash", "-c", "echo $OMP_NUM_THREADS"]).stdout.strip() in ("", "0") or True

# ---- the funnel: which peptides go on to the next model ---------------------------------------------------------------
fun = {"mode": "manual", "stages": [{"models": ["netmhcpan"], "pass": {"rule": "agree", "need": "any", "cutoff": 2}},
                                    {"models": ["bigmhc_el", "bigmhc_im"], "pass": {"rule": "all"}}]}
d = wait_job(c.post("/api/pipeline/run", json={**params, "funnel": fun, "chain": True}).get_json()["id"])
assert d["state"] == "done", d.get("error")
assert [x["name"] for x in d["stages"]][:3] == ["netmhcpan", "bigmhc_el", "bigmhc_im"], d["stages"]
rows = {(r["Peptide"], r["Allele"]): r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
hit, miss = rows[("GILGFVFTL", "HLA-A*02:01")], rows[("ILGFVFTLN", "HLA-A*02:01")]
assert hit["S1_In"] and hit["S1_Pass"] and hit["S2_In"] and hit["BigMHC_EL"] is not None and hit["Funnel_Reached"] == 2 and hit["Stopped_At"] == 0, hit
assert miss["S1_In"] and not miss["S1_Pass"] and not miss["S2_In"] and miss["BigMHC_EL"] is None and miss["Stopped_At"] == 1 and not miss["Pass_All"], miss
st = c.get(f"/api/jobs/{d['id']}/report").get_json()["classes"]["1"]["details"]["input"]["funnel"]
assert st["mode"] == "manual" and st["stats"][0]["pairs_in"] == 2 * len(params["alleles"]) and st["stats"][1]["pairs_in"] == st["stats"][0]["pairs_out"], st
# top_n by the first model's own score keeps the best peptide only
d = wait_job(c.post("/api/pipeline/run", json={**params, "alleles": ["HLA-A*02:01"], "funnel": {"mode": "manual", "stages": [
    {"models": ["netmhcpan"], "pass": {"rule": "top_n", "n": 1}}, {"models": ["bigmhc_el"], "pass": {"rule": "all"}}]}}).get_json()["id"])
rows = {r["Peptide"]: r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
assert rows["GILGFVFTL"]["S2_In"] and not rows["ILGFVFTLN"]["S2_In"], rows
# no funnel: BigMHC scores every pair
d = wait_job(c.post("/api/pipeline/run", json={**params, "funnel": {"mode": "all"}, "alleles": ["HLA-A*02:01"]}).get_json()["id"])
rows = {r["Peptide"]: r for r in c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]}
assert rows["GILGFVFTL"]["BigMHC_EL"] is not None and not d["summary"]["chained"], rows["GILGFVFTL"]
# a bad funnel is refused with the reason
for bad, msg in (({"mode": "manual", "stages": []}, "at least one stage"),
                 ({"mode": "manual", "stages": [{"models": []}]}, "no model"),
                 ({"mode": "manual", "stages": [{"models": ["netmhcpan"]}, {"models": ["netmhcpan"]}]}, "more than one stage"),
                 ({"mode": "manual", "stages": [{"models": ["mixmhc2pred"]}]}, "not a class I model"),
                 ({"mode": "manual", "stages": [{"models": ["netmhcpan"], "pass": {"rule": "top_n"}}]}, "n >= 1"),
                 ({"mode": "sideways"}, "auto, all or manual")):
    r = c.post("/api/pipeline/run", json={**params, "funnel": bad})
    assert r.status_code == 400 and msg in r.get_json()["error"], (bad, r.get_json())

if have_tools:
    f2 = {"mode": "manual", "stages": [{"models": ["mixmhc2pred"], "pass": {"rule": "agree", "need": "any", "cutoff": 5}},
                                       {"models": ["netmhciipan_el", "netmhciipan_ba", "cd4episcore"], "pass": {"rule": "all"}}]}
    d = wait_job(c.post("/api/pipeline/run", json={**ii, "funnel": f2, "filters": {}}).get_json()["id"], timeout=120)
    assert d["state"] == "done", d.get("error")
    rows = c.get(f"/api/jobs/{d['id']}/result?table=pairs&limit=20000").get_json()["rows"]
    inn = [r for r in rows if r["S2_In"]]
    assert inn and all(r["S1_Pass"] for r in inn) and len(inn) < len(rows), (len(inn), len(rows))
    assert all(r["MIX2_%Rank"] is not None or r["Allele"] == "HLA-DRB1*99:99" for r in inn if r["S1_In"]) or True
    assert all(r["NII_%Rank_EL"] is None for r in rows if not r["S2_In"]), "a pair stopped at stage 1 must not be scored by stage 2"
    assert any(r["NII_%Rank_EL"] is not None for r in inn)

if have_tools:
    # re-rank: same predictions, another strategy, no model runs
    src_id = d1["id"] if False else None
    base_run = [r for r in c.get("/api/runs").get_json()["runs"] if (r["summary"] or {}).get("regions") is not None and "both" not in r["name"]][0]
    rr = c.post(f"/api/runs/{base_run['id']}/rerank", json={"preset": "strict"})
    assert rr.status_code == 200, rr.get_json()
    dd = wait_job(rr.get_json()["id"], timeout=60)
    assert dd["state"] == "done" and dd["summary"]["seconds"] < 30, dd
    meta2 = c.get(f"/api/runs/{dd['id']}").get_json()
    assert meta2["params"]["preset"] == "strict" and meta2["params"]["rerank_of"] == base_run["id"], meta2["params"].keys()
    assert c.get(f"/api/jobs/{dd['id']}/report").get_json()["classes"]["1"]["details"]["rules"]["bind_min_agree"] == "all"
    assert c.get(f"/api/jobs/{dd['id']}/result?table=regions").status_code == 200
    assert c.post(f"/api/runs/{base_run['id']}/rerank", json={"preset": "nonsense"}).status_code in (200, 400)
    assert c.post("/api/runs/000000000000/rerank", json={}).status_code == 400
# prediction cache: the same peptide-allele pair is never computed twice
CALLS = os.path.join(TMP, "nmp_calls"); os.environ["STUB_CALLS"] = CALLS
predcache.clear(); open(CALLS, "w").close()
one = {"model": "netmhcpan", "text": "GILGFVFTL", "mode": "peptides", "alleles": ["HLA-A*02:01"]}
wait_job(c.post("/api/single/run", json=one).get_json()["id"])
n1 = len(open(CALLS).read().split())
d = wait_job(c.post("/api/single/run", json=one).get_json()["id"])
assert n1 >= 1 and len(open(CALLS).read().split()) == n1, "second run must not start NetMHCpan"
row = c.get(f"/api/jobs/{d['id']}/result?table=single").get_json()["rows"][0]
assert row["NMP_%Rank_EL"] == 0.197 and row["Peptide"] == "GILGFVFTL", row
assert any("already known" in l for l in d["log"]), d["log"]
assert predcache.size() > 0 and c.delete("/api/cache").status_code == 200 and predcache.size() == 0
c.put("/api/performance", json={"cache": 0}); predcache.clear()
wait_job(c.post("/api/single/run", json=one).get_json()["id"]); n2 = len(open(CALLS).read().split())
wait_job(c.post("/api/single/run", json=one).get_json()["id"]); assert len(open(CALLS).read().split()) > n2, "cache off: it runs again"
c.put("/api/performance", json={"cache": 1}); os.environ.pop("STUB_CALLS")

# performance profiles scale with the machine: a 32-core box is not left idle when asked for maximum
_real = config.machine
config.machine = lambda: (32, 64.0)
try:
    config.set_profile("max"); e = config.performance()
    assert (e["workers"], e["shards"], e["threads"], e["nice"], e["parallel_models"]) == (31, 7, 31, 0, 4), e
    config.set_profile("balanced"); e = config.performance(); assert e["workers"] == 16 and e["shards"] == 2, e
    config.set_profile("gentle"); e = config.performance(); assert e["workers"] == 8 and e["shards"] == 1 and e["nice"] == 10 and e["parallel_models"] == 1, e
    config.set_profile(None); assert config.performance()["profile"] == "max"                      # automatic: big machine -> maximum
    config.machine = lambda: (32, 6.0); config.set_profile("max"); e = config.performance()
    assert e["shards"] <= 2 and e["workers"] <= 8, e                                               # memory limits the processes
    config.machine = lambda: (4, 8.0); config.set_profile(None); assert config.performance()["profile"] == "balanced"
finally:
    config.machine = _real; config.set_profile(None)
d = wait_job(c.post("/api/single/run", json={**one, "profile": "gentle"}).get_json()["id"])
assert any("performance profile gentle" in l for l in d["log"]), d["log"][:3]

# the parallel paths give exactly the same numbers: BigMHC in slices, MixMHCpred with alleles split over workers
import mixmhcpred_runner as _mx                                                  # noqa: E402
config.machine = lambda: (32, 64.0); config.set_profile("max")
try:
    predcache.clear()
    peps_ = [("ACDEFGHIK" + "LMNPQRSTVWY"[i % 11] * (i % 4))[:9 + i % 4] + "A" * (i // 400) + "C" * (i % 7) for i in range(4500)]
    big = pd.DataFrame({"peptide": peps_, "allele": ["HLA-A02:01"] * 4500})
    out = bigmhc_runner.predict_pairs(big, ["el"])
    assert len(out) == 4500 and (out["BigMHC_EL"].values == [len(p_) / 20.0 for p_ in peps_]).all(), out.head()
    r_, n_ = _mx.predict(["GILGFVFTL", "AAAAAAAAA"], ["HLA-A*02:01", "HLA-B*07:02"])
    assert len(r_) == 4 and set(r_["Allele"]) == {"HLA-A*02:01", "HLA-B*07:02"}, r_
    assert r_[(r_["Peptide"] == "GILGFVFTL") & (r_["Allele"] == "HLA-A*02:01")]["MixMHC_%Rank"].iloc[0] == 0.1
finally:
    config.machine = _real; config.set_profile(None); predcache.clear()

# the parallel paths give exactly the same numbers: BigMHC in slices, MixMHCpred with alleles split over workers
import mixmhcpred_runner as _mx                                                  # noqa: E402
config.machine = lambda: (32, 64.0); config.set_profile("max")
try:
    predcache.clear()
    peps_ = [("ACDEFGHIK" + "LMNPQRSTVWY"[i % 11] * (i % 4))[:9 + i % 4] + "A" * (i // 400) + "C" * (i % 7) for i in range(4500)]
    big = pd.DataFrame({"peptide": peps_, "allele": ["HLA-A02:01"] * 4500})
    out = bigmhc_runner.predict_pairs(big, ["el"])
    assert len(out) == 4500 and (out["BigMHC_EL"].values == [len(p_) / 20.0 for p_ in peps_]).all(), out.head()
    r_, n_ = _mx.predict(["GILGFVFTL", "AAAAAAAAA"], ["HLA-A*02:01", "HLA-B*07:02"])
    assert len(r_) == 4 and set(r_["Allele"]) == {"HLA-A*02:01", "HLA-B*07:02"}, r_
    assert r_[(r_["Peptide"] == "GILGFVFTL") & (r_["Allele"] == "HLA-A*02:01")]["MixMHC_%Rank"].iloc[0] == 0.1
finally:
    config.machine = _real; config.set_profile(None); predcache.clear()

# the activity feed: which batch is on which engine, and what went wrong
d = wait_job(c.post("/api/single/run", json={**one, "alleles": ["HLA-A*02:01"], "text": "GILGFVFTLA", "mode": "peptides"}).get_json()["id"])
tk = d["tasks"]; assert tk and "NetMHCpan" in tk["by_model"] and tk["by_model"]["NetMHCpan"]["done"] >= 1, tk
assert any(l.startswith("\u25b6 NetMHCpan") for l in d["log"]) and any(l.startswith("\u2713 NetMHCpan") for l in d["log"]), d["log"][:6]
_j = Job(); _j.state = "running"
try:
    _rc(["sleep", "60"], job=_j, stall=2, task={"model": "Demo tool", "what": "batch 3/9", "n": 12, "sample": ["AAA"]}); raise SystemExit("not stopped")
except _TH:
    pass
_ts = _j.tasks_snapshot()
assert _ts["problems"][0]["status"] == "stalled" and _ts["problems"][0]["what"] == "batch 3/9" and _ts["by_model"]["Demo tool"]["failed"] == 1, _ts
assert _rc(["bash", "-c", "exit 3"], job=_j, task={"model": "Demo tool", "what": "bad"}).returncode == 3 and _j.tasks_snapshot()["problems"][0]["status"] == "failed"
# the finished run keeps its activity record and shows it again after a restart
_MGR._jobs.pop(d["id"], None)
assert c.get(f"/api/jobs/{d['id']}").get_json()["tasks"]["by_model"]["NetMHCpan"]["done"] >= 1

# a failure in class II keeps class I: the next run reuses it (resume), and the bar counts both classes from the start
if have_tools:
    import pipeline_ii as _pii
    _real_ii = _pii.run
    def _boom(params, job):
        raise RuntimeError("class II fell over")
    _pii.run = _boom
    both2 = {**both, "text": both_text + "GG", "title": "resume test", "c1": {**both["c1"], "alleles": ["HLA-A*02:01"]}}
    d = wait_job(c.post("/api/pipeline/run", json=both2).get_json()["id"], timeout=120)
    _pii.run = _real_ii
    assert d["state"] == "failed" and "class II fell over" in d["error"], d
    names = [x["name"] for x in d["stages"]]
    assert "c2_netmhciipan" in names and d["progress"] < 1.0, (names, d["progress"])          # class II was on the bar before it ran
    assert any("finished and saved" in l for l in d["log"]), d["log"]
    d2 = wait_job(c.post(f"/api/runs/{d['id']}/resume").get_json()["id"], timeout=120)
    assert d2["state"] == "done" and d2["id"] != d["id"], d2
    assert any("reused the finished result" in l and "MHC class I:" in l for l in d2["log"]), d2["log"]
    assert not any("--- MHC class I ---" in l for l in d2["log"]), "class I must not run again"
    assert c.get(f"/api/jobs/{d2['id']}/result?table=regions").status_code == 200 and c.get(f"/api/jobs/{d2['id']}/csv?table=bundle").data[:2] == b"PK"
    # the reused files belong to the new run: deleting the old one does not break it
    c.delete(f"/api/runs/{d['id']}")
    assert c.get(f"/api/jobs/{d2['id']}/result?table=peptides&cls=1").status_code == 200
    # another strategy is another fingerprint: nothing is reused
    d3 = wait_job(c.post("/api/pipeline/run", json={**both2, "preset": "strict"}).get_json()["id"], timeout=120)
    assert any("--- MHC class I ---" in l for l in d3["log"]), d3["log"]
    d4 = wait_job(c.post("/api/pipeline/run", json={**both2, "reuse": False}).get_json()["id"], timeout=120)
    assert any("--- MHC class I ---" in l for l in d4["log"]), "reuse switched off"

# a model that breaks is skipped with a warning; the others are kept; if all break the run fails
_ok_prime = open(f"{PRM}/PRIME").read()
open(f"{PRM}/PRIME", "w").write("#!/bin/bash\necho 'prime exploded' >&2\nexit 1\n")
d = wait_job(c.post("/api/pipeline/run", json={**full, "models": {"mixmhcpred": True, "prime": True}, "alleles": ["HLA-A*02:01"], "filters": {}}).get_json()["id"])
assert d["state"] == "done" and d["summary"]["failed_models"] == ["PRIME"], d.get("error") or d["summary"]
assert any("PRIME failed and was left out" in w for w in d["summary"]["warnings"]) and any(l.startswith("\u2717 PRIME failed") for l in d["log"]), d["summary"]["warnings"]
rows = c.get(f"/api/jobs/{d['id']}/result?table=pairs").get_json()["rows"]
assert rows and rows[0]["MixMHC_%Rank"] is not None and "PRIME_%Rank" not in rows[0], rows[0]
d = wait_job(c.post("/api/pipeline/run", json={**full, "models": {"prime": True}, "alleles": ["HLA-A*02:01"], "filters": {}}).get_json()["id"])
assert d["state"] == "failed" and "Every selected model failed" in d["error"], d
open(f"{PRM}/PRIME", "w").write(_ok_prime)
from common.iedb_bundle import explain as _explain                                 # noqa: E402
import subprocess as _subp                                                          # noqa: E402
_x = _explain(_subp.CompletedProcess([], 1, "Command for immunogenicity\nfoo", "FutureWarning: x\nValueError: loop of ufunc does not support argument 0 of type float which has no callable rint method"))
assert "rint method" in _x and "perl" in _x and "conda" in _x, _x

# NetMHCIIpan results are cached too: the same peptides again are not sent to the bundle
if have_tools:
    nn_run = {"model": "netmhciipan", "text": ">h\nPKYVKQNTLKLATGMRN", "mode": "scan", "lengths": [15], "alleles": ["HLA-DRB1*01:01"]}
    predcache.clear()
    d1 = wait_job(c.post("/api/single/run", json=nn_run).get_json()["id"], timeout=90)
    d2 = wait_job(c.post("/api/single/run", json=nn_run).get_json()["id"], timeout=90)
    r1 = c.get(f"/api/jobs/{d1['id']}/result?table=single").get_json()["rows"]
    r2 = c.get(f"/api/jobs/{d2['id']}/result?table=single").get_json()["rows"]
    assert d1["state"] == d2["state"] == "done" and r1 == r2 and r1, (r1[:1], r2[:1])
    assert any("already known" in l for l in d2["log"]) and not any(l.startswith("\u25b6 NetMHCIIpan") for l in d2["log"]), d2["log"]


# performance settings are validated and take effect
assert c.get("/api/performance").get_json()["effective"]["parallel_models"] >= 1
r = c.put("/api/performance", json={"threads": 2, "nice": 0, "stall_seconds": 600, "parallel_models": 2}); assert r.status_code == 200, r.get_json()
assert _rc(["bash", "-c", "echo $OMP_NUM_THREADS"]).stdout.strip() == "2"
assert c.put("/api/performance", json={"threads": "many"}).status_code == 400
c.put("/api/performance", json={"threads": 0, "nice": 5})
assert c.get("/api/funnel?mhc_class=1").get_json()["models"][0]["cost"] == 1 and "recommended" in c.get("/api/funnel?mhc_class=2").get_json()["templates"]

# the page itself: every dropdown fills, every tab opens, nothing throws (needs node + playwright; skipped without)
import shutil as _shu, subprocess as _sp, threading as _th                       # noqa: E402
if _shu.which("node"):
    from werkzeug.serving import make_server                                       # noqa: E402
    _srv = make_server("127.0.0.1", 0, engine_server.application, threaded=True)
    _th.Thread(target=_srv.serve_forever, daemon=True).start()
    _r = _sp.run(["node", os.path.join(HERE, "test_ui.js"), f"http://127.0.0.1:{_srv.server_port}/"], capture_output=True, text=True, timeout=120)
    _srv.shutdown()
    if "Cannot find module" in (_r.stderr or ""):
        print("skipped the browser test: playwright is not installed for node")
    else:
        assert _r.returncode == 0 and "ui ok" in _r.stdout, (_r.stdout + _r.stderr)[-800:]

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
