"""Epitope Filter Engine server: one page with a tab per model plus a Settings tab
whose JSON (engine_config.json) can be edited live.

    python3 engine_server.py          # http://127.0.0.1:5000

Each model keeps its own interpreter/install (see common/config.py); this process
only needs flask and pandas. Settings are re-read on every request, so saving in
the GUI takes effect immediately.
"""
import json
import os, subprocess, sys

from flask import Flask, jsonify, request, send_file
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.serving import run_simple

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ("MHC1/BA_EL_scripts", "MHC1/Immunogenicity_scripts", "MHC1/Processing_scripts",
            "MHC2/BA_EL_scripts", "MHC2/Immunogenicity_scripts"):
    sys.path.insert(0, os.path.join(HERE, sub))
sys.path.insert(0, HERE)

from common import alleles as allele_lib, config, population            # noqa: E402
from common.jobs import MANAGER                              # noqa: E402
from common import runs, predcache                           # noqa: E402
import netmhcpan_app, bigmhc_app, bigmhc_runner as bm       # noqa: E402
import mixmhcpred_runner as mx, prime_runner as pr           # noqa: E402
import mhcflurry_runner as mf, iedb_processing as ip         # noqa: E402
import funnel, models as model_registry
import pipeline, pipeline_ii, pipeline_dual, pipeline_single                                  # noqa: E402
import mixmhc2pred_runner as mx2, netmhciipan_runner as nn2, cd4episcore_runner as cd4, mhcnuggets_runner as mn2   # noqa: E402
import pandas as pd                                          # noqa: E402

main = Flask(__name__)


def _check(name, ok, detail=""):
    return {"name": name, "ok": bool(ok), "detail": str(detail)[:600]}


def _auto(m, key):
    return "  (auto-detected)" if key in m.get("auto", []) else ""


def health_netmhcpan(deep):
    m = config.model("netmhcpan")
    script = netmhcpan_app.netmhc_path()
    checks = [_check("install dir", os.path.isdir(m.get("home") or ""),
                     (m.get("home") or f"not found under {config.search_root()}") + _auto(m, "home"))]
    checks.append(_check("netMHCpan script executable", os.access(script, os.X_OK), script))
    if not checks[-1]["ok"]:
        return checks
    netmhcpan_app._list_alleles.cache_clear()
    n = len(netmhcpan_app.list_alleles())
    checks.append(_check("-listMHC", n > 0, f"{n} alleles" if n else "no output: is tcsh installed and NMHOME set in the script?"))
    if deep and n:
        res, _ = netmhcpan_app.run_netmhcpan(["GILGFVFTL"], ["HLA-A02:01"],
                                             dict(rankS=0.5, rankW=2.0, ba=True, t=None, raw=""))
        hit = res[res["Peptide"] == "GILGFVFTL"] if not res.empty else res
        got = (not hit.empty) and "SB" in str(hit["NMP_BindLevel"].iloc[0])
        checks.append(_check("GILGFVFTL / A02:01 -> strong binder", got,
                             "ok" if got else "unexpected or empty result"))
    return checks


def health_bigmhc(deep):
    m = config.model("bigmhc")
    py, home = m.get("python", ""), m.get("home", "")
    checks = [_check("python (has torch)" if "python" in m["auto"] else "python", os.path.isfile(py or ""),
                     (py or "no env with torch found; create the bigmhc env or set models.bigmhc.python") + _auto(m, "python")),
              _check("predict.py", os.path.isfile(os.path.join(home, "src", "predict.py")),
                     (home or f"BigMHC clone not found under {config.search_root()}") + _auto(m, "home"))]
    if not all(c["ok"] for c in checks):
        return checks
    r = subprocess.run([py, "-c", "import torch, numpy, pandas, psutil; print('torch', torch.__version__)"],
                       capture_output=True, text=True)
    checks.append(_check("imports (torch, numpy, pandas, psutil)", r.returncode == 0,
                         r.stdout.strip() or r.stderr.strip().splitlines()[-1:]))
    if deep and checks[-1]["ok"]:
        try:
            g = bm.predict(bm.given_peptides(["GILGFVFTL"]), ["HLA-A*02:01"], ["el", "im"])
            checks.append(_check("GILGFVFTL / A*02:01 prediction", True,
                                 f"EL {g['BigMHC_EL'].iloc[0]:.3f}  IM {g['BigMHC_IM'].iloc[0]:.3f}"))
        except Exception as e:
            checks.append(_check("GILGFVFTL / A*02:01 prediction", False, e))
    return checks


def _py_imports(py, mods):
    r = subprocess.run([py, "-c", "import " + ", ".join(mods)], capture_output=True, text=True)
    return r.returncode == 0, (r.stderr.strip().splitlines() or ["ok"])[-1]


def health_mixmhcpred(deep):
    m = config.model("mixmhcpred")
    home, py = m.get("home") or "", m.get("python") or ""
    checks = [_check("install dir", os.path.isfile(os.path.join(home, "MixMHCpred")),
                     (home or f"not found under {config.search_root()} (clone GfellerLab/MixMHCpred there)") + _auto(m, "home")),
              _check("python (numpy, pandas)", bool(py), (py or "none found") + _auto(m, "python"))]
    if not all(c["ok"] for c in checks):
        return checks
    ok, detail = _py_imports(py, ("numpy", "pandas"))
    checks.append(_check("imports (numpy, pandas)", ok, detail if not ok else "ok"))
    if deep and ok:
        try:
            df, _ = mx.predict(["GILGFVFTL"], ["HLA-A*02:01"])
            good = len(df) == 1 and float(df["MixMHC_%Rank"].iloc[0]) < 1.0
            checks.append(_check("GILGFVFTL / A*02:01 prediction", good,
                                 f"%Rank {float(df['MixMHC_%Rank'].iloc[0]):.3f}" if len(df) else "no result"))
        except Exception as e:
            checks.append(_check("GILGFVFTL / A*02:01 prediction", False, e))
    return checks


def health_prime(deep):
    m, mix = config.model("prime"), config.model("mixmhcpred")
    home, py = m.get("home") or "", m.get("python") or ""
    checks = [_check("install dir", os.path.isfile(os.path.join(home, "lib", "PRIME.cc")),
                     (home or f"not found under {config.search_root()} (clone GfellerLab/PRIME there)") + _auto(m, "home")),
              _check("MixMHCpred 3.0 (required by PRIME)", bool(mix.get("home")),
                     mix.get("home") or "not found (clone GfellerLab/MixMHCpred beside PRIME)")]
    if not all(c["ok"] for c in checks):
        return checks
    try:
        elf = open(os.path.join(home, "lib", "PRIME.x"), "rb").read(4) == b"\x7fELF"
    except OSError:
        elf = False
    import shutil
    checks.append(_check("C++ program", elf or bool(shutil.which("g++")),
                         "built for Linux" if elf else ("will be compiled on first use" if shutil.which("g++")
                                                         else "needs compiling and g++ is missing: sudo apt install g++")))
    if not py:
        return checks + [_check("python (numpy, pandas)", False, "none found")]
    if deep and checks[-1]["ok"]:
        try:
            df, _ = pr.predict(["GILGFVFTL"], ["HLA-A*02:01"])
            good = len(df) == 1 and float(df["PRIME_%Rank"].iloc[0]) < 5.0
            checks.append(_check("GILGFVFTL / A*02:01 prediction", good,
                                 f"PRIME %Rank {float(df['PRIME_%Rank'].iloc[0]):.3f}" if len(df) else "no result"))
        except Exception as e:
            checks.append(_check("GILGFVFTL / A*02:01 prediction", False, e))
    return checks


def health_mhcflurry(deep):
    m = config.model("mhcflurry")
    py = m.get("python") or ""
    checks = [_check("python (has mhcflurry)", bool(py),
                     (py or "no environment with mhcflurry found: see README, or set models.mhcflurry.python") + _auto(m, "python"))]
    if not py:
        return checks
    cmd = m.get("command")
    cmd = os.path.join(os.path.dirname(py), "mhcflurry-predict") if not cmd or cmd == "auto" else cmd
    checks.append(_check("mhcflurry-predict", os.path.isfile(cmd), cmd))
    r = subprocess.run([py, "-c", "import os, mhcflurry; from mhcflurry.downloads import get_path; "
                        "p = get_path('models_class1_presentation', test_exists=False); "
                        "print(p); raise SystemExit(0 if os.path.isdir(p) else 3)"], capture_output=True, text=True)
    checks.append(_check("models_class1_presentation downloaded", r.returncode == 0,
                         "ok" if r.returncode == 0 else "run: mhcflurry-downloads fetch models_class1_presentation"
                         if r.returncode == 3 else (r.stderr.strip().splitlines() or ["import mhcflurry failed"])[-1]))
    if deep and all(c["ok"] for c in checks):
        try:
            df, _ = mf.predict(pd.DataFrame({"Protein": ["t"], "Start": [1], "Peptide": ["GILGFVFTL"]}), ["HLA-A*02:01"])
            good = len(df) == 1 and float(df["MHCF_Affinity_nM"].iloc[0]) < 100
            checks.append(_check("GILGFVFTL / A*02:01 prediction", good,
                                 f"{float(df['MHCF_Affinity_nM'].iloc[0]):.1f} nM" if len(df) else "no result"))
        except Exception as e:
            checks.append(_check("GILGFVFTL / A*02:01 prediction", False, e))
    return checks


def health_iedb_bundle(deep):
    import shutil
    m = config.model("iedb_bundle")
    home, py = m.get("home") or "", m.get("python") or ""
    checks = [_check("IEDB bundle folder", os.path.isfile(os.path.join(home, "src", "tcell_mhci.py")),
                     (home or f"not found under {config.search_root()}: download IEDB_NG_TC1-*.tar.gz (see README)") + _auto(m, "home"))]
    if checks[0]["ok"]:
        checks.append(_check("bundle environment (.venv)", bool(py), (py or "run `uv sync` then `./configure` inside the bundle") + _auto(m, "python")))
        miss = [t for t in ("tcsh", "gawk") if not shutil.which(t, path="/usr/bin:/bin")]
        checks.append(_check("tcsh and gawk", not miss, "ok" if not miss else "missing: sudo apt install " + " ".join(miss)))
        if len(home) > 51:
            checks.append(_check("path length", True, f"{len(home)} characters; NetCTL (not NetCTLpan) needs 51 or fewer"))
    if deep and all(c["ok"] for c in checks) and py:
        try:
            df, _ = ip.netctlpan([("t", "MSLLTEVETYVLSIIPSGPLKAEIAQRLEDVFAGKNTDLEVLMEWLKTRPILSPLTKGILGFVFTLTVPSERG")], [9], ["HLA-A*02:01"])
            hit = df[df["Peptide"] == "GILGFVFTL"]
            good = len(hit) == 1 and float(hit["CTLpan_%Rank"].iloc[0]) < 5
            checks.append(_check("NetCTLpan: GILGFVFTL / A*02:01", good,
                                 f"%Rank {float(hit['CTLpan_%Rank'].iloc[0]):.2f}" if len(hit) else "no result"))
        except Exception as e:
            checks.append(_check("NetCTLpan: GILGFVFTL / A*02:01", False, e))
    return checks


def health_mixmhc2pred(deep):
    m = config.model("mixmhc2pred")
    home = m.get("home") or ""
    checks = [_check("MixMHC2pred folder", os.path.isfile(os.path.join(home, "MixMHC2pred_unix")),
                     (home or f"not found under {config.search_root()} (clone GfellerLab/MixMHC2pred there)") + _auto(m, "home"))]
    if not checks[0]["ok"]:
        return checks
    pwm = os.path.join(home, "PWMdef")
    have = os.path.isdir(pwm) and any(f.endswith(".txt") for f in os.listdir(pwm))
    checks.append(_check("allele definitions (PWMdef)", have,
                         "ok" if have else mx2.HINT))
    if deep and have:
        try:
            import pandas as _pd
            df, _ = mx2.predict(_pd.DataFrame({"Protein": ["t"], "Start": [1], "Peptide": ["PKYVKQNTLKLAT"]}), ["HLA-DRB1*01:01"])
            good = len(df) == 1 and float(df["MIX2_%Rank"].iloc[0]) < 20
            checks.append(_check("PKYVKQNTLKLAT / DRB1*01:01 prediction", good,
                                 f"%Rank {float(df['MIX2_%Rank'].iloc[0]):.2f}" if len(df) else "no result"))
        except Exception as e:
            checks.append(_check("PKYVKQNTLKLAT / DRB1*01:01 prediction", False, e))
    return checks


def health_iedb_bundle_ii(deep):
    import shutil
    m = config.model("iedb_bundle_ii")
    home, py = m.get("home") or "", m.get("python") or ""
    checks = [_check("IEDB class II bundle folder", os.path.isfile(os.path.join(home, "src", "tcell_mhcii.py")),
                     (home or f"not found under {config.search_root()}: download IEDB_NG_TC2-*.tar.gz (see README)") + _auto(m, "home"))]
    if checks[0]["ok"]:
        checks.append(_check("bundle environment (.venv)", bool(py), (py or "run `uv sync` inside the bundle") + _auto(m, "python")))
        miss = [t for t in ("tcsh", "gawk", "perl", "R") if not shutil.which(t, path="/usr/bin:/bin")]
        checks.append(_check("tcsh, gawk, perl, R", not miss, "ok" if not miss else
                             "missing: " + " ".join(miss) + "  (sudo apt install tcsh gawk perl r-base-core)"))
    if deep and all(c["ok"] for c in checks) and py:
        import pandas as _pd
        rows = _pd.DataFrame({"Protein": ["t"], "Start": [1], "Peptide": ["PKYVKQNTLKLAT"]})
        try:
            df, _ = nn2.predict(rows, ["HLA-DRB1*01:01"], [13], None, kinds=("el",), version="4.3")
            good = len(df) == 1 and float(df["NII_%Rank_EL"].iloc[0]) < 2
            checks.append(_check("NetMHCIIpan: PKYVKQNTLKLAT / DRB1*01:01", good,
                                 f"EL %Rank {float(df['NII_%Rank_EL'].iloc[0]):.2f}" if len(df) else "no result"))
        except Exception as e:
            checks.append(_check("NetMHCIIpan: PKYVKQNTLKLAT / DRB1*01:01", False, e))
        try:
            e = cd4.predict(["PKYVKQNTLKLAT"])
            checks.append(_check("CD4episcore", len(e) == 1, f"score {float(e['CD4_Episcore'].iloc[0]):.2f}" if len(e) else "no result"))
        except Exception as ex:
            checks.append(_check("CD4episcore", False, ex))
    return checks


def health_mhcnuggets(deep):
    m = config.model("mhcnuggets")
    py = m.get("python") or ""
    checks = [_check("Python with mhcnuggets (TensorFlow)", bool(py) and os.path.isfile(py),
                     (py or "not found: " + mn2.HINT) + _auto(m, "python"))]
    if deep and checks[0]["ok"]:
        try:
            df, _ = mn2.predict(["PKYVKQNTLKLATGM"], ["HLA-DRB1*01:01"])
            checks.append(_check("MHCnuggets: PKYVKQNTLKLATGM / DRB1*01:01", len(df) == 1, f"IC50 {float(df['MHCnuggets_IC50_nM'].iloc[0]):.0f} nM" if len(df) else "no result"))
        except Exception as e:
            checks.append(_check("MHCnuggets: PKYVKQNTLKLATGM / DRB1*01:01", False, e))
    return checks


HEALTH = {"mhcnuggets": health_mhcnuggets, "netmhcpan": health_netmhcpan, "mixmhcpred": health_mixmhcpred, "mhcflurry": health_mhcflurry,
          "prime": health_prime, "iedb_bundle": health_iedb_bundle, "bigmhc": health_bigmhc,
          "mixmhc2pred": health_mixmhc2pred, "iedb_bundle_ii": health_iedb_bundle_ii}


@main.get("/api/config")
def get_config():
    try:
        config.load()
        err = None
    except config.ConfigError as e:
        err = str(e)
    eff = None
    if not err:
        eff = config.effective()
    return jsonify(path=config.CONFIG_PATH, text=config.raw_text(), error=err,
                   models=list(HEALTH), effective=eff)


@main.put("/api/config")
def put_config():
    try:
        config.save_text(request.get_json(force=True).get("text", ""))
    except config.ConfigError as e:
        return jsonify(error=str(e)), 400
    return jsonify(ok=True, path=config.CONFIG_PATH)


@main.post("/api/health/<name>")
def health(name):
    if name not in HEALTH:
        return jsonify(error="unknown model"), 404
    try:
        cfg = config.model(name)
        if not cfg.get("enabled", True):
            return jsonify(ok=False, checks=[_check("enabled", False, "disabled in settings")])
        checks = HEALTH[name](request.args.get("deep") == "1")
    except config.ConfigError as e:
        checks = [_check("settings", False, e)]
    except Exception as e:                                  # health must never 500 the page
        checks = [_check("unexpected error", False, e)]
    return jsonify(ok=all(c["ok"] for c in checks), checks=checks)


# ---- pipeline API ----------------------------------------------------------------------
_all_cache = {}


@main.get("/api/alleles")
def api_alleles():
    """Presets (canonical names) and every allele NetMHCpan knows, for the pipeline form."""
    key = netmhcpan_app.netmhc_path()
    if key not in _all_cache:
        seen = {allele_lib.canonical(a) for a in netmhcpan_app.list_alleles()}
        seen |= {a for v in allele_lib.PRESETS.values() for a in v["alleles"]}     # searchable even if NetMHCpan is missing
        _all_cache.clear()
        _all_cache[key] = sorted(seen)
    pre1, pre2 = dict(allele_lib.PRESETS), dict(allele_lib.PRESETS_II)
    pop, msg = _population()
    c1, c2, c3 = population.presets(pop) if pop else (None, None, None)
    if c1:
        pre1 = {f"{pop.label} population (common alleles)": c1, **pre1}
        pre2 = {f"{pop.label} population (common DRB1)": c2, f"{pop.label} population (DR + DQ)": c3, **pre2}
    return jsonify(presets=pre1, all=_all_cache[key], presets_ii=pre2)


def _population():
    return population.load(config.load().get("population"), config.search_root())


@main.get("/api/population")
def api_population():
    pop, msg = _population()
    return jsonify(ok=pop is not None, message=msg, info=pop.info() if pop else None)


def _needed_models(params):
    if _is_dual(params):
        out = []
        for c, p in pipeline_dual.split(params).items():
            out += _needed_models(p) if p else []
        return list(dict.fromkeys(out))
    m = params.get("models") or {}
    fun = params.get("funnel") or {}
    if fun.get("mode") == "manual":
        m = {k: True for st in fun.get("stages") or [] for k in st.get("models", [])}
    need = []
    if int(params.get("mhc_class", 1)) == 2:
        if m.get("netmhciipan_el") or m.get("netmhciipan_ba") or m.get("cd4episcore"): need.append("iedb_bundle_ii")
        if m.get("mixmhc2pred"): need.append("mixmhc2pred")
        if m.get("mhcnuggets"): need.append("mhcnuggets")
        return need
    if m.get("netmhcpan"): need.append("netmhcpan")
    if m.get("mixmhcpred") or m.get("prime"): need.append("mixmhcpred")      # PRIME calls MixMHCpred
    if m.get("prime"): need.append("prime")
    if m.get("mhcflurry"): need.append("mhcflurry")
    if m.get("netctlpan") or m.get("netchop"): need.append("iedb_bundle")
    if m.get("bigmhc_el") or m.get("bigmhc_im"): need.append("bigmhc")
    return need


def _is_dual(params):
    return "c1" in params or "c2" in params


def _pipe(params):
    if _is_dual(params):
        return pipeline_dual
    return pipeline_ii if int(params.get("mhc_class", 1)) == 2 else pipeline


@main.get("/api/alleles/search")
def api_allele_search():
    """Class II allele search (MixMHC2pred's table has ~350,000 names, so it is searched here, not in the browser)."""
    q = request.args.get("q", "").strip().lower().replace("*", "").replace(":", "").replace("_", "")
    if len(q) < 2:
        return jsonify(hits=[])
    pool = {a for v in allele_lib.PRESETS_II.values() for a in v["alleles"]}
    try:
        home = config.model("mixmhc2pred").get("home")
        if home and os.path.isdir(os.path.join(home, "PWMdef")):
            pool |= set(mx2.allele_map(home))
    except OSError:
        pass
    norm = lambda a: a.lower().replace("*", "").replace(":", "").replace("_", "")
    hits = sorted((a for a in pool if q in norm(a)), key=lambda a: (len(a), a))[:80]
    return jsonify(hits=hits)


@main.post("/api/pipeline/preview")
def pipeline_preview():
    try:
        params = request.get_json(force=True)
        return jsonify(_pipe(params).preview(params))
    except (ValueError, config.ConfigError) as e:
        return jsonify(hint=str(e))            # normal while the form is half filled: not an HTTP error


@main.post("/api/pipeline/run")
def pipeline_run():
    params = request.get_json(force=True)
    try:
        _pipe(params).normalise(params)
        for name in _needed_models(params):
            if not config.model(name).get("enabled", True):
                return jsonify(error=f"{name} is disabled in Settings."), 400
            bad = [c for c in HEALTH[name](False) if not c["ok"]]
            if bad:
                return jsonify(error=f"{name} is not ready: " + "; ".join(f"{c['name']} ({c['detail']})" for c in bad)), 400
    except (ValueError, config.ConfigError) as e:
        return jsonify(error=str(e)), 400
    job = _submit("pipeline", params, lambda j, mod=_pipe(params): mod.run(params, j))
    return jsonify(id=job.id)


def _submit(kind, params, fn):
    """Queue a job and keep a record of it in the workspace (input, settings, dates, outputs)."""
    prof = params.get("profile")

    def with_profile(j):
        config.set_profile(prof)                  # "max" / "balanced" / "gentle" for this run, else the one in Settings
        e = config.performance()
        j.log(f"performance profile {e['profile']}: {e['workers']} programs side by side, {e['threads']} threads per tool, "
              f"{e['shards']} BigMHC/MHCflurry processes ({e['cpus']} CPU threads, {e['mem_gb']} GB free)")
        try:
            return fn(j)
        finally:
            config.set_profile(None)
    return MANAGER.submit(with_profile, kind, on_create=lambda j: runs.create(j, kind, params), on_finish=runs.finish)


def _job(job_id):
    """A running/recent job, or an old one rebuilt from the workspace (also after a server restart)."""
    return MANAGER.get(job_id) or runs.load_job(job_id)


@main.post("/api/single/preview")
def single_preview():
    try:
        return jsonify(pipeline_single.preview(request.get_json(force=True)))
    except (ValueError, config.ConfigError) as e:
        return jsonify(hint=str(e))


@main.post("/api/single/run")
def single_run():
    params = request.get_json(force=True)
    try:
        p = pipeline_single.normalise(params)
        for name in pipeline_single.needed(p):
            if not config.model(name).get("enabled", True):
                return jsonify(error=f"{name} is disabled in Settings."), 400
            bad = [c for c in HEALTH[name](False) if not c["ok"]]
            if bad:
                return jsonify(error=f"{name} is not ready: " + "; ".join(f"{c['name']} ({c['detail']})" for c in bad)), 400
    except (ValueError, config.ConfigError) as e:
        return jsonify(error=str(e)), 400
    job = _submit("single", params, lambda j: pipeline_single.run(params, j))
    return jsonify(id=job.id)


@main.get("/api/runs")
def runs_list():
    live = {j: MANAGER.get(j) for j in [r["id"] for r in runs.listing()]}
    out = runs.listing()
    for r in out:                                       # a live job knows its progress better than the file
        j = live.get(r["id"])
        if j and j.state in ("queued", "running"):
            r["state"], r["progress"] = j.state, j.snapshot(False)["progress"]
    return jsonify(runs=out, workspace=config.workspace(), total_size=sum(r["size"] for r in out))


@main.post("/api/runs/<run_id>/rerank")
def runs_rerank(run_id):
    """The same predictions under another strategy: seconds, no model runs; the result is a new run in Projects."""
    body = request.get_json(force=True) or {}
    src, meta = _job(run_id), runs.read(run_id) if run_id.isalnum() else None
    if not src or not meta or src.state != "done":
        return jsonify(error="Only a finished run can be re-ranked."), 400
    preset = body.get("preset") or (meta["params"].get("preset") or "balanced")
    params = {**meta["params"], "preset": preset, "rules": body.get("rules") or {}, "rerank_of": run_id,
              "title": body.get("title") or f"{meta['name']} (re-ranked: {preset})"}
    try:
        for c in (1, 2):
            if params.get(f"c{c}"):
                import consensus
                consensus.resolve_rules(c, preset, params["rules"])
    except ValueError as e:
        return jsonify(error=str(e)), 400
    job = _submit("pipeline", params, lambda j: pipeline_dual.rerank(src, params, j))
    return jsonify(id=job.id)


@main.post("/api/runs/<run_id>/resume")
def runs_resume(run_id):
    """Run a stopped, failed or cancelled run again: the same settings; classes that already finished are reused, not recomputed."""
    meta = runs.read(run_id) if run_id.isalnum() else None
    if not meta:
        return jsonify(error="unknown run"), 404
    params = {**meta["params"], "reuse": True, "resumed_from": run_id}
    if meta["kind"] == "single":
        job = _submit("single", params, lambda j: pipeline_single.run(params, j))
        return jsonify(id=job.id)
    try:
        pipeline_dual.normalise(params)
        for name in _needed_models(params):
            bad = [c for c in HEALTH[name](False) if not c["ok"]]
            if bad:
                return jsonify(error=f"{name} is not ready: " + "; ".join(f"{c['name']} ({c['detail']})" for c in bad)), 400
    except (ValueError, config.ConfigError) as e:
        return jsonify(error=str(e)), 400
    job = _submit("pipeline", params, lambda j: pipeline_dual.run(params, j))
    return jsonify(id=job.id)


@main.get("/api/runs/<run_id>")
def runs_get(run_id):
    m = runs.read(run_id) if run_id.isalnum() else None
    if not m:
        return jsonify(error="unknown run"), 404
    m.pop("log", None)
    return jsonify(m)


@main.put("/api/runs/<run_id>")
def runs_put(run_id):
    body = request.get_json(force=True)
    m = runs.update(run_id, **{k: body[k] for k in ("name", "notes", "starred") if k in body}) if run_id.isalnum() else None
    return (jsonify(ok=True) if m else (jsonify(error="unknown run"), 404))


@main.delete("/api/runs/<run_id>")
def runs_delete(run_id):
    j = MANAGER.get(run_id)
    if j and j.state in ("queued", "running"):
        return jsonify(error="This run is still going: cancel it first."), 400
    try:
        runs.delete(run_id)
    except ValueError as e:
        return jsonify(error=str(e)), 400
    return jsonify(ok=True)


@main.get("/api/jobs/<job_id>")
def job_status(job_id):
    job = _job(job_id)
    return (jsonify(job.snapshot()) if job else (jsonify(error="unknown job (server restarted?)"), 404))


@main.post("/api/jobs/<job_id>/cancel")
def job_cancel(job_id):
    job = _job(job_id)
    if not job:
        return jsonify(error="unknown job"), 404
    job.cancel()
    return jsonify(ok=True)


def _res(job, cls=None):
    """The result dict to read tables from: the job's own, or one class of a combined job."""
    if not job or not job.result:
        return None
    r = job.result
    if r.get("dual"):
        return r["classes"].get(str(cls or min(r["classes"])))
    return r


def _table_path(job, table, cls=None):
    if job and job.result and table == "model":
        return job.result.get("permodel", {}).get(f"{cls}:{request.args.get('key', '')}")
    if job and job.result and table in ("regions", "report", "bundle"):
        return job.result.get({"regions": "regions_csv", "report": "report_html", "bundle": "bundle_zip"}[table])
    r = _res(job, cls)
    if not r:
        return None
    key = {"summary": "summary_csv", "pairs": "pairs_csv", "details": "details_json", "single": "table_csv"}.get(table, f"{table}_csv")
    return r.get(key)


@main.get("/api/jobs/<job_id>/result")
def job_result(job_id):
    path = _table_path(_job(job_id), request.args.get("table", "summary"), request.args.get("cls"))
    if not path:
        return jsonify(error="no result"), 404
    limit = min(int(request.args.get("limit", 5000)), 20000)
    df = pd.read_csv(path, nrows=limit)
    total = sum(1 for _ in open(path)) - 1
    return jsonify(columns=list(df.columns), rows=json.loads(df.to_json(orient="records")), total=total, shown=len(df))


@main.get("/api/jobs/<job_id>/csv")
def job_csv(job_id):
    table = request.args.get("table", "summary")
    path = _table_path(_job(job_id), table, request.args.get("cls"))
    if not path:
        return jsonify(error="no result"), 404
    if table == "single":
        return send_file(path, mimetype="text/csv", as_attachment=True, download_name=os.path.basename(path))
    if table == "details":
        return send_file(path, mimetype="application/json", as_attachment=True, download_name="epitope_run_details.json")
    if table == "report":
        return send_file(path, mimetype="text/html", as_attachment=request.args.get("inline") != "1", download_name="epitope_report.html")
    if table == "bundle":
        return send_file(path, mimetype="application/zip", as_attachment=True, download_name="epitope_all_outputs.zip")
    return send_file(path, mimetype="text/csv", as_attachment=True, download_name=f"epitopes_{table}.csv")


@main.get("/api/jobs/<job_id>/detail")
def job_detail(job_id):
    """Every tool value for one peptide window, per allele (the expandable row in the results)."""
    r = _res(_job(job_id), request.args.get("cls"))
    if not r:
        return jsonify(error="no result"), 404
    df = pd.read_csv(r["master_csv"], low_memory=False)
    df = df[(df["Peptide"] == request.args.get("peptide", "")) & (df["Start"] == int(request.args.get("start", 0)))]
    return jsonify(columns=list(df.columns), rows=json.loads(df.to_json(orient="records")))


@main.get("/api/funnel")
def api_funnel():
    cls = int(request.args.get("mhc_class", 1))
    cat = [{"key": k, "label": sp.label, "cost": sp.cost, "protein": sp.protein, "alleles": sp.alleles}
           for k, sp in model_registry.SPECS.items() if sp.cls == cls]
    return jsonify(models=sorted(cat, key=lambda x: (x["cost"], x["label"])), templates={k: {"label": v[0], "stages": v[1]} for k, v in funnel.templates(cls).items()})


@main.get("/api/cache")
def api_cache_get():
    return jsonify(size=predcache.size())


@main.delete("/api/cache")
def api_cache_clear():
    predcache.clear()
    return jsonify(ok=True)


def _profiles():
    out = {}
    for name in ("gentle", "balanced", "max"):
        config.set_profile(name)
        out[name] = config.performance()
    config.set_profile(None)
    return out


@main.get("/api/performance")
def api_performance_get():
    return jsonify(performance=config.load()["performance"], effective=config.performance(), profiles=_profiles())


@main.put("/api/performance")
def api_performance_put():
    body = request.get_json(force=True)
    try:
        cur = json.loads(config.raw_text())
        perf = {**config.DEFAULT["performance"], **cur.get("performance", {})}
        for k in perf:
            if k in body:
                perf[k] = str(body[k]) if k == "profile" else float(body[k]) if k == "timeout_minutes" else int(body[k])
        cur["performance"] = perf
        config.save_text(json.dumps(cur))
    except (ValueError, config.ConfigError) as e:
        return jsonify(error=str(e)), 400
    return jsonify(performance=perf, effective=config.performance())


@main.get("/api/rules")
def api_rules():
    import consensus
    cls = int(request.args.get("mhc_class", 1))
    return jsonify(presets={k: {"label": consensus.PRESET_INFO[k][0], "why": consensus.PRESET_INFO[k][1], "rules": v}
                            for k, v in consensus.preset_rules(cls).items()}, rank_by=list(consensus.RANK_BY))


@main.get("/api/jobs/<job_id>/report")
def job_report(job_id):
    job = _job(job_id)
    if not job or not job.result:
        return jsonify(error="no result"), 404
    r = job.result
    if r.get("dual"):
        out = {}
        for c, cr in r["classes"].items():
            with open(cr["details_json"]) as f:
                out[c] = {"consensus": cr["consensus"], "details": json.load(f)}
        return jsonify(dual=True, classes=out, summary=r["summary"])
    with open(r["details_json"]) as f:
        details = json.load(f)
    return jsonify(dual=False, classes={"1" if details.get("mhc_class", 1) == 1 else "2": {"consensus": r["consensus"], "details": details}}, summary=r["summary"])


@main.get("/")
def index():
    with open(os.path.join(HERE, "web", "index.html")) as f:        # re-read so UI edits show on refresh
        return f.read()


runs.mark_interrupted()
application = DispatcherMiddleware(main, {"/netmhcpan": netmhcpan_app.app, "/bigmhc": bigmhc_app.app})

if __name__ == "__main__":
    run_simple("127.0.0.1", int(os.environ.get("PORT", 5000)), application, threaded=True)
