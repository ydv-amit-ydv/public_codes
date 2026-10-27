"""Epitope Filter Engine server: one page with a tab per model plus a Settings tab
whose JSON (engine_config.json) can be edited live.

    python3 engine_server.py          # http://127.0.0.1:5000

Each model keeps its own interpreter/install (see common/config.py); this process
only needs flask and pandas. Settings are re-read on every request, so saving in
the GUI takes effect immediately.
"""
import os, subprocess, sys

from flask import Flask, jsonify, request, send_file
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.serving import run_simple

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ("MHC1/BA_EL_scripts", "MHC1/Immunogenicity_scripts"):
    sys.path.insert(0, os.path.join(HERE, sub))
sys.path.insert(0, HERE)

from common import alleles as allele_lib, config            # noqa: E402
from common.jobs import MANAGER                              # noqa: E402
import netmhcpan_app, bigmhc_app, bigmhc_runner as bm       # noqa: E402
import mixmhcpred_runner as mx, prime_runner as pr           # noqa: E402
import pipeline                                              # noqa: E402
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


HEALTH = {"netmhcpan": health_netmhcpan, "mixmhcpred": health_mixmhcpred, "prime": health_prime, "bigmhc": health_bigmhc}


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
    return jsonify(presets=allele_lib.PRESETS, all=_all_cache[key])


def _needed_models(params):
    m = params.get("models") or {}
    need = []
    if m.get("netmhcpan"): need.append("netmhcpan")
    if m.get("mixmhcpred") or m.get("prime"): need.append("mixmhcpred")      # PRIME calls MixMHCpred
    if m.get("prime"): need.append("prime")
    if m.get("bigmhc_el") or m.get("bigmhc_im"): need.append("bigmhc")
    return need


@main.post("/api/pipeline/preview")
def pipeline_preview():
    try:
        return jsonify(pipeline.preview(request.get_json(force=True)))
    except (ValueError, config.ConfigError) as e:
        return jsonify(hint=str(e))            # normal while the form is half filled: not an HTTP error


@main.post("/api/pipeline/run")
def pipeline_run():
    params = request.get_json(force=True)
    try:
        pipeline.normalise(params)
        for name in _needed_models(params):
            if not config.model(name).get("enabled", True):
                return jsonify(error=f"{name} is disabled in Settings."), 400
            bad = [c for c in HEALTH[name](False) if not c["ok"]]
            if bad:
                return jsonify(error=f"{name} is not ready: " + "; ".join(f"{c['name']} ({c['detail']})" for c in bad)), 400
    except (ValueError, config.ConfigError) as e:
        return jsonify(error=str(e)), 400
    job = MANAGER.submit(lambda j: pipeline.run(params, j))
    return jsonify(id=job.id)


@main.get("/api/jobs/<job_id>")
def job_status(job_id):
    job = MANAGER.get(job_id)
    return (jsonify(job.snapshot()) if job else (jsonify(error="unknown job (server restarted?)"), 404))


@main.post("/api/jobs/<job_id>/cancel")
def job_cancel(job_id):
    job = MANAGER.get(job_id)
    if not job:
        return jsonify(error="unknown job"), 404
    job.cancel()
    return jsonify(ok=True)


def _table_path(job, table):
    if not job or not job.result:
        return None
    return job.result["summary_csv" if table == "summary" else "pairs_csv"]


@main.get("/api/jobs/<job_id>/result")
def job_result(job_id):
    path = _table_path(MANAGER.get(job_id), request.args.get("table", "summary"))
    if not path:
        return jsonify(error="no result"), 404
    limit = min(int(request.args.get("limit", 5000)), 20000)
    df = pd.read_csv(path, nrows=limit)
    total = sum(1 for _ in open(path)) - 1
    df = df.astype(object).where(df.notna(), None)
    return jsonify(columns=list(df.columns), rows=df.to_dict("records"), total=total, shown=len(df))


@main.get("/api/jobs/<job_id>/csv")
def job_csv(job_id):
    table = request.args.get("table", "summary")
    path = _table_path(MANAGER.get(job_id), table)
    if not path:
        return jsonify(error="no result"), 404
    return send_file(path, mimetype="text/csv", as_attachment=True, download_name=f"epitopes_{table}.csv")


@main.get("/")
def index():
    with open(os.path.join(HERE, "web", "index.html")) as f:        # re-read so UI edits show on refresh
        return f.read()


application = DispatcherMiddleware(main, {"/netmhcpan": netmhcpan_app.app, "/bigmhc": bigmhc_app.app})

if __name__ == "__main__":
    run_simple("127.0.0.1", int(os.environ.get("PORT", 5000)), application, threaded=True)
