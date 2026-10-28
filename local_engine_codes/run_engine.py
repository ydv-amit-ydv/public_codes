#!/usr/bin/env python3
"""Run the engine without the web page: python3 run_engine.py settings.json [-o outdir]

settings.json is exactly what the page sends (use "Save settings" on the page to get one): {"text": ..., "mode": "scan",
"preset": "balanced", "rules": {}, "c1": {...}, "c2": {...}}. The results are written to outdir (default ./epitope_run):
report.html, regions.csv, the per-class tables and a zip of everything. Needs the same models as the page; progress is printed.
"""
import argparse, json, os, shutil, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import engine_server                                        # noqa: E402,F401  (puts the tool folders on sys.path)
import pipeline_dual                                        # noqa: E402
from common import config                                   # noqa: E402
from common.jobs import Job                                 # noqa: E402
from common import runs                                     # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("settings")
    ap.add_argument("-o", "--out", default="epitope_run")
    ap.add_argument("--fasta", help="protein FASTA to use instead of the text in settings.json")
    a = ap.parse_args()
    params = json.load(open(a.settings))
    if a.fasta:
        params["text"], params["mode"] = open(a.fasta).read(), "scan"
    problems = []
    for name in engine_server._needed_models(params):
        bad = [c for c in engine_server.HEALTH[name](False) if not c["ok"]]
        if bad:
            problems.append(f"{name} is not ready: " + "; ".join(f"{c['name']} ({c['detail']})" for c in bad))
    if problems:
        sys.exit("\n".join(problems))
    job = Job()
    runs.create(job, "pipeline", params)
    job.state, job.started = "running", time.time()
    last = [-1]

    def watch():
        while job.state == "running":
            f = job.fraction()
            if int(f * 100) != last[0]:
                last[0] = int(f * 100)
                cur = next((s["label"] for s in job.stages if s["state"] == "running"), "")
                print(f"\r{last[0]:3d}%  {cur[:70]:70s}", end="", file=sys.stderr, flush=True)
            time.sleep(0.5)
    import threading
    threading.Thread(target=watch, daemon=True).start()
    try:
        pipeline_dual.run(params, job)
        job.state = "done"
    except Exception as e:                                   # noqa: BLE001
        job.state, job.error, job.finished = "failed", str(e), time.time()
        runs.finish(job)
        sys.exit(f"\nfailed: {type(e).__name__}: {e}")
    job.finished = time.time()
    runs.finish(job)                                         # the run is kept in the workspace too, so the page lists it
    os.makedirs(a.out, exist_ok=True)
    r = job.result
    shutil.copy(r["report_html"], os.path.join(a.out, "report.html"))
    shutil.copy(r["regions_csv"], os.path.join(a.out, "regions.csv"))
    shutil.copy(r["bundle_zip"], os.path.join(a.out, "all_outputs.zip"))
    print(f"\n{r['summary']['regions']} regions ({r['summary']['regions_dual']} CTL+HTL) in {r['summary']['seconds']} s -> {a.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
