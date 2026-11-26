"""Immune simulation: the engine's dynamic-response module.

The engine ranks epitopes. This module takes those epitopes - or a protein
sequence, or a finished engine run - and simulates what the immune system
actually does with them over months: T cells, B cells, class switching,
affinity maturation, memory, and the secondary response to a boost.

It lives inside the engine and behaves like the rest of it: the same job
manager (one job at a time), the same stage-by-stage progress, the same
health-check cards, and an Install / Update pair in Settings. Install trains
the epitope models from the datasets committed under immune_sim/data, and
needs no network. Update re-downloads those datasets from IEDB, CMI-PB and
ClinicalTrials.gov.

Everything long runs as the project's own CLI script through run_cmd, so the
GUI cannot drift from what VALIDATION.md measured.
"""
import json
import os
import sys

from flask import jsonify, request, send_file

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, "immune_sim")
if SIM not in sys.path:
    sys.path.insert(0, SIM)

from common.jobs import MANAGER, run_cmd              # noqa: E402

PY_EXE = sys.executable or "python3"


# ---------------------------------------------------------------------------
# what is present
# ---------------------------------------------------------------------------
def _p(*parts):
    return os.path.join(SIM, *parts)


def installed():
    """True when the epitope models have been trained at least once."""
    return os.path.exists(_p("models", "report.json"))


def _dataio():
    sys.path.insert(0, SIM)
    import dataio
    return dataio


def datasets():
    try:
        return _dataio().inventory()
    except Exception as e:                                   # noqa: BLE001
        return [{"name": "error", "what": str(e), "present": False}]


def _report(name):
    f = _p("models", name)
    if not os.path.exists(f):
        return None
    try:
        with open(f) as fh:
            return json.load(fh)
    except (ValueError, OSError):
        return None


# The Simulation module runs in THIS interpreter, not in a predictor's own
# environment, so it needs its packages here. An environment created before
# the simulator was merged in has only flask and pandas, and the first thing
# Install did was crash with "No module named 'scipy'" several steps in.
# Check up front instead, and offer to fix it.
NEEDED = [("numpy", "numpy"), ("scipy", "scipy"),
          ("sklearn", "scikit-learn"), ("pandas", "pandas")]


def missing_packages():
    """[(import name, pip name)] that this interpreter cannot import."""
    import importlib.util
    return [(mod, pip) for mod, pip in NEEDED
            if importlib.util.find_spec(mod) is None]


def health(deep=False):
    """Health card lines, in the engine's own format."""
    from server import _check                          # late: avoids a cycle
    out = []
    miss = missing_packages()
    out.append(_check(
        "python packages", not miss,
        "all present in " + os.path.basename(PY_EXE) if not miss
        else ("missing: " + ", ".join(p for _m, p in miss)
              + " - press Install, or run: "
              + f"{PY_EXE} -m pip install " + " ".join(p for _m, p in miss))))
    data = datasets()
    missing = [d["name"] for d in data if not d.get("present")]
    bundled = sum(1 for d in data if d.get("copy") == "bundled")
    out.append(_check(
        "datasets", not missing,
        f"{len(data) - len(missing)} of {len(data)} present"
        + (f", {bundled} from the copy committed in this repository" if bundled else "")
        + (f"; missing: {', '.join(missing)}" if missing else "")))

    rep = _report("report.json")
    out.append(_check("per-allele binding models", bool(rep),
                      f"{len(rep)} alleles trained" if rep
                      else "not trained yet - press Install"))

    pan = _report("pan_report.json")
    detail = "not trained - the per-allele models are used instead"
    if pan:
        bits = []
        for cls, r in pan.items():
            auc = (r.get("metrics", {}).get("ensemble", {}) or {}).get("auc")
            if auc is not None:
                bits.append(f"class {cls} AUC {auc}")
        detail = ("; ".join(bits) or "trained") + " (random split: optimistic)"
    out.append(_check("pan-allele core-scan model", bool(pan), detail, warn=True))

    cal = os.path.exists(_p("models", "calibrated_wetlab.json"))
    out.append(_check("calibrated dynamic model", cal,
                      "models/calibrated_wetlab.json - validated against "
                      "held-out human subjects" if cal
                      else "missing; the simulation would use unvalidated defaults"))

    wl = os.path.exists(_p("results", "wetlab_validation_test.json"))
    out.append(_check("held-out human check has been run", wl,
                      "results/wetlab_validation_test.json" if wl
                      else "not run here yet - Install ends with it", warn=True))
    return out


# ---------------------------------------------------------------------------
# jobs: install and update
# ---------------------------------------------------------------------------
INSTALL_STEPS = [
    ("train_allele", "Per-allele binding models",
     [PY_EXE, "-u", "model.py"]),
    ("train_pan", "Pan-allele core-scan ensemble",
     [PY_EXE, "-u", "-c", "import model_pan; model_pan.main()", "--realign", "4"]),
    ("train_immuno", "Immunogenicity model",
     [PY_EXE, "-u", "immunogenicity.py"]),
    ("wetlab", "Check against held-out human subjects",
     [PY_EXE, "-u", "wetlab_validation.py", "--split", "test",
      "--params", "models/calibrated_wetlab.json"]),
    ("benchmarks", "13 platform benchmarks",
     [PY_EXE, "-u", "test_platforms.py"]),
]

QUICK_STEPS = [INSTALL_STEPS[0], INSTALL_STEPS[4]]

UPDATE_STEPS = [
    ("harvest_iedb", "IEDB binding + immunogenicity", [PY_EXE, "-u", "harvest_iedb.py"]),
    ("harvest_cmipb", "CMI-PB human boost data", [PY_EXE, "-u", "harvest_cmipb.py"]),
    ("harvest_eluted", "Mass-spectrometry eluted ligands", [PY_EXE, "-u", "harvest_eluted.py"]),
    ("harvest_trials", "ClinicalTrials.gov fold-rises", [PY_EXE, "-u", "harvest_trials.py"]),
]

SENSITIVITY_STEPS = [
    ("sensitivity", "Controls, scrambles and point mutations",
     [PY_EXE, "-u", "sensitivity.py", "--reps", "10"]),
]

VALIDATE_STEPS = [
    ("benchmarks", "13 platform benchmarks", [PY_EXE, "-u", "test_platforms.py"]),
    ("wetlab", "Held-out human subjects",
     [PY_EXE, "-u", "wetlab_validation.py", "--split", "test",
      "--params", "models/calibrated_wetlab.json"]),
    ("lopo", "Leave-one-pathogen-out, per-allele", [PY_EXE, "-u", "external_validation.py"]),
    ("lopo_pan", "Leave-one-pathogen-out, pan-allele",
     [PY_EXE, "-u", "external_validation_pan.py", "--realign", "4"]),
    ("eluted", "Mass-spectrometry eluted ligands", [PY_EXE, "-u", "validate_eluted.py"]),
    ("trials", "Clinical trial fold-rises", [PY_EXE, "-u", "validate_trials.py"]),
    ("scorecard", "Evidence scorecard + robustness",
     [PY_EXE, "-u", "scorecard.py", "--robustness", "100"]),
    ("sensitivity", "Controls, scrambles and point mutations",
     [PY_EXE, "-u", "sensitivity.py", "--reps", "10"]),
]

PLANS = {
    "install": dict(steps=INSTALL_STEPS, label="Install",
                    about="Trains all three epitope models from the datasets "
                          "committed in this repository, then checks the "
                          "dynamic model against held-out human subjects and "
                          "runs the 13 platform benchmarks. No network needed.",
                    est="about 10 minutes"),
    "quick": dict(steps=QUICK_STEPS, label="Quick install",
                  about="Only the fast per-allele models and the benchmarks, "
                        "so the Simulation tab works straight away. The "
                        "pan-allele model is the better predictor on an "
                        "unseen pathogen (0.873 against 0.766 class I); run "
                        "the full install when you have ten minutes.",
                  est="about 1 minute"),
    "update": dict(steps=UPDATE_STEPS, label="Update data",
                   about="Re-downloads every dataset from IEDB, CMI-PB and "
                         "ClinicalTrials.gov, replacing the committed copies. "
                         "Every number in VALIDATION.md came from the "
                         "committed copies, so results may move after this.",
                   est="15-40 minutes, needs the network"),
    "sensitivity": dict(
        steps=SENSITIVITY_STEPS, label="Check antigen sensitivity",
        about="Runs the simulation against positive controls, self proteins, "
              "low-complexity sequences, shuffled versions of each antigen, "
              "and single point substitutions at anchor and non-anchor "
              "positions. This is the standing check that the model reads "
              "the antigen rather than returning the same response for "
              "anything - the defect reported for C-ImmSim.",
        est="1-2 hours (86 arms, each a full protein scan)"),
    "validate": dict(steps=VALIDATE_STEPS, label="Re-run the validation",
                     about="Everything in VALIDATION.md end to end. Any of "
                           "it can fail, and a failure is a result.",
                     est="about 2 hours"),
}


def _install_packages(job):
    """Install the Simulation module's dependencies into this interpreter.

    Only ever installs what is actually missing, and names them in the log.
    This is what "Install" means for a module that runs in the server's own
    process, and it is shown as a step rather than done silently.
    """
    miss = missing_packages()
    if not miss:
        job.log("all required packages are already present")
        return
    pkgs = [pip for _m, pip in miss]
    job.log(f"installing into {PY_EXE}: {', '.join(pkgs)}")
    r = run_cmd([PY_EXE, "-m", "pip", "install", *pkgs], job=job, merge=True,
                on_line=lambda l: job.log(l.rstrip()), cwd=SIM)
    if r.returncode != 0:
        raise RuntimeError(
            "Could not install " + ", ".join(pkgs) + ". Install them yourself "
            "and press Install again:  " + PY_EXE + " -m pip install "
            + " ".join(pkgs))
    still = missing_packages()
    if still:
        raise RuntimeError(
            "pip reported success but " + ", ".join(p for _m, p in still)
            + " still cannot be imported. If this is a conda environment, "
            "install them with conda instead.")


def _run_plan(job, steps):
    """Run the steps in order as one engine job, one stage each."""
    steps = [("deps", "Python packages the simulator needs", None)] + list(steps)
    for name, label, _cmd in steps:
        job.stage(name, label, weight=1.0)
    for name, label, cmd in steps:
        job.check()
        if cmd is None:                      # the dependency step
            job.update(name, state="running", msg="")
            job.log(f"=== {label} ===")
            _install_packages(job)
            job.finish_stage(name)
            continue
        job.update(name, state="running", msg="")
        job.log(f"=== {label} ===")

        def on_line(line, _n=name):
            line = line.rstrip()
            if line:
                job.log(line)
                job.update(_n, msg=line[:160])

        r = run_cmd(cmd, job=job, merge=True, on_line=on_line, cwd=SIM)
        if r.returncode != 0:
            job.update(name, state="failed",
                       msg=f"exited with code {r.returncode}")
            raise RuntimeError(
                f"{label} failed (exit code {r.returncode}). The log above is "
                f"the script's own output.")
        job.finish_stage(name)
    return {"ok": True, "steps": [s[0] for s in steps]}


# ---------------------------------------------------------------------------
# routes, registered onto the engine's own Flask app
# ---------------------------------------------------------------------------
def register(app):
    @app.get("/vendor/<path:name>")
    def sim_vendor(name):
        """Serve the charting library from disk, never from a CDN.

        The engine is run on machines with no internet - that is the whole
        reason its predictors and this module's datasets are local - so a
        CDN-loaded chart library would silently draw nothing.
        """
        safe = os.path.basename(name)
        path = os.path.join(HERE, "web", "vendor", safe)
        if not os.path.exists(path):
            return jsonify(error="not found"), 404
        kind = ("text/javascript" if safe.endswith(".js")
                else "text/markdown" if safe.endswith(".md") else "text/plain")
        return send_file(path, mimetype=kind)

    @app.get("/api/sim/runs")
    def sim_runs():
        """Finished engine runs that have epitopes worth simulating."""
        from common import runs
        out = []
        for r in runs.listing():
            if r.get("state") != "done":
                continue
            # runs.listing() gives `input` as a dict (proteins, residues,
            # mode, sha1). Sending it raw put "[object Object]" in the
            # dropdown, so the readable line is built here, once.
            inp = r.get("input") or {}
            bits = []
            prots = inp.get("proteins") or []
            if prots:
                bits.append(", ".join(prots[:2])
                            + (f" +{len(prots) - 2} more" if len(prots) > 2 else ""))
            elif inp.get("n_proteins"):
                bits.append(f"{inp['n_proteins']} proteins")
            if inp.get("residues"):
                bits.append(f"{inp['residues']:,} aa")
            out.append({"id": r["id"], "name": r.get("name") or r["id"][:8],
                        "kind": r.get("kind"), "created": r.get("created"),
                        "detail": " \u00b7 ".join(bits)})
        return jsonify(runs=out[:60])

    @app.get("/api/sim/status")
    def sim_status():
        rep = _report("report.json")
        return jsonify({
            "installed": installed(),
            "datasets": datasets(),
            "alleles": [{"allele": a, "mhc_class": v.get("mhc_class"),
                         "auc": v.get("auc")}
                        for a, v in sorted((rep or {}).items())],
            "plans": {k: {"name": k, "label": v["label"], "about": v["about"],
                          "est": v["est"],
                          "steps": [{"name": n, "label": l}
                                    for n, l, _c in v["steps"]]}
                      for k, v in PLANS.items()},
            "health": health(),
        })

    @app.post("/api/sim/run/<plan>")
    def sim_run(plan):
        spec = PLANS.get(plan)
        if spec is None:
            return jsonify(error=f"unknown plan {plan!r}"), 400
        if plan != "update" and not os.path.exists(_p("data")):
            return jsonify(error="immune_sim/data is missing"), 400
        job = MANAGER.submit(lambda j, s=spec["steps"]: _run_plan(j, s),
                             kind=f"sim_{plan}")
        return jsonify(id=job.id)

    @app.get("/api/sim/meta")
    def sim_meta():
        import platforms as PF
        return jsonify({
            "platforms": [{"key": k, "label": p.label, "notes": p.notes}
                          for k, p in PF.PLATFORMS.items()],
            "adjuvants": [{"key": k, "label": a.label, "notes": a.notes}
                          for k, a in PF.ADJUVANTS.items()],
        })

    @app.post("/api/sim/design")
    def sim_design():
        """A sequence, an engine run, or the demo antigen -> a response."""
        from immune_sim_design import design
        try:
            return jsonify(design(request.get_json(force=True) or {}))
        except ValueError as e:
            return jsonify(error=str(e)), 400
        except SystemExit as e:
            # SystemExit inherits from BaseException, so "except Exception"
            # does not catch it: it unwinds the worker thread and the browser
            # sees nothing but "Failed to fetch". These modules are libraries
            # here and CLIs elsewhere, so this stays as a backstop even though
            # the ones we call now raise ValueError.
            return jsonify(error=f"the simulator stopped: {e}"), 400
        except Exception as e:                               # noqa: BLE001
            import traceback
            if os.environ.get("IMMUNOFORGE_DEBUG"):
                traceback.print_exc()
            return jsonify(error=f"{type(e).__name__}: {e}"), 500

    @app.post("/api/sim/compare")
    def sim_compare():
        """Several formulations, antigens or schedules on one axis.

        Each arm goes through the same code path as a single simulation, so
        a curve here and a curve on the Simulation tab are the same
        computation. An arm that fails returns its reason instead of
        removing itself: a self protein giving nothing is a RESULT, and
        dropping it silently would hide the most informative arm in the
        comparison.
        """
        from immune_sim_design import design
        body = request.get_json(force=True) or {}
        arms = body.get("arms") or []
        if not 1 <= len(arms) <= 8:
            return jsonify(error="give between 1 and 8 arms"), 400
        shared = {k: v for k, v in body.items() if k != "arms"}
        out = []
        for i, a in enumerate(arms):
            cfg = dict(shared, **a)
            label = a.get("label") or f"arm {i + 1}"
            try:
                d = design(cfg)
            except ValueError as e:
                out.append({"label": label, "error": str(e),
                            "no_response": "tolerance" in str(e).lower()})
                continue
            except Exception as e:                           # noqa: BLE001
                out.append({"label": label, "error": f"{type(e).__name__}: {e}"})
                continue
            keep = ("IgG", "IgM", "TC_active", "TH_active", "B_active",
                    "B_memory", "TC_memory", "affinity", "antigen")
            out.append({
                "label": label,
                "antigen": d["antigen"],
                "day": d["day"],
                "series": {k: v for k, v in d["series"].items() if k in keep},
                "summary": d["summary"],
                "epitopes": {k: d["epitopes"][k]
                             for k in ("n", "n_alleles", "n_clones", "by_class")
                             if k in d["epitopes"]},
                "info": d["info"],
            })
        return jsonify(arms=out)

    @app.get("/api/sim/sensitivity")
    def sim_sensitivity():
        """The standing controls: positives, negatives, scrambles, mutations."""
        f = os.path.join(SIM, "results", "sensitivity.json")
        if not os.path.exists(f):
            return jsonify(available=False,
                           hint="Run it with: python3 immune_sim/sensitivity.py")
        with open(f) as fh:
            return jsonify(available=True, **json.load(fh))

    @app.get("/api/sim/validation")
    def sim_validation():
        sys.path.insert(0, os.path.join(SIM, "webapp"))
        import verdicts
        vs = verdicts.all_verdicts()
        return jsonify(verdicts=vs, headline=verdicts.headline(vs),
                       charts=verdicts.charts())
