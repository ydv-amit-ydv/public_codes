#!/usr/bin/env python3
"""
Flask front end for immune_sim.

    cd immune_sim/webapp && python3 app.py
    open http://127.0.0.1:5000

What it does:

  /            platform + adjuvant + schedule -> simulated T/B/antibody course
  /epitopes    paste or fetch a protein -> ranked epitopes -> feed the sim
  /compare     several formulations on one axis
  /data        download every dataset; see what is present and what is stale
  /train       train the epitope models, with their held-out metrics
  /calibrate   refit the ODE parameters against measured human data
  /validation  the evidence behind every claim, and re-run any of it

Deliberately a thin layer: every route calls the same functions the CLI does,
and every long operation runs the same script the README documents, as a
subprocess with its output streamed to the browser. The GUI therefore cannot
drift from the validated pipeline, and a number on screen is always something
a reader could reproduce from a terminal.

Deliberately a thin layer: every route calls the same functions the CLI does,
so the GUI cannot drift from the validated pipeline. Nothing is cached that
would let a figure outlive the parameters it came from.
"""
import io
import json
import sys
import traceback
from pathlib import Path

import numpy as np
import pandas as pd
from flask import (Flask, abort, jsonify, render_template, request,
                   send_file, Response)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT))

import platforms as PF                                      # noqa: E402
from dynamics import P, Sim, load_epitopes                   # noqa: E402

sys.path.insert(0, str(HERE))
import jobs as JOBS                                          # noqa: E402
import tasks as TASKS                                        # noqa: E402
import verdicts as VERDICTS                                  # noqa: E402
sys.path.insert(0, str(ROOT))
import dataio as DATAIO                                      # noqa: E402

MANAGER = JOBS.JobManager()

RESULTS = ROOT / "results"
MODELS = ROOT / "models"
DEMO_EPITOPES = RESULTS / "demo_epitopes.csv"

app = Flask(__name__)
app.config["JSON_SORT_KEYS"] = False

# Series the browser plots. Anything else in the frame is still downloadable
# as CSV, but sending 60 columns of per-clone detail to a chart is noise.
SERIES = [
    ("antigen", "Antigen (free)", "antigen"),
    ("antigen_depot", "Antigen on FDC", "antigen"),
    ("template", "Template (mRNA/vector)", "antigen"),
    ("IgM", "IgM", "antibody"),
    ("IgG", "IgG", "antibody"),
    ("affinity", "Antibody affinity", "affinity"),
    ("TH_active", "CD4 T effector", "tcell"),
    ("TH_memory", "CD4 T memory", "tcell"),
    ("TC_active", "CD8 T effector (CTL)", "tcell"),
    ("TC_memory", "CD8 T memory", "tcell"),
    ("B_active", "B cell (activated)", "bcell"),
    ("B_memory", "B cell (memory)", "bcell"),
    ("plasma", "Plasma cells (short-lived)", "bcell"),
    ("plasma_longlived", "Plasma cells (long-lived)", "bcell"),
    ("danger", "Danger signal", "innate"),
    ("DC_mature", "Mature DC", "innate"),
    ("IL2", "IL-2", "innate"),
]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def _params(which="wetlab"):
    """
    Which calibration to run with. `wetlab` is the set validated against
    CMI-PB and is the default everywhere; `literature` is kept so the GUI can
    show what changed, and `default` is the unfitted starting point.
    """
    f = {"wetlab": MODELS / "calibrated_wetlab.json",
         "literature": MODELS / "calibrated.json"}.get(which)
    if f and f.exists():
        return json.loads(f.read_text())["params"], f.name
    return dict(P), "defaults (uncalibrated)"


def _epitope_table(payload):
    """Epitopes for a run: an uploaded CSV, a prior prediction, or the demo."""
    src = (payload.get("epitope_source") or "demo").strip()
    if src == "demo" or not src:
        if not DEMO_EPITOPES.exists():
            abort(503, "results/demo_epitopes.csv is missing; run predict.py")
        return load_epitopes(DEMO_EPITOPES), "demo (SARS-CoV-2 spike)"
    if src == "inline":
        rows = payload.get("epitope_rows") or []
        if not rows:
            abort(400, "no epitopes were given for epitope_source=inline")
        d = pd.DataFrame(rows)
        return load_epitopes(d), "pasted sequence"
    f = RESULTS / Path(src).name                    # never leave results/
    if not f.exists():
        abort(400, f"no such epitope file: {f.name}")
    return load_epitopes(f), f.name


def _schedule(payload):
    """(step, dose, adjuvant) tuples from the form's day numbers."""
    adj = payload.get("adjuvant") or None
    days = payload.get("inject_days") or [0, 90]
    out = []
    for d in days:
        step = max(2, int(round(float(d) * 24.0 / 8.0))) if float(d) > 0 else 2
        out.append((step, None, adj))
    return sorted(set(out)) or [(2, None, adj)]


def run_sim(payload):
    ep, ep_name = _epitope_table(payload)
    params, params_name = _params(payload.get("calibration", "wetlab"))
    if payload.get("replicating"):
        params = dict(params, ag_growth=0.55)   # same value simulate.py uses
    inj = _schedule(payload)
    days = float(payload.get("days", 240))
    n_steps = int(max(30, min(days * 3, 3000)))      # 8 h per step

    sim = Sim(ep, p=params, injections=inj,
              platform=payload.get("platform", "protein"),
              adjuvant=payload.get("adjuvant") or None,
              dose_ug=(float(payload["dose_ug"])
                       if payload.get("dose_ug") else None),
              heterologous=payload.get("heterologous") or None)
    res_h = float(payload.get("resolution_h", 8.0))
    df = sim.run(n_steps=n_steps, out_hours=max(0.5, min(res_h, 48.0)))
    return df, sim, dict(epitopes=ep_name, params=params_name,
                         n_clones=int(len(ep)),
                         injection_days=[round(s * 8 / 24, 2) for s, _, _ in inj])


def _summary(df, sim, inj_days):
    """The handful of numbers a user actually reads off a run."""
    boost = max(inj_days) if len(inj_days) > 1 else None
    out = {
        "peak_IgG": float(df.IgG.max()),
        "peak_IgG_day": float(df.day[df.IgG.idxmax()]),
        "peak_IgM": float(df.IgM.max()),
        "peak_CD4": float(df.TH_active.max()),
        "peak_CTL": float(df.TC_active.max()),
        "CD4_memory_end": float(df.TH_memory.iloc[-1]),
        "CTL_memory_end": float(df.TC_memory.iloc[-1]),
        "B_memory_end": float(df.B_memory.iloc[-1]),
        "peak_B": float(df.B_active.max()),
        "LLPC_end": float(df.get("plasma_longlived", df.plasma).iloc[-1]),
        "igg_durability": float(df.IgG.iloc[-1] / max(df.IgG.max(), 1e-9)),
        "affinity_fold": float(df.affinity.max() / max(df.affinity.iloc[0], 1e-9)),
        "IgG_over_IgM": float(df.IgG.max() / max(df.IgM.max(), 1e-9)),
        "CTL_over_IgG": float(df.TC_active.max() / max(df.IgG.max(), 1e-9)),
    }
    if boost:
        pre, post = df[df.day < boost], df[df.day >= boost]
        if not pre.empty and not post.empty:
            base = float(pre.IgG.iloc[-1])
            out["boost_fold_over_baseline"] = float(post.IgG.max()
                                                    / max(base, 1e-9))
            out["boost_peak_day_after"] = float(
                post.day[post.IgG.idxmax()] - boost)
            out["primary_peak"] = float(pre.IgG.max())
            out["secondary_peak"] = float(post.IgG.max())
    pf = sim.pf
    out["platform_label"] = pf.label
    out["adjuvant_label"] = sim.adj.label
    out["dose_ug"] = sim.dose_ug
    return out


# --------------------------------------------------------------------------
# pages
# --------------------------------------------------------------------------
@app.route("/")
def page_design():
    """The merged pipeline: antigen in, predicted response out."""
    return render_template("design.html", meta=_meta())


@app.route("/simulate")
def page_simulate():
    """The original epitope-file-driven simulator, kept for the CSV route."""
    return render_template("simulate.html", meta=_meta())


@app.route("/epitopes")
def page_epitopes():
    return render_template("epitopes.html", meta=_meta())


@app.route("/compare")
def page_compare():
    return render_template("compare.html", meta=_meta())


@app.route("/setup")
def page_setup():
    return render_template("setup.html", meta=_meta(), group="setup",
                           title="Setup", groups=TASKS.GROUPS)


@app.route("/data")
def page_data():
    return render_template("run.html", meta=_meta(), group="data",
                           title="Data", groups=TASKS.GROUPS)


@app.route("/train")
def page_train():
    return render_template("run.html", meta=_meta(), group="train",
                           title="Models", groups=TASKS.GROUPS,
                           reports={
                               "report.json": _read_json("models/report.json"),
                               "pan_report.json": _read_json("models/pan_report.json"),
                           })


@app.route("/calibrate")
def page_calibrate():
    return render_template("run.html", meta=_meta(), group="calibrate",
                           title="Calibration", groups=TASKS.GROUPS)


@app.route("/validation")
def page_validation():
    vs = VERDICTS.all_verdicts()
    return render_template("validation.html", meta=_meta(), group="validate",
                           verdicts=vs, headline=VERDICTS.headline(vs),
                           title="Validation", groups=TASKS.GROUPS,
                           ev=_read_csv("external_validation.csv"),
                           evp=_read_csv("external_validation_pan.csv"),
                           el=_read_csv("eluted_validation.csv"),
                           tr=_read_csv("trial_validation.csv"),
                           wl=_read_csv("wetlab_validation_test.csv"),
                           sc=_read_csv("scorecard.csv"))


def _read_json(rel):
    f = ROOT / rel
    if not f.exists():
        return None
    try:
        return json.loads(f.read_text())
    except (ValueError, OSError):
        return None


def _meta():
    _, pname = _params("wetlab")
    return {
        "platforms": [{"key": k, "label": p.label, "template": p.template,
                       "half_life_h": p.half_life_h,
                       "mhc1": p.mhc1_routing, "mhc2": p.mhc2_routing,
                       "innate": p.innate, "anti_vector": p.anti_vector,
                       "reacto": p.reactogenicity,
                       "dose": list(p.typical_dose_ug),
                       "self_adj": p.self_adjuvanting, "notes": p.notes}
                      for k, p in PF.PLATFORMS.items()],
        "adjuvants": [{"key": k, "label": a.label, "th1_bias": a.th1_bias,
                       "magnitude": a.magnitude, "duration": a.duration_days,
                       "ctl": a.ctl_boost, "ab": a.ab_boost,
                       "reacto": a.reactogenicity, "notes": a.notes}
                      for k, a in PF.ADJUVANTS.items()],
        "epitope_files": sorted(f.name for f in RESULTS.glob("*epitopes*.csv")),
        "series": [{"key": k, "label": l, "group": g} for k, l, g in SERIES],
        "calibration": pname,
    }


def _read_csv(name, limit=400):
    f = RESULTS / name
    if not f.exists():
        return None
    d = pd.read_csv(f).head(limit)
    return {"columns": list(d.columns),
            "rows": d.where(pd.notna(d), None).values.tolist()}


# --------------------------------------------------------------------------
# api
# --------------------------------------------------------------------------
@app.route("/api/meta")
def api_meta():
    return jsonify(_meta())


@app.route("/api/simulate", methods=["POST"])
def api_simulate():
    payload = request.get_json(force=True) or {}
    try:
        df, sim, info = run_sim(payload)
    except SystemExit as e:
        return jsonify({"error": str(e)}), 400
    except Exception as e:
        app.logger.error(traceback.format_exc())
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500

    keep = [k for k, _, _ in SERIES if k in df.columns]
    # thin the grid for the browser; the CSV download keeps every step
    step = max(1, len(df) // 900)
    d = df.iloc[::step]
    clones = {c: d[c].round(4).tolist() for c in df.columns
              if c.startswith(("TH_clone_", "TC_clone_"))}
    return jsonify({
        "day": d.day.round(3).tolist(),
        "series": {k: np.nan_to_num(d[k].values, posinf=0, neginf=0)
                    .round(6).tolist() for k in keep},
        "clones": clones,
        "summary": _summary(df, sim, info["injection_days"]),
        "info": info,
        "describe": PF.describe(sim.platform_key,
                                sim.adj.key, sim.dose_ug),
    })


@app.route("/api/compare", methods=["POST"])
def api_compare():
    """Same antigen and schedule, several formulations."""
    payload = request.get_json(force=True) or {}
    arms = payload.get("arms") or []
    if not 1 <= len(arms) <= 6:
        return jsonify({"error": "give between 1 and 6 arms"}), 400
    out = []
    for a in arms:
        cfg = dict(payload, **a)
        cfg.pop("arms", None)
        try:
            df, sim, info = run_sim(cfg)
        except Exception as e:
            out.append({"label": a.get("label") or "?", "error": str(e)})
            continue
        step = max(1, len(df) // 500)
        d = df.iloc[::step]
        out.append({
            "label": (a.get("label")
                      or f"{sim.pf.label} + {sim.adj.label}"),
            "platform": cfg.get("platform"), "adjuvant": cfg.get("adjuvant"),
            "day": d.day.round(3).tolist(),
            "IgG": d.IgG.round(4).tolist(),
            "TC_active": d.TC_active.round(4).tolist(),
            "TH_active": d.TH_active.round(4).tolist(),
            "summary": _summary(df, sim, info["injection_days"]),
        })
    return jsonify({"arms": out})


@app.route("/api/predict", methods=["POST"])
def api_predict():
    """Scan a protein sequence for epitopes with the trained models."""
    payload = request.get_json(force=True) or {}
    seq = "".join((payload.get("sequence") or "").split()).upper()
    seq = "".join(c for c in seq if c.isalpha())
    if len(seq) < 9:
        return jsonify({"error": "need at least 9 residues"}), 400
    if len(seq) > 4000:
        return jsonify({"error": "sequence longer than 4000 residues"}), 400
    try:
        import predict as PR
        models = PR.load_models(payload.get("alleles") or None)
        if not models:
            return jsonify({"error": "no trained models in models/; "
                                     "run model.py first"}), 503
        df = PR.scan(seq, models)
    except Exception as e:
        app.logger.error(traceback.format_exc())
        return jsonify({"error": f"{type(e).__name__}: {e}"}), 500

    df = df[df.call != "non-binder"].sort_values("pred_IC50_nM")
    name = (payload.get("name") or "query").strip()[:40] or "query"
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
    out = RESULTS / f"{safe}_epitopes.csv"
    df.to_csv(out, index=False)
    top = df.head(int(payload.get("limit", 60)))
    return jsonify({
        "n_binders": int(len(df)),
        "n_strong": int((df.call == "strong").sum()),
        "saved_as": out.name,
        "columns": ["allele", "mhc_class", "start", "peptide", "length",
                    "pred_IC50_nM", "percentile_rank", "call"],
        "rows": top[["allele", "mhc_class", "start", "peptide", "length",
                     "pred_IC50_nM", "percentile_rank", "call"]]
                 .round({"pred_IC50_nM": 2, "percentile_rank": 3})
                 .values.tolist(),
        "by_class": df.groupby("mhc_class").size().to_dict(),
    })


@app.route("/api/alleles")
def api_alleles():
    """Alleles the trained models cover, so the form offers only those.

    An allele the models were NOT trained on is deliberately not offered:
    zero-shot transfer to a new allele does not work here (PPV@1% 0.087 to
    0.121, Wilcoxon p = 0.45), so letting someone pick one would be offering
    a number that means nothing.
    """
    f = MODELS / "report.json"
    if not f.exists():
        return jsonify({"alleles": [], "trained": False})
    try:
        rep = json.loads(f.read_text())
    except (ValueError, OSError):
        return jsonify({"alleles": [], "trained": False})
    out = [{"allele": a, "mhc_class": v.get("mhc_class"),
            "auc": v.get("auc"), "n_train": v.get("n_train")}
           for a, v in sorted(rep.items())]
    return jsonify({"alleles": out, "trained": bool(out)})


@app.route("/api/design", methods=["POST"])
def api_design():
    """The whole pipeline in one call: a sequence in, a predicted response out.

    This exists because the two halves of this project were only ever joined
    by a CSV on disk: predict.py wrote epitopes, simulate.py read them, and
    the GUI made you do the same thing by hand. Here the antigen goes in once
    and the answer comes back once - the epitope prediction and the dynamic
    simulation are two stages of one request, not two tools.

    Accepts a pasted sequence, a UniProt accession, or the committed demo
    antigen, and every simulation setting. Returns the epitopes it found AND
    the simulated response, so the page can show which epitopes drove it.
    """
    payload = request.get_json(force=True) or {}
    source = (payload.get("antigen_source") or "sequence").strip()

    # --- 1. get the antigen ------------------------------------------------
    name = (payload.get("name") or "").strip()[:40]
    seq = ""
    if source == "uniprot":
        acc = "".join(c for c in str(payload.get("accession", ""))
                      if c.isalnum() or c in "-_")[:20]
        if not acc:
            abort(400, "give a UniProt accession")
        try:
            import predict as PR
            fetched, seq = PR.fetch_uniprot(acc)
        except Exception as e:
            abort(502, f"could not fetch {acc}: {e}")
        name = name or fetched[:40] or acc
    elif source == "demo":
        name = name or "demo (SARS-CoV-2 spike)"
    else:
        raw = payload.get("sequence") or ""
        if ">" in raw:                       # FASTA: drop the header line(s)
            raw = "\n".join(l for l in raw.splitlines()
                             if not l.startswith(">"))
            if not name:
                hdr = next((l for l in (payload.get("sequence") or "").splitlines()
                            if l.startswith(">")), "")
                name = hdr.lstrip(">").split("|")[-1].strip()[:40]
        seq = "".join(c for c in raw.upper() if c.isalpha())
        if len(seq) < 9:
            abort(400, "need at least 9 amino acid residues")
        if len(seq) > 4000:
            abort(400, f"sequence is {len(seq)} residues; the limit is 4000")
        bad = sorted(set(seq) - set("ACDEFGHIKLMNPQRSTVWY"))
        if bad:
            abort(400, "not amino acids: " + ", ".join(bad)
                       + ". Is this a DNA or RNA sequence?")
        # A, C, G, T and U are ALL valid amino acid letters (Ala, Cys, Gly,
        # Thr, Sec), so a pasted gene passes the check above and is silently
        # scanned as a protein. Catch it by composition instead: a real
        # protein of any length is essentially never built from four or five
        # residues alone.
        if len(seq) >= 24 and set(seq) <= set("ACGTU"):
            kind = "RNA" if "U" in seq and "T" not in seq else "DNA"
            abort(400,
                  f"This looks like a {kind} sequence, not a protein: every "
                  f"one of its {len(seq)} residues is A, C, G, T or U. Those "
                  f"are also valid amino acid letters, so nothing here can "
                  f"tell them apart for you. Translate it to protein first, "
                  f"then paste that.")
    name = name or "query"

    # --- 2. epitopes -------------------------------------------------------
    stages = {}
    if source == "demo":
        ep_src, ep_rows, ep_cols, counts = "demo", None, None, {}
        stages["predict"] = "used the committed demo epitope set"
    else:
        try:
            import predict as PR
            models = PR.load_models(payload.get("alleles") or None)
        except Exception as e:
            app.logger.error(traceback.format_exc())
            abort(500, f"{type(e).__name__}: {e}")
        if not models:
            abort(503, "No trained models yet. Press 'Set up everything' on "
                       "the Setup page - it needs no network, the data is "
                       "already in this repository.")
        try:
            pred = PR.scan(seq, models)
        except Exception as e:
            app.logger.error(traceback.format_exc())
            abort(500, f"{type(e).__name__}: {e}")
        binders = pred[pred.call != "non-binder"].sort_values("pred_IC50_nM")
        if binders.empty:
            abort(400, "No binders predicted for this sequence and these "
                       "alleles, so there is nothing to simulate. Try more "
                       "alleles, or a longer sequence.")
        safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
        (RESULTS / f"{safe}_epitopes.csv").write_text(binders.to_csv(index=False))
        ep_src = "inline"
        ep_rows = binders.to_dict("records")
        ep_cols = ["allele", "mhc_class", "start", "peptide", "length",
                   "pred_IC50_nM", "percentile_rank", "call"]
        counts = {"n_binders": int(len(binders)),
                  "n_strong": int((binders.call == "strong").sum()),
                  "by_class": binders.groupby("mhc_class").size().to_dict(),
                  "n_alleles": int(binders.allele.nunique()),
                  "saved_as": f"{safe}_epitopes.csv"}
        stages["predict"] = (f"{len(binders):,} binders over "
                             f"{binders.allele.nunique()} alleles")

    # --- 3. simulate -------------------------------------------------------
    cfg = dict(payload)
    cfg["epitope_source"] = ep_src
    if ep_rows is not None:
        cfg["epitope_rows"] = ep_rows
    try:
        df, sim, info = run_sim(cfg)
    except ValueError as e:
        abort(400, str(e))
    except Exception as e:
        app.logger.error(traceback.format_exc())
        abort(500, f"{type(e).__name__}: {e}")

    keep = [k for k, _, _ in SERIES if k in df.columns]
    step = max(1, len(df) // 900)
    d = df.iloc[::step]
    top = (pd.DataFrame(ep_rows)[ep_cols]
           .head(int(payload.get("limit", 80)))
           .round({"pred_IC50_nM": 2, "percentile_rank": 3})
           .values.tolist()) if ep_rows is not None else []
    return jsonify({
        "antigen": {"name": name, "length": len(seq) or None,
                    "source": source},
        "epitopes": {"columns": ep_cols, "rows": top, **counts},
        "day": d.day.round(3).tolist(),
        "series": {k: np.nan_to_num(d[k].values, posinf=0, neginf=0)
                    .round(6).tolist() for k in keep},
        "summary": _summary(df, sim, info["injection_days"]),
        "info": info,
        "stages": stages,
        "describe": PF.describe(sim.platform_key, sim.adj.key, sim.dose_ug),
    })


@app.route("/api/uniprot/<acc>")
def api_uniprot(acc):
    """Fetch a sequence by accession, so a user need not paste FASTA."""
    acc = "".join(c for c in acc if c.isalnum() or c in "-_")[:20]
    try:
        import predict as PR
        name, seq = PR.fetch_uniprot(acc)
    except Exception as e:
        return jsonify({"error": f"could not fetch {acc}: {e}"}), 502
    return jsonify({"accession": acc, "name": name,
                    "length": len(seq), "sequence": seq})


@app.route("/api/timecourse.csv", methods=["POST"])
def api_timecourse_csv():
    payload = request.get_json(force=True) or {}
    df, sim, _ = run_sim(payload)
    buf = io.StringIO()
    df.to_csv(buf, index=False)
    return Response(buf.getvalue(), mimetype="text/csv",
                    headers={"Content-Disposition":
                             "attachment; filename=timecourse.csv"})


@app.route("/api/results/<name>")
def api_results(name):
    f = RESULTS / Path(name).name
    if not f.exists() or f.suffix not in (".csv", ".png", ".json", ".txt"):
        abort(404)
    return send_file(f)


# --------------------------------------------------------------------------
# jobs: run the project's own scripts and stream their output
# --------------------------------------------------------------------------
def _file_info(rel):
    f = ROOT / rel
    if not f.exists():
        return {"path": rel, "exists": False}
    st = f.stat()
    info = {"path": rel, "exists": True, "bytes": st.st_size,
            "mtime": st.st_mtime, "stale": TASKS.stale_because(rel)}
    if f.suffix == ".csv":
        try:                            # row count without loading the file
            with f.open("rb") as fh:
                info["rows"] = max(sum(1 for _ in fh) - 1, 0)
        except OSError:
            pass
    return info


@app.route("/api/status")
def api_status():
    """What exists on disk, what each task needs, and what is running.

    The browser polls this so a page never shows a button as available when
    its inputs are missing, and never shows a result as current when
    something upstream of it has been rebuilt since.
    """
    specs = {}
    for name, spec in TASKS.TASKS.items():
        specs[name] = {
            "name": name, "group": spec["group"], "label": spec["label"],
            "about": spec["about"], "est": spec["est"],
            "network": bool(spec.get("network")),
            "requires": spec["requires"], "produces": spec["produces"],
            "fields": spec.get("fields", []),
            "missing": TASKS.missing(name),
            "outputs": [_file_info(p) for p in spec["produces"]],
        }
    return jsonify({
        "tasks": specs,
        "groups": [{"key": k, "label": l, "about": a} for k, l, a in TASKS.GROUPS],
        "plans": {n: {"name": n, "group": p["group"], "label": p["label"],
                      "about": p["about"], "est": p["est"],
                      "steps": [{"name": sn, "label": sl}
                                for sn, sl in p["steps"]],
                      "missing": TASKS.plan_missing(n),
                      "outputs": [_file_info(o) for o in p["produces"]]}
                  for n, p in TASKS.PLANS.items()},
        "data": DATAIO.inventory(),
        "inventory": [dict(_file_info(rel), label=label)
                      for rel, label in TASKS.INVENTORY],
        "jobs": [j.snapshot(log_lines=0) for j in MANAGER.recent(24)],
        "param_files": sorted(f"models/{f.name}"
                              for f in MODELS.glob("calibrated*.json")),
        "candidates": _candidates(),
        "shipped": _param_score("models/calibrated_wetlab.json"),
    })


def _param_score(rel):
    """The training score stored in a fitted parameter file, if any."""
    f = ROOT / rel
    if not f.exists():
        return None
    try:
        d = json.loads(f.read_text())
    except (ValueError, OSError):
        return None
    return {"path": rel, "score": d.get("score"),
            "mtime": f.stat().st_mtime}


def _candidates():
    return sorted(
        (c for c in (_param_score(f"models/{f.name}")
                     for f in MODELS.glob("calibrated_candidate_*.json")) if c),
        key=lambda c: -(c["mtime"] or 0))


@app.route("/api/promote", methods=["POST"])
def api_promote():
    """Make a calibration candidate the parameter set everything else uses.

    Deliberately explicit and deliberately separate from running the search.
    The previous set is kept as models/calibrated_wetlab_prev.json so this is
    reversible, and a WORSE training score is allowed through with a warning
    rather than blocked: a lower training score is not the same as a better
    model, and the held-out split on the Validation page is the thing that
    decides. Refusing on training score would be the GUI pretending to know
    better than the out-of-sample test.
    """
    body = request.get_json(force=True) or {}
    rel = str(body.get("file", ""))
    src = (ROOT / rel).resolve()
    allowed = (ROOT / "models").resolve()
    if (src.parent != allowed or not src.name.startswith("calibrated_candidate_")
            or src.suffix != ".json" or not src.exists()):
        abort(400, "not a calibration candidate in models/")
    cand = _param_score(rel)
    cur = _param_score("models/calibrated_wetlab.json")
    dest = ROOT / "models" / "calibrated_wetlab.json"
    if dest.exists():
        (ROOT / "models" / "calibrated_wetlab_prev.json").write_text(
            dest.read_text())
    dest.write_text(src.read_text())
    warn = None
    if cand and cur and cand["score"] is not None and cur["score"] is not None:
        if cand["score"] > cur["score"]:
            warn = (f"the promoted candidate's training score "
                    f"({cand['score']:.4g}) is WORSE than the set it replaced "
                    f"({cur['score']:.4g}); the previous set is kept as "
                    f"models/calibrated_wetlab_prev.json")
    return jsonify({"promoted": rel, "previous": cur, "warning": warn,
                    "note": "Re-run the held-out validation and the platform "
                            "benchmarks: every result computed from the old "
                            "parameters is now stale."})


@app.route("/api/jobs", methods=["POST"])
def api_job_start():
    body = request.get_json(force=True) or {}
    name = str(body.get("task", ""))
    if name not in TASKS.TASKS:
        abort(400, f"unknown task {name!r}")
    miss = TASKS.missing(name)
    if miss:
        abort(400, "missing input: " + ", ".join(miss)
                   + ". Download it on the Data page first.")
    try:
        label, argv, produces, parser = TASKS.build(name, body.get("args") or {})
    except ValueError as e:
        abort(400, str(e))
    try:
        job = MANAGER.submit(name, label, argv, ROOT, produces, parser)
    except RuntimeError as e:
        abort(409, str(e))
    return jsonify(job.snapshot(log_lines=0))


@app.route("/api/plans", methods=["POST"])
def api_plan_start():
    """Start a one-click plan: several steps, one job, one progress view."""
    body = request.get_json(force=True) or {}
    name = str(body.get("plan", ""))
    if name not in TASKS.PLANS:
        abort(400, f"unknown plan {name!r}")
    miss = TASKS.plan_missing(name)
    if miss:
        abort(400, "missing input: " + ", ".join(miss))
    label, steps, produces = TASKS.build_plan(name)
    try:
        job = MANAGER.submit_plan(name, label, steps, ROOT, produces)
    except RuntimeError as e:
        abort(409, str(e))
    return jsonify(job.snapshot(log_lines=0))


@app.route("/api/jobs/<job_id>")
def api_job(job_id):
    job = MANAGER.get(job_id)
    if job is None:
        abort(404, "no such job (the server may have restarted)")
    try:
        n = min(max(int(request.args.get("log", 200)), 0), 600)
    except ValueError:
        n = 200
    return jsonify(job.snapshot(log_lines=n))


@app.route("/api/jobs/<job_id>/cancel", methods=["POST"])
def api_job_cancel(job_id):
    job = MANAGER.get(job_id)
    if job is None:
        abort(404, "no such job")
    job.cancel()
    return jsonify(job.snapshot(log_lines=0))


@app.errorhandler(400)
@app.errorhandler(404)
@app.errorhandler(409)
@app.errorhandler(500)
@app.errorhandler(502)
@app.errorhandler(503)
def _err(e):
    return jsonify({"error": getattr(e, "description", str(e))}), e.code


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", default="127.0.0.1")
    ap.add_argument("--port", type=int, default=5000)
    ap.add_argument("--debug", action="store_true")
    a = ap.parse_args()
    print(f"immune_sim GUI -> http://{a.host}:{a.port}")
    app.run(host=a.host, port=a.port, debug=a.debug, threaded=True)
