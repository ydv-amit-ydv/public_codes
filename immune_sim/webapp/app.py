#!/usr/bin/env python3
"""
Flask front end for immune_sim.

    cd immune_sim/webapp && python3 app.py
    open http://127.0.0.1:5000

Four things it does:

  /            platform + adjuvant + schedule -> simulated T/B/antibody course
  /epitopes    paste or fetch a protein -> ranked epitopes -> feed the sim
  /compare     several formulations on one axis
  /validation  the evidence behind every claim, read from results/

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
    inj = _schedule(payload)
    days = float(payload.get("days", 240))
    n_steps = int(max(30, min(days * 3, 3000)))      # 8 h per step

    sim = Sim(ep, p=params, injections=inj,
              platform=payload.get("platform", "protein"),
              adjuvant=payload.get("adjuvant") or None,
              dose_ug=(float(payload["dose_ug"])
                       if payload.get("dose_ug") else None),
              heterologous=payload.get("heterologous") or None)
    df = sim.run(n_steps=n_steps)
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
def page_simulate():
    return render_template("simulate.html", meta=_meta())


@app.route("/epitopes")
def page_epitopes():
    return render_template("epitopes.html", meta=_meta())


@app.route("/compare")
def page_compare():
    return render_template("compare.html", meta=_meta())


@app.route("/validation")
def page_validation():
    return render_template("validation.html", meta=_meta(),
                           ev=_read_csv("external_validation.csv"),
                           wl=_read_csv("wetlab_validation_test.csv"),
                           sc=_read_csv("scorecard.csv"))


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


@app.errorhandler(400)
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
