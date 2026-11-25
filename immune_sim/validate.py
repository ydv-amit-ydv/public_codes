#!/usr/bin/env python3
"""
Full validation of the PREDICTION layer: ROC and precision-recall curves
with bootstrap confidence bands, for both models.

  binding        - peptide/MHC affinity, held-out peptide-disjoint split
  immunogenicity - T cell response, cross-validated out-of-fold predictions

IMPORTANT SCOPE NOTE
--------------------
ROC/AUC here certifies the *prediction* layer only: which peptides bind, and
which elicit a T cell response. It says nothing about the dynamic model's
kinetics, which are calibrated against qualitative immunological benchmarks
(see calibrate.py) and are NOT fitted to measured time-course data, because
no titer time courses were harvested. Do not quote these AUCs as validation
of the simulated curves.
"""
import json
import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (auc, average_precision_score,
                             precision_recall_curve, roc_auc_score, roc_curve)

from model import CANON_LEN, encode, inv_log50k, load, split_by_peptide

HERE = Path(__file__).parent
MODELS = HERE / "models"
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

BINDER_NM = 500.0
N_BOOT = 500
SEED = 0

plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "axes.spines.top": False,
                     "axes.spines.right": False})


def boot_band(y, p, n=N_BOOT, seed=SEED):
    """Bootstrap ROC band on a common FPR grid."""
    rng = np.random.default_rng(seed)
    grid = np.linspace(0, 1, 101)
    curves, aucs = [], []
    idx = np.arange(len(y))
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[s])) < 2:
            continue
        f, t, _ = roc_curve(y[s], p[s])
        curves.append(np.interp(grid, f, t))
        aucs.append(roc_auc_score(y[s], p[s]))
    if not curves:
        return grid, None, None, (np.nan, np.nan)
    C = np.vstack(curves)
    return (grid, np.percentile(C, 2.5, axis=0), np.percentile(C, 97.5, axis=0),
            (float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5))))


def binding_holdout():
    """Re-derive the held-out predictions for every binding model."""
    df = load()
    out = {}
    for (allele, mclass), sub in df.groupby(["allele", "mhc_class"]):
        f = MODELS / f"{allele.replace('*','_').replace(':','_')}.pkl"
        if not f.exists():
            continue
        with open(f, "rb") as fh:
            m = pickle.load(fh)
        L = m["length"]
        sub = sub[sub.length == L]
        if len(sub) < 50:
            continue
        _, te = split_by_peptide(sub)
        if te.empty:
            continue
        p = m["model"].predict(encode(te.peptide.values, L))
        y = (inv_log50k(te.y.values) < BINDER_NM).astype(int)
        if len(np.unique(y)) < 2:
            continue
        out[allele] = {"y": y, "p": p, "mhc_class": mclass}
    return out


def immuno_oof():
    f = MODELS / "immuno_curves.json"
    if not f.exists():
        return {}
    raw = json.load(open(f))
    rep = json.load(open(MODELS / "immuno_report.json"))
    return {a: {"y": np.array(v["y"]), "p": np.array(v["p"]),
                "affinity": np.array(v["affinity"]),
                "mhc_class": rep[a]["mhc_class"]}
            for a, v in raw.items() if a in rep}


def panel(data, title, fname, baseline_key=None):
    """ROC (with bootstrap band) + PR curve, coloured by MHC class."""
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    rows = []
    order = sorted(data, key=lambda a: (data[a]["mhc_class"], a))
    cmap = plt.get_cmap("tab20")

    for i, allele in enumerate(order):
        d = data[allele]
        y, p = d["y"], d["p"]
        a_roc = roc_auc_score(y, p)
        a_pr = average_precision_score(y, p)
        f, t, _ = roc_curve(y, p)
        ls = "-" if d["mhc_class"] == "I" else "--"
        ax[0].plot(f, t, ls, lw=1.3, color=cmap(i % 20),
                   label=f"{allele}  {a_roc:.3f}")
        pr, rc, _ = precision_recall_curve(y, p)
        ax[1].plot(rc, pr, ls, lw=1.3, color=cmap(i % 20),
                   label=f"{allele}  {a_pr:.3f}")

        _, lo, hi, ci = boot_band(y, p)
        row = {"allele": allele, "mhc_class": d["mhc_class"], "n": len(y),
               "positive_rate": round(float(y.mean()), 3),
               "roc_auc": round(float(a_roc), 3),
               "roc_auc_lo": round(ci[0], 3), "roc_auc_hi": round(ci[1], 3),
               "pr_auc": round(float(a_pr), 3),
               "prevalence_baseline_pr": round(float(y.mean()), 3)}
        if baseline_key and baseline_key in d:
            row["baseline_auc"] = round(float(roc_auc_score(y, d[baseline_key])), 3)
        rows.append(row)

    ax[0].plot([0, 1], [0, 1], color="k", lw=0.8, alpha=0.4)
    ax[0].set_xlabel("false positive rate"); ax[0].set_ylabel("true positive rate")
    ax[0].set_title("ROC"); ax[0].legend(fontsize=6.2, loc="lower right", ncol=1)
    ax[1].set_xlabel("recall"); ax[1].set_ylabel("precision")
    ax[1].set_title("Precision-recall (solid = class I, dashed = class II)")
    ax[1].legend(fontsize=6.2, loc="upper right", ncol=1)
    fig.suptitle(title, y=1.0)
    fig.tight_layout(); fig.savefig(RESULTS / fname, bbox_inches="tight")
    plt.close(fig)
    return pd.DataFrame(rows)


def main():
    print("Validating the prediction layer\n")

    b = binding_holdout()
    tb = panel(b, "Peptide-MHC binding (held-out, peptide-disjoint)",
               "validation_binding.png")
    print("--- binding ---")
    print(tb.to_string(index=False))

    im = immuno_oof()
    if im:
        ti = panel(im, "T cell immunogenicity (cross-validated out-of-fold)",
                   "validation_immunogenicity.png", baseline_key="affinity")
        print("\n--- immunogenicity ---")
        print(ti.to_string(index=False))
        ti.to_csv(RESULTS / "validation_immunogenicity.csv", index=False)
    tb.to_csv(RESULTS / "validation_binding.csv", index=False)

    print("\n--- summary ---")
    for nm, t in (("binding", tb), ("immunogenicity", ti if im else None)):
        if t is None:
            continue
        for cls in ("I", "II"):
            s = t[t.mhc_class == cls]
            if s.empty:
                continue
            print(f"  {nm:15s} class {cls:2s}  mean ROC-AUC {s.roc_auc.mean():.3f}"
                  f"   mean PR-AUC {s.pr_auc.mean():.3f}"
                  f"   (prevalence {s.positive_rate.mean():.3f})")

    print("\n  SCOPE: these curves validate epitope PREDICTION only.")
    print("  The dynamic model's kinetics are calibrated to qualitative")
    print("  benchmarks, not fitted to measured titer data - see calibrate.py.")
    print(f"\n  wrote {RESULTS}/validation_*.png and .csv")


if __name__ == "__main__":
    main()
