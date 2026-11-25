#!/usr/bin/env python3
"""
T cell immunogenicity model: will a peptide actually elicit a T cell response?

This is a different question from binding. A peptide can bind MHC tightly and
still not be immunogenic (no T cell in the repertoire, tolerance, processing).

DATA HYGIENE: the IEDB tcell table contains 'qualitative binding' assays
alongside functional ones. Training on those would be circular - they measure
binding, not response. Only FUNCTIONAL assays are used here:
    IFNg release, activation, cytotoxicity, proliferation, IL-2 release

Features: BLOSUM62 peptide encoding + predicted binding affinity from the
stage-1 model (so immunogenicity is modelled on top of binding, not instead).

Validation: stratified CV with peptide-disjoint folds, ROC-AUC and PR-AUC
with bootstrap confidence intervals.
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from model import CANON_LEN, encode, inv_log50k

HERE = Path(__file__).parent
DATA = HERE / "data"
MODELS = HERE / "models"
MODELS.mkdir(exist_ok=True)

FUNCTIONAL = {"IFNg release", "activation", "cytotoxicity",
              "proliferation", "IL-2 release", "IFNg"}
MIN_ROWS = 300
N_FOLDS = 5
N_BOOT = 1000
SEED = 0


def load_binding_models():
    out = {}
    for p in MODELS.glob("HLA*.pkl"):
        if p.name == "immunogenicity.pkl":
            continue
        with open(p, "rb") as f:
            d = pickle.load(f)
        out[d["allele"]] = d
    return out


def load_tcell():
    d = pd.read_csv(DATA / "tcell_immunogenicity.csv")
    before = len(d)
    d = d[d.assay.isin(FUNCTIONAL)].copy()
    print(f"  kept {len(d):,}/{before:,} rows from functional assays only "
          f"(dropped {before - len(d):,} binding-type rows)")

    # Collapse repeat measurements of the same peptide/allele by majority call.
    g = (d.groupby(["allele", "mhc_class", "peptide", "length"], as_index=False)
           .agg(pos_frac=("positive", "mean"), n_assays=("positive", "size")))
    g = g[(g.pos_frac <= 0.25) | (g.pos_frac >= 0.75)]       # drop ambiguous
    g["label"] = (g.pos_frac >= 0.75).astype(int)
    print(f"  {len(g):,} unique peptide-allele pairs after consensus filter")
    return g


def add_binding_feature(df, bmodels):
    """Predicted log50k affinity as an extra feature."""
    df = df.copy()
    df["pred_affinity"] = np.nan
    for allele, sub in df.groupby("allele"):
        bm = bmodels.get(allele)
        if bm is None:
            continue
        L = bm["length"]
        ok = sub[sub.length == L]
        if ok.empty:
            continue
        df.loc[ok.index, "pred_affinity"] = bm["model"].predict(
            encode(ok.peptide.values, L))
    return df


def bootstrap_ci(y, p, metric, n=N_BOOT, seed=SEED):
    rng = np.random.default_rng(seed)
    vals = []
    idx = np.arange(len(y))
    for _ in range(n):
        s = rng.choice(idx, size=len(idx), replace=True)
        if len(np.unique(y[s])) < 2:
            continue
        vals.append(metric(y[s], p[s]))
    if not vals:
        return (float("nan"), float("nan"))
    return (float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5)))


def train_all():
    bmodels = load_binding_models()
    df = add_binding_feature(load_tcell(), bmodels)
    report, curves = {}, {}

    for (allele, mclass), sub in df.groupby(["allele", "mhc_class"]):
        L = CANON_LEN.get(mclass)
        sub = sub[(sub.length == L) & sub.pred_affinity.notna()]
        if len(sub) < MIN_ROWS or sub.label.nunique() < 2:
            print(f"  {allele:16s} skipped ({len(sub)} usable rows)")
            continue

        X = np.hstack([encode(sub.peptide.values, L),
                       sub.pred_affinity.values.reshape(-1, 1)])
        y = sub.label.values
        groups = sub.peptide.values        # peptide-disjoint folds

        oof = np.zeros(len(y))
        cv = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True,
                                  random_state=SEED)
        for tr, te in cv.split(X, y, groups):
            m = HistGradientBoostingClassifier(
                max_iter=300, learning_rate=0.08, max_leaf_nodes=31,
                l2_regularization=1.0, early_stopping=True,
                validation_fraction=0.1, random_state=SEED)
            m.fit(X[tr], y[tr])
            oof[te] = m.predict_proba(X[te])[:, 1]

        auc = roc_auc_score(y, oof)
        ap = average_precision_score(y, oof)
        auc_lo, auc_hi = bootstrap_ci(y, oof, roc_auc_score)
        ap_lo, ap_hi = bootstrap_ci(y, oof, average_precision_score)

        # affinity alone, as a baseline to beat
        base_auc = roc_auc_score(y, sub.pred_affinity.values)

        report[allele] = {
            "mhc_class": mclass, "n": int(len(sub)),
            "positive_rate": round(float(y.mean()), 3),
            "roc_auc": round(float(auc), 3),
            "roc_auc_ci": [round(auc_lo, 3), round(auc_hi, 3)],
            "pr_auc": round(float(ap), 3),
            "pr_auc_ci": [round(ap_lo, 3), round(ap_hi, 3)],
            "baseline_affinity_auc": round(float(base_auc), 3),
            "lift_over_affinity": round(float(auc - base_auc), 3),
        }
        curves[allele] = {"y": y.tolist(), "p": oof.tolist(),
                          "affinity": sub.pred_affinity.tolist()}

        # Final model on all data for downstream prediction.
        m = HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.08, max_leaf_nodes=31,
            l2_regularization=1.0, early_stopping=True,
            validation_fraction=0.1, random_state=SEED)
        m.fit(X, y)
        with open(MODELS / f"immuno_{allele.replace('*','_').replace(':','_')}.pkl",
                  "wb") as f:
            pickle.dump({"model": m, "length": L, "allele": allele,
                         "mhc_class": mclass}, f)

        print(f"  {allele:16s} n={len(sub):6d} pos={y.mean():.2f}  "
              f"ROC-AUC={auc:.3f} [{auc_lo:.3f}-{auc_hi:.3f}]  "
              f"PR-AUC={ap:.3f}  (affinity alone {base_auc:.3f})")

    with open(MODELS / "immuno_report.json", "w") as f:
        json.dump(report, f, indent=2)
    with open(MODELS / "immuno_curves.json", "w") as f:
        json.dump(curves, f)
    return report


if __name__ == "__main__":
    print("T cell immunogenicity - functional assays only")
    train_all()
    print(f"\n  report -> {MODELS/'immuno_report.json'}")
