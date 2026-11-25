#!/usr/bin/env python3
"""
External validation of the binding model: LEAVE-ONE-PATHOGEN-OUT.

A random held-out split overstates how well a predictor will do on a NEW
pathogen, because peptides from the same organism share composition and are
often near-duplicates of training peptides. The realistic question is:

    train with an organism entirely absent, then predict its epitopes.

That is the actual use case - a novel pathogen arrives and you must rank its
epitopes with a model that has never seen it. This file measures that, with
confidence intervals, significance tests against baselines, probability
calibration, and multiple-testing correction across pathogens.

Run:  python3 external_validation.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from model import BINDER_NM, CANON_LEN, SEED, encode, inv_log50k, load  # noqa: E402
from stats import (auc_ci_delong, bh_fdr, brier, calibration,           # noqa: E402
                   delong_test, fmt_p, permutation_auc)

RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

MIN_TEST = 150          # usable test peptides needed to report a pathogen
MIN_TRAIN = 800
N_PATHOGENS = 8


def organism(df):
    return (df.source_organism.fillna("unknown")
              .str.split("(").str[0].str.strip())


def length_baseline(train, test, L):
    """
    A deliberately weak baseline: predict every peptide's affinity as the
    training mean. Any real model must beat this, and comparing against it
    makes the AUC gain interpretable rather than absolute.
    """
    return np.full(len(test), train.y.mean())


def aa_composition_baseline(train, test, L):
    """
    A stronger baseline: amino-acid composition only, ignoring position.
    Beating this shows the model has learned POSITIONAL motif information,
    not just that hydrophobic peptides bind better.
    """
    AAS = "ACDEFGHIKLMNPQRSTVWY"
    def comp(peps):
        X = np.zeros((len(peps), 20), dtype=np.float32)
        for i, p in enumerate(peps):
            for a in p:
                j = AAS.find(a)
                if j >= 0:
                    X[i, j] += 1
        return X / max(L, 1)
    m = HistGradientBoostingRegressor(max_iter=200, learning_rate=0.08,
                                      random_state=SEED)
    m.fit(comp(train.peptide.values), train.y.values)
    return m.predict(comp(test.peptide.values))


def run():
    df = load()
    raw = pd.read_csv(HERE / "data" / "binding_affinity.csv")
    raw["org"] = organism(raw)
    pep_org = (raw.groupby("peptide").org
                  .agg(lambda s: s.value_counts().index[0]))
    df["org"] = df.peptide.map(pep_org)
    df = df[df.org.notna() & (df.org != "unknown")]

    counts = df.groupby("org").size().sort_values(ascending=False)
    pathogens = [o for o in counts.index if counts[o] >= MIN_TEST][:N_PATHOGENS]
    print(f"Leave-one-pathogen-out, {len(pathogens)} pathogens\n")

    rows, curves = [], {}
    for org in pathogens:
        for mclass in ("I", "II"):
            L = CANON_LEN[mclass]
            sub = df[(df.mhc_class == mclass) & (df.length == L)]
            test = sub[sub.org == org]
            train = sub[sub.org != org]
            # drop any peptide that also appears in training under another
            # organism label, or the split is not clean
            test = test[~test.peptide.isin(set(train.peptide))]
            if len(test) < MIN_TEST or len(train) < MIN_TRAIN:
                continue
            y = (inv_log50k(test.y.values) < BINDER_NM).astype(int)
            if len(np.unique(y)) < 2:
                continue

            m = HistGradientBoostingRegressor(
                max_iter=400, learning_rate=0.08, max_leaf_nodes=31,
                l2_regularization=1.0, early_stopping=True,
                validation_fraction=0.1, random_state=SEED)
            m.fit(encode(train.peptide.values, L), train.y.values)
            p = m.predict(encode(test.peptide.values, L))

            base = aa_composition_baseline(train, test, L)
            auc, (lo, hi) = auc_ci_delong(y, p)
            _, _, diff, p_delong = delong_test(y, p, base)
            _, p_perm = permutation_auc(y, p, n=2000)
            ap = average_precision_score(y, p)

            # Probability calibration: map the score to P(binder) by its rank
            # among training predictions, then check the reliability curve.
            tr_scores = m.predict(encode(train.peptide.values, L))
            tr_y = (inv_log50k(train.y.values) < BINDER_NM).astype(int)
            order = np.argsort(tr_scores)
            cal_x, cal_y = tr_scores[order], tr_y[order]
            k = max(50, len(cal_x) // 20)
            smooth = np.convolve(cal_y, np.ones(k) / k, mode="same")
            prob = np.interp(p, cal_x, smooth)
            b = brier(y, prob)
            c = calibration(y, prob)

            rows.append({
                "pathogen": org, "mhc_class": mclass,
                "n_test": int(len(test)), "n_train": int(len(train)),
                "prevalence": round(float(y.mean()), 3),
                "roc_auc": round(float(auc), 3),
                "auc_lo": round(lo, 3), "auc_hi": round(hi, 3),
                "pr_auc": round(float(ap), 3),
                "baseline_auc": round(float(roc_auc_score(y, base)), 3),
                "auc_gain": round(float(diff), 3),
                "p_vs_baseline": float(p_delong),
                "p_vs_chance": float(p_perm),
                "brier": round(b["brier"], 4),
                "reliability": round(b["reliability"], 4),
                "ece": round(c["ece"], 3),
            })
            curves[f"{org}|{mclass}"] = {"y": y.tolist(), "p": p.tolist(),
                                         "prob": prob.tolist()}
            print(f"  {org[:34]:34s} {mclass:2s} n={len(test):5d} "
                  f"AUC {auc:.3f} [{lo:.3f}-{hi:.3f}]  "
                  f"vs comp baseline {roc_auc_score(y, base):.3f} "
                  f"(p={fmt_p(p_delong)})  ECE {c['ece']:.3f}")

    if not rows:
        sys.exit("no pathogen had enough data for a clean split")

    tab = pd.DataFrame(rows)
    for col, newcol in (("p_vs_baseline", "q_vs_baseline"),
                        ("p_vs_chance", "q_vs_chance")):
        adj, sig = bh_fdr(tab[col].values)
        tab[newcol] = adj.round(4)
        tab[newcol.replace("q_", "sig_")] = sig

    tab.to_csv(RESULTS / "external_validation.csv", index=False)
    (RESULTS / "external_validation_curves.json").write_text(json.dumps(curves))

    print("\n" + "=" * 78)
    print("SUMMARY - generalisation to an unseen pathogen")
    print("=" * 78)
    for mclass in ("I", "II"):
        s = tab[tab.mhc_class == mclass]
        if s.empty:
            continue
        print(f"\n  class {mclass}:  {len(s)} pathogens")
        print(f"    ROC-AUC         median {s.roc_auc.median():.3f}  "
              f"range {s.roc_auc.min():.3f}-{s.roc_auc.max():.3f}")
        print(f"    vs composition  median gain {s.auc_gain.median():+.3f}  "
              f"significant after BH in "
              f"{int(s.sig_vs_baseline.sum())}/{len(s)}")
        print(f"    above chance    "
              f"{int(s.sig_vs_chance.sum())}/{len(s)} after BH")
        print(f"    calibration     median ECE {s.ece.median():.3f} "
              f"(0 = perfectly calibrated probabilities)")
    print(f"\n  wrote {RESULTS/'external_validation.csv'}")
    return tab


if __name__ == "__main__":
    run()
