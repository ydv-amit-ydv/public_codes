#!/usr/bin/env python3
"""
Validate the binding models against MASS-SPECTROMETRY ELUTED LIGANDS.

This is the most independent test in the repository. The models are trained on
in-vitro binding affinity - purified MHC, competitive assay, IC50 in nM. The
test set is peptides recovered from the surface of live cells by
immunoprecipitation and mass spectrometry, from 131 separate publications.
Different technology, different laboratories, and different biology: natural
presentation also needs proteasomal cleavage, TAP transport and loading, none
of which an affinity assay measures.

Reported per allele, for both `model.py` and `model_pan.py`:

  ROC-AUC with a DeLong confidence interval
  PR-AUC (average precision), which is the honest summary at 5% prevalence
  PPV in the top 1% and top 0.1% of scores - the metric the eluted-ligand
      benchmarks use, because it is what a vaccine designer actually does:
      take the top few peptides and test them
  a paired DeLong test between the two models on identical cases
  Benjamini-Hochberg correction across alleles

    python3 harvest_eluted.py && python3 validate_eluted.py
"""
import argparse
import json
import pickle
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
warnings.filterwarnings("ignore")

from model import CANON_LEN, encode                     # noqa: E402
from stats import (auc_ci_delong, bh_fdr, brier,         # noqa: E402
                   calibration, delong_test, fmt_p)

DATA = HERE / "data"
MODELS = HERE / "models"
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

MIN_POS = 40


def ppv_at(y, p, frac):
    """Precision among the top `frac` of scores: the benchmark's own metric."""
    k = max(1, int(round(len(p) * frac)))
    idx = np.argsort(p)[::-1][:k]
    return float(np.mean(y[idx])), k


def load_old_models():
    out = {}
    for f in MODELS.glob("*.pkl"):
        if f.name.startswith(("immuno_", "pan_")):
            continue
        try:
            with open(f, "rb") as fh:
                d = pickle.load(fh)
            out[d["allele"]] = d
        except Exception:
            continue
    return out


def load_pan():
    from model_pan import load_pan_model
    out = {}
    for mclass in ("I", "II"):
        m = load_pan_model(mclass)
        if m is not None:
            out[mclass] = m
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-per-allele", type=int, default=60000,
                    help="cap rows per allele so one allele cannot dominate")
    args = ap.parse_args()

    f = DATA / "eluted_ligands.csv"
    if not f.exists():
        sys.exit(f"missing {f} - run harvest_eluted.py first")
    df = pd.read_csv(f)
    meta = {}
    mf = DATA / "eluted_ligands_meta.json"
    if mf.exists():
        meta = json.loads(mf.read_text())

    old = load_old_models()
    pan = load_pan()
    if not old and not pan:
        sys.exit("no trained models found - run model.py / model_pan.py")

    print("Validation against mass-spectrometry eluted ligands")
    print(f"  {int((df.label == 1).sum()):,} ligands, "
          f"{int((df.label == 0).sum()):,} decoys, "
          f"{df.allele.nunique()} alleles, "
          f"{meta.get('publications', '?')} publications")
    print(f"  prevalence {df.label.mean():.3f}  "
          f"(decoys are unobserved, not measured, negatives)\n")

    rows = []
    for (allele, mclass), grp in df.groupby(["allele", "mhc_class"]):
        if grp.label.sum() < MIN_POS:
            continue
        if len(grp) > args.max_per_allele:
            pos = grp[grp.label == 1]
            neg = grp[grp.label == 0].sample(
                min((grp.label == 0).sum(), args.max_per_allele - len(pos)),
                random_state=0)
            grp = pd.concat([pos, neg])
        y = grp.label.values.astype(int)
        peps = grp.peptide.values
        rec = {"allele": allele, "mhc_class": mclass,
               "n": int(len(grp)), "n_ligands": int(y.sum()),
               "prevalence": round(float(y.mean()), 4)}
        # How much of this allele the model has actually seen. Reporting one
        # median over all alleles hides the only thing that matters here:
        # whether the allele was in the training set at all.
        cov = "unknown"
        if mclass in pan and getattr(pan[mclass], "allele_enc", None):
            c = pan[mclass].allele_coverage(allele)
            cov = ("trained" if c["allele"] else
                   "group_only" if c["group"] else
                   "locus_only" if c["locus"] else "unknown")
        rec["coverage"] = cov

        # --- pan model: every length, pooled alleles ---------------------
        p_pan = None
        if mclass in pan:
            try:
                p_pan = pan[mclass].predict(peps, [allele] * len(peps))
                auc, (lo, hi) = auc_ci_delong(y, p_pan)
                rec.update({
                    "pan_auc": round(float(auc), 3),
                    "pan_lo": round(lo, 3), "pan_hi": round(hi, 3),
                    "pan_pr_auc": round(
                        float(average_precision_score(y, p_pan)), 3),
                    "pan_ppv_1pct": round(ppv_at(y, p_pan, 0.01)[0], 3),
                    "pan_ppv_01pct": round(ppv_at(y, p_pan, 0.001)[0], 3)})
            except Exception as e:
                rec["pan_error"] = f"{type(e).__name__}: {e}"

        # --- model.py: canonical length only -----------------------------
        p_old = None
        L = CANON_LEN[mclass]
        m = grp.length.values == L
        if allele in old and m.sum() > 50 and 0 < y[m].sum() < m.sum():
            try:
                pred = old[allele]["model"].predict(encode(peps[m], L))
                p_old = pred
                auc, (lo, hi) = auc_ci_delong(y[m], pred)
                rec.update({
                    "old_auc": round(float(auc), 3),
                    "old_lo": round(lo, 3), "old_hi": round(hi, 3),
                    "old_pr_auc": round(
                        float(average_precision_score(y[m], pred)), 3),
                    "old_ppv_1pct": round(ppv_at(y[m], pred, 0.01)[0], 3),
                    "n_canonical": int(m.sum())})
            except Exception as e:
                rec["old_error"] = f"{type(e).__name__}: {e}"

        # --- paired test on identical cases ------------------------------
        if p_pan is not None and p_old is not None:
            try:
                ys = y[m]
                if 0 < ys.sum() < len(ys):
                    _, _, diff, pv = delong_test(ys, p_pan[m], p_old)
                    rec.update({"auc_diff": round(float(diff), 3),
                                "p_pan_vs_old": float(pv),
                                "pan_auc_canonical": round(
                                    float(roc_auc_score(ys, p_pan[m])), 3)})
            except Exception:
                pass

        rows.append(rec)
        print(f"  {allele:16s} {mclass:2s} [{rec['coverage']:10s}] "
              f"n={rec['n']:7,d} "
              f"lig={rec['n_ligands']:6,d}  "
              f"pan AUC {rec.get('pan_auc', float('nan')):.3f}  "
              f"PR {rec.get('pan_pr_auc', float('nan')):.3f}  "
              f"PPV@1% {rec.get('pan_ppv_1pct', float('nan')):.3f}"
              + (f"   old {rec['old_auc']:.3f}"
                 f"  diff {rec.get('auc_diff', float('nan')):+.3f}"
                 f" (p={fmt_p(rec.get('p_pan_vs_old', float('nan')))})"
                 if "old_auc" in rec else "   old n/a"), flush=True)

    if not rows:
        sys.exit("no allele had enough ligands")

    tab = pd.DataFrame(rows)
    if "p_pan_vs_old" in tab:
        ok = tab.p_pan_vs_old.notna()
        if ok.sum():
            adj, sig = bh_fdr(tab.loc[ok, "p_pan_vs_old"].values)
            tab.loc[ok, "q_pan_vs_old"] = adj.round(4)
            tab.loc[ok, "sig_pan_vs_old"] = sig
    tab.to_csv(RESULTS / "eluted_validation.csv", index=False)

    print("\n" + "=" * 78)
    print("SUMMARY - predicting what a live cell actually presents")
    print("=" * 78)
    for mclass in ("I", "II"):
        s = tab[tab.mhc_class == mclass]
        if s.empty:
            continue
        print(f"\n  class {mclass}: {len(s)} alleles, "
              f"{int(s.n_ligands.sum()):,} ligands")
        if "pan_auc" in s and s.pan_auc.notna().any():
            print(f"    pan  ROC-AUC  median {s.pan_auc.median():.3f}  "
                  f"range {s.pan_auc.min():.3f}-{s.pan_auc.max():.3f}")
            print(f"    pan  PR-AUC   median {s.pan_pr_auc.median():.3f}  "
                  f"(prevalence {s.prevalence.median():.3f})")
            print(f"    pan  PPV@1%   median {s.pan_ppv_1pct.median():.3f}  "
                  f"PPV@0.1% median {s.pan_ppv_01pct.median():.3f}")
        print("\n    stratified by how much of the allele the model saw:")
        for cov in ("trained", "group_only", "locus_only", "unknown"):
            g = s[s.coverage == cov]
            if g.empty or g.pan_auc.notna().sum() == 0:
                continue
            print(f"      {cov:11s} {len(g):2d} alleles  "
                  f"ROC-AUC median {g.pan_auc.median():.3f}  "
                  f"PPV@1% median {g.pan_ppv_1pct.median():.3f}")
        c = s[s.get("auc_diff", pd.Series(dtype=float)).notna()] \
            if "auc_diff" in s else pd.DataFrame()
        if not c.empty:
            nsig = int(c.get("sig_pan_vs_old", pd.Series(dtype=bool)).sum())
            print(f"    old  ROC-AUC  median {c.old_auc.median():.3f}  "
                  f"(canonical length only)")
            print(f"    pan - old     median {c.auc_diff.median():+.3f}  "
                  f"pan better in {int((c.auc_diff > 0).sum())}/{len(c)}  "
                  f"significant after BH in {nsig}/{len(c)}")
    print(f"\n  wrote {RESULTS/'eluted_validation.csv'}")
    print("\n  Caveat, stated because it bounds the numbers above: decoys are")
    print("  peptides that mass spectrometry did not detect, not peptides")
    print("  measured as unpresented. Absolute ROC-AUC therefore depends on")
    print("  the decoy ratio; PPV in the top 1% is far less sensitive to it")
    print("  and is the figure to quote.")
    return tab


if __name__ == "__main__":
    main()
