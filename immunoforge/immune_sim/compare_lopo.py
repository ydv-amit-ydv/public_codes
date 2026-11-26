#!/usr/bin/env python3
"""Paired comparison of two leave-one-pathogen-out runs.

The two arms differ in one thing only - whether eluted ligands were added to
the training set - and every fold holds out the same pathogen, so the folds
pair exactly. A paired test is the right one: the pathogens differ from each
other far more than the arms differ from each other, and an unpaired
comparison would drown the effect in that.

The bar to clear is not zero. Rebuilding this model from an unchanged tree
once moved class II by 0.029 (see VALIDATION.md), so that is the floor an
improvement has to beat before it means anything.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FLOOR = 0.029          # the reproducibility floor measured earlier


def main(a="lopoII_base.csv", b="lopoII_eluted.csv"):
    fa, fb = RESULTS / a, RESULTS / b
    for f in (fa, fb):
        if not f.exists():
            sys.exit(f"missing {f.name} - that arm has not finished")
    A = pd.read_csv(fa).set_index("pathogen")
    B = pd.read_csv(fb).set_index("pathogen")
    shared = [p for p in A.index if p in B.index]
    if not shared:
        sys.exit("the two runs share no pathogens")
    A, B = A.loc[shared], B.loc[shared]

    print(f"{'pathogen':34s} {'baseline':>9s} {'+eluted':>9s} {'diff':>8s}")
    for p in shared:
        d = B.auc_pan[p] - A.auc_pan[p]
        print(f"{p[:34]:34s} {A.auc_pan[p]:9.3f} {B.auc_pan[p]:9.3f} {d:+8.3f}")

    diff = (B.auc_pan - A.auc_pan).values
    med = float(np.median(diff))
    print(f"\n  folds compared        {len(shared)}")
    print(f"  median AUC, baseline  {A.auc_pan.median():.3f}")
    print(f"  median AUC, +eluted   {B.auc_pan.median():.3f}")
    print(f"  median paired change  {med:+.3f}")
    print(f"  better in             {int((diff > 0).sum())}/{len(diff)} folds")
    if len(diff) >= 6:
        try:
            st, p = wilcoxon(diff)
            print(f"  Wilcoxon signed-rank  p = {p:.3f}")
        except ValueError:
            pass
    for col, label in (("sig_vs_comp", "beats the composition baseline"),):
        if col in A.columns and col in B.columns:
            print(f"  {label:34s} {int(A[col].sum())}/{len(A)} -> "
                  f"{int(B[col].sum())}/{len(B)}")

    print()
    if med > FLOOR:
        print(f"  VERDICT: +{med:.3f} clears the {FLOOR} reproducibility "
              f"floor. Worth shipping; change the --eluted default.")
    elif med < -FLOOR:
        print(f"  VERDICT: {med:+.3f} is a real loss. Do not ship.")
    else:
        print(f"  VERDICT: {med:+.3f} is inside the {FLOOR} floor that "
              f"rebuilding this model from an unchanged tree already moves "
              f"it by.\n           Not distinguishable from noise, so it "
              f"does not ship: more training\n           data that does not "
              f"measurably help is still a change to defend.")


if __name__ == "__main__":
    main(*sys.argv[1:3] if len(sys.argv) > 2 else ())
