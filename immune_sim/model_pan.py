#!/usr/bin/env python3
"""
Pan-allele, core-scanning peptide-MHC binding models.

Three changes from `model.py`, each with a reason to expect it to help on an
UNSEEN pathogen rather than only on a random split:

1. **Core scan with iterative realignment (NN-align).** The binding register
   is a latent variable. Start every peptide at its central 9-mer, fit, then
   re-assign each peptide to the register the model scores highest, and refit.
   A few rounds of this sharpen the motif instead of smearing it across the
   peptide. Class II gains most, because its groove is open-ended.

2. **Pan-allele pooling.** One model per class over all alleles, with the
   allele as an input. A per-allele model sees only its own rows; a pooled
   model learns that a hydrophobic P9 anchor is a hydrophobic P9 anchor
   whichever allele it is, then specialises. This is why NetMHCpan beats
   allele-specific predictors, and it matters most for low-data alleles.

3. **All peptide lengths.** `model.py` keeps only 9-mers / 15-mers and drops a
   quarter of the data. A core at an offset handles any length.

The ensemble averages gradient-boosted trees and a neural network over the
same features. They fail differently - trees on smooth extrapolation, the net
on sharp motif thresholds - so the mean is better calibrated than either.

    python3 model_pan.py                 # train and report
    python3 model_pan.py --quick         # smaller, for a smoke test
"""
import argparse
import json
import pickle
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler

import encoding as E
from model import BINDER_NM, SEED, inv_log50k, load, split_by_peptide

warnings.filterwarnings("ignore")
HERE = Path(__file__).parent
MODELS = HERE / "models"
MODELS.mkdir(exist_ok=True)

MIN_LEN = {"I": 8, "II": 11}      # below this the register is not meaningful
MAX_LEN = {"I": 14, "II": 25}


# --------------------------------------------------------------------------
# core-scanning estimator
# --------------------------------------------------------------------------
class CoreScan:
    """
    Regressor over 9-mer binding cores with the register treated as latent.

    fit()      alternates (fit on current registers) / (re-assign registers)
    predict()  scores every register and returns the best, which is the
               biological claim: a peptide binds if ANY register binds.
    """

    def __init__(self, base="gbm", realign=4, seed=SEED, n_alleles=0,
                 max_iter=400, hidden=(256, 128)):
        self.base_kind = base
        self.realign = realign
        self.seed = seed
        self.n_alleles = n_alleles
        self.max_iter = max_iter
        self.hidden = hidden
        self.model = None
        self.scaler = None

    # -- allele one-hot block -------------------------------------------
    def _with_allele(self, X, allele_idx):
        if not self.n_alleles:
            return X
        A = np.zeros((X.shape[0], self.n_alleles), dtype=np.float32)
        A[np.arange(X.shape[0]), np.asarray(allele_idx, dtype=np.int64)] = 1.0
        return np.hstack([X, A])

    def _new_base(self):
        if self.base_kind == "gbm":
            return HistGradientBoostingRegressor(
                max_iter=self.max_iter, learning_rate=0.06, max_leaf_nodes=63,
                min_samples_leaf=20, l2_regularization=1.0,
                early_stopping=True, validation_fraction=0.1,
                n_iter_no_change=25, random_state=self.seed)
        return MLPRegressor(
            hidden_layer_sizes=self.hidden, activation="relu",
            alpha=1e-4, learning_rate_init=1e-3, batch_size=256,
            max_iter=self.max_iter, early_stopping=True,
            validation_fraction=0.1, n_iter_no_change=12,
            random_state=self.seed)

    def _fit_base(self, X, y):
        if self.base_kind == "mlp":
            self.scaler = StandardScaler().fit(X)
            X = self.scaler.transform(X)
        m = self._new_base()
        m.fit(X, y)
        self.model = m

    def _predict_base(self, X):
        if self.scaler is not None:
            X = self.scaler.transform(X)
        return self.model.predict(X)

    # -- public ----------------------------------------------------------
    def fit(self, peptides, y, allele_idx=None, verbose=False):
        peptides = list(peptides)
        y = np.asarray(y, dtype=np.float64)
        if allele_idx is None:
            allele_idx = np.zeros(len(peptides), dtype=np.int64)
        allele_idx = np.asarray(allele_idx, dtype=np.int64)

        # register 0: the central 9-mer, the best guess before any model
        offs = np.array([max(0, (len(p) - E.CORE) // 2) for p in peptides])

        # all registers, encoded once and reused every round
        ridx, roffs, RX = E.expand_all_offsets(peptides)
        RXa = self._with_allele(RX, allele_idx[ridx])

        for r in range(max(1, self.realign)):
            X = E.encode_cores(peptides, offs)
            self._fit_base(self._with_allele(X, allele_idx), y)
            if r == self.realign - 1:
                break
            s = self._predict_base(RXa)
            _, arg = E.pool_max(ridx, s, len(peptides))
            new = roffs[arg]
            moved = float(np.mean(new != offs))
            offs = new
            if verbose:
                print(f"      realign {r + 1}: {moved:5.1%} of peptides "
                      f"changed register", flush=True)
            if moved < 0.01:
                break
        self.offsets_ = offs
        return self

    def predict(self, peptides, allele_idx=None):
        peptides = list(peptides)
        if allele_idx is None:
            allele_idx = np.zeros(len(peptides), dtype=np.int64)
        allele_idx = np.asarray(allele_idx, dtype=np.int64)
        ridx, roffs, RX = E.expand_all_offsets(peptides)
        s = self._predict_base(self._with_allele(RX, allele_idx[ridx]))
        best, arg = E.pool_max(ridx, s, len(peptides))
        self.last_offsets_ = roffs[arg]
        return best


class PanEnsemble:
    """
    Blend of a core-scanning GBM and a core-scanning network.

    Weights are FITTED, not equal. An unweighted mean is only right when the
    members are comparably good: in the first run here the network scored
    0.866 against the tree's 0.909 and dragged the mean down to 0.904, so the
    "ensemble" was worse than its best member. Non-negative least squares on
    an inner, peptide-disjoint blend split chooses the weights, which lets a
    weak member be given zero rather than a third of the vote.
    """

    def __init__(self, alleles, realign=4, seed=SEED, quick=False):
        self.alleles = list(alleles)
        self.aidx = {a: i for i, a in enumerate(self.alleles)}
        self.members = [
            CoreScan("gbm", realign=realign, seed=seed,
                     n_alleles=len(self.alleles),
                     max_iter=120 if quick else 500),
            CoreScan("mlp", realign=max(2, realign - 1), seed=seed,
                     n_alleles=len(self.alleles),
                     max_iter=40 if quick else 180,
                     hidden=(128,) if quick else (256, 128)),
        ]
        self.names = ["corescan_gbm", "corescan_mlp"]
        self.weights = None
        self.blend_auc = {}
        self.seed = seed

    def _ai(self, alleles):
        return np.array([self.aidx.get(a, 0) for a in alleles], dtype=np.int64)

    def _fit_members(self, peptides, y, ai, verbose):
        for name, m in zip(self.names, self.members):
            t = time.time()
            m.fit(peptides, y, ai, verbose=verbose)
            if verbose:
                print(f"    {name:14s} fitted in {time.time() - t:5.1f}s",
                      flush=True)

    def fit(self, peptides, y, alleles, verbose=True, blend_frac=0.15):
        """
        Two passes. First on a sub-split to learn the blend weights, then on
        everything with those weights kept. The blend split is
        peptide-disjoint, so a weight is not chosen on a peptide the members
        memorised.
        """
        from scipy.optimize import nnls
        from sklearn.metrics import roc_auc_score

        peptides = np.asarray(list(peptides), dtype=object)
        y = np.asarray(y, dtype=np.float64)
        ai = self._ai(alleles)

        peps = np.unique(peptides)
        rng = np.random.default_rng(self.seed)
        rng.shuffle(peps)
        held = set(peps[:max(1, int(len(peps) * blend_frac))])
        bm = np.array([p in held for p in peptides])

        if bm.sum() > 50 and (~bm).sum() > 200:
            if verbose:
                print(f"    blend split: {int((~bm).sum())} fit / "
                      f"{int(bm.sum())} weight-fitting rows", flush=True)
            self._fit_members(peptides[~bm], y[~bm], ai[~bm], verbose=False)
            Pm = np.vstack([m.predict(peptides[bm], ai[bm])
                            for m in self.members]).T
            w, _ = nnls(Pm, y[bm])
            if w.sum() <= 0:
                w = np.ones(len(self.members))
            self.weights = w / w.sum()
            yb = (inv_log50k(y[bm]) < BINDER_NM).astype(int)
            if len(np.unique(yb)) > 1:
                for i, n in enumerate(self.names):
                    self.blend_auc[n] = round(
                        float(roc_auc_score(yb, Pm[:, i])), 3)
            if verbose:
                print("    blend weights  " + "  ".join(
                    f"{n}={w_:.2f}" for n, w_ in zip(self.names,
                                                     self.weights)),
                      flush=True)
        else:
            self.weights = np.ones(len(self.members)) / len(self.members)

        # refit on everything, weights fixed
        self._fit_members(peptides, y, ai, verbose)
        return self

    def predict_members(self, peptides, alleles):
        ai = self._ai(alleles)
        return {n: m.predict(peptides, ai)
                for n, m in zip(self.names, self.members)}

    def predict(self, peptides, alleles):
        P = self.predict_members(peptides, alleles)
        M = np.vstack([P[n] for n in self.names])
        w = (self.weights if self.weights is not None
             else np.ones(len(self.names)) / len(self.names))
        return np.asarray(w, dtype=np.float64) @ M


# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
def load_wide(mhc_class):
    """All usable lengths for one MHC class, not just the canonical one."""
    df = load()
    df = df[(df.mhc_class == mhc_class)
            & (df.length >= MIN_LEN[mhc_class])
            & (df.length <= MAX_LEN[mhc_class])]
    return df.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="small models, for checking the plumbing")
    ap.add_argument("--realign", type=int, default=4)
    ap.add_argument("--classes", default="I,II")
    args = ap.parse_args()

    from sklearn.metrics import roc_auc_score
    from scipy.stats import spearmanr
    from stats import auc_ci_delong, delong_test, fmt_p

    report = {}
    for mclass in args.classes.split(","):
        df = load_wide(mclass)
        if df.empty:
            print(f"class {mclass}: no data")
            continue
        alleles = sorted(df.allele.unique())
        print(f"\nclass {mclass}: {len(df)} rows, {len(alleles)} alleles, "
              f"lengths {df.length.min()}-{df.length.max()}")
        print(f"  ({(df.length == {'I': 9, 'II': 15}[mclass]).sum()} at the "
              f"canonical length - model.py would have used only those)")

        tr, te = split_by_peptide(df)
        ens = PanEnsemble(alleles, realign=args.realign, quick=args.quick)
        ens.fit(tr.peptide.values, tr.y.values, tr.allele.values)

        y = (inv_log50k(te.y.values) < BINDER_NM).astype(int)
        members = ens.predict_members(te.peptide.values, te.allele.values)
        pred = ens.predict(te.peptide.values, te.allele.values)

        rows = {}
        for n, p in list(members.items()) + [("ensemble", pred)]:
            auc, (lo, hi) = auc_ci_delong(y, p)
            rows[n] = {"auc": round(float(auc), 3),
                       "lo": round(lo, 3), "hi": round(hi, 3),
                       "spearman": round(float(spearmanr(te.y.values, p)
                                                .statistic), 3)}
            print(f"    {n:14s} AUC {auc:.3f} [{lo:.3f}-{hi:.3f}]  "
                  f"rho {rows[n]['spearman']:.3f}")

        # per-allele breakdown, so a good mean cannot hide a dead allele
        per = {}
        for a in alleles:
            m = te.allele.values == a
            if m.sum() < 60 or len(np.unique(y[m])) < 2:
                continue
            per[a] = round(float(roc_auc_score(y[m], pred[m])), 3)
            print(f"      {a:16s} n={int(m.sum()):5d}  AUC {per[a]:.3f}")

        report[mclass] = {"n_rows": int(len(df)), "n_alleles": len(alleles),
                          "metrics": rows, "per_allele": per,
                          "blend_weights": dict(zip(ens.names,
                                                    [round(float(x), 3)
                                                     for x in ens.weights])),
                          "blend_split_auc": ens.blend_auc}
        with open(MODELS / f"pan_{mclass}.pkl", "wb") as f:
            pickle.dump({"ensemble": ens, "mhc_class": mclass,
                         "alleles": alleles}, f)

    (MODELS / "pan_report.json").write_text(json.dumps(report, indent=2))
    print(f"\n  wrote {MODELS/'pan_report.json'} and pan_*.pkl")
    print("\n  NOTE: these are random-split numbers and are optimistic. "
          "Run external_validation_pan.py for leave-one-pathogen-out.")


if __name__ == "__main__":
    main()
