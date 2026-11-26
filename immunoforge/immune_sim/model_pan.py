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

import dataio

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
                 max_iter=400, hidden=(256, 128), allele_enc=None):
        self.base_kind = base
        self.realign = realign
        self.seed = seed
        self.n_alleles = n_alleles
        # Hierarchical allele features (locus / group / exact allele). With a
        # flat one-hot an unseen allele had no column and fell back to index
        # 0, i.e. was scored as another allele entirely.
        self.allele_enc = allele_enc
        self.max_iter = max_iter
        self.hidden = hidden
        self.model = None
        self.scaler = None

    # -- allele one-hot block -------------------------------------------
    def _with_allele(self, X, allele_idx):
        if self.allele_enc is not None:
            names = self.allele_enc.alleles
            labels = [names[i] if 0 <= i < len(names) else "?"
                      for i in np.asarray(allele_idx, dtype=np.int64)]
            return np.hstack([X, self.allele_enc.transform(labels)])
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

    def _fit_base(self, X, y, sw=None):
        if self.base_kind == "mlp":
            self.scaler = StandardScaler().fit(X)
            X = self.scaler.transform(X)
        m = self._new_base()
        # MLPRegressor takes no sample_weight; the GBM does. Eluted ligands
        # are down-weighted relative to measured affinities because they are
        # a different quantity (what a cell presented, not what binds in a
        # tube), so they should inform the model without overruling it.
        if sw is not None and self.base_kind != "mlp":
            m.fit(X, y, sample_weight=sw)
        else:
            m.fit(X, y)
        self.model = m

    def _predict_base(self, X):
        if self.scaler is not None:
            X = self.scaler.transform(X)
        return self.model.predict(X)

    # -- public ----------------------------------------------------------
    def fit(self, peptides, y, allele_idx=None, verbose=False, sw=None):
        peptides = list(peptides)
        y = np.asarray(y, dtype=np.float64)
        sw = None if sw is None else np.asarray(sw, dtype=np.float64)
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
            self._fit_base(self._with_allele(X, allele_idx), y, sw)
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

    def predict(self, peptides, allele_idx=None, chunk=40000):
        """
        Scored in chunks. Expanding every register of a large set at once is
        what makes this explode: 850,000 eluted-ligand candidates at ~4
        registers each is 3.4 million rows of 268 float32 columns, about
        3.6 GB for the feature matrix alone. Chunking bounds it.
        """
        peptides = list(peptides)
        n = len(peptides)
        if allele_idx is None:
            allele_idx = np.zeros(n, dtype=np.int64)
        allele_idx = np.asarray(allele_idx, dtype=np.int64)

        best = np.empty(n, dtype=np.float64)
        offs = np.zeros(n, dtype=np.int64)
        for lo in range(0, n, chunk):
            hi = min(lo + chunk, n)
            sub = peptides[lo:hi]
            ridx, roffs, RX = E.expand_all_offsets(sub)
            sc = self._predict_base(
                self._with_allele(RX, allele_idx[lo:hi][ridx]))
            b, arg = E.pool_max(ridx, sc, hi - lo)
            best[lo:hi] = b
            offs[lo:hi] = roffs[arg]
        self.last_offsets_ = offs
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
        self.allele_enc = E.AlleleEncoder(self.alleles)
        self.members = [
            CoreScan("gbm", realign=realign, seed=seed,
                     allele_enc=self.allele_enc,
                     max_iter=120 if quick else 500),
            CoreScan("mlp", realign=max(2, realign - 1), seed=seed,
                     allele_enc=self.allele_enc,
                     max_iter=40 if quick else 180,
                     hidden=(128,) if quick else (256, 128)),
        ]
        self.names = ["corescan_gbm", "corescan_mlp"]
        self.weights = None
        self.blend_auc = {}
        self.seed = seed

    def _ai(self, alleles):
        """
        Index into the known-allele list, or -1 for an allele never trained
        on. -1 must NOT silently become 0: that is the bug this encoding
        replaces. `_with_allele` turns -1 into "?", which matches no exact
        allele and falls back to the group and locus features.
        """
        return np.array([self.aidx.get(a, -1) for a in alleles],
                        dtype=np.int64)

    def allele_coverage(self, allele):
        return self.allele_enc.coverage(allele)

    def _fit_members(self, peptides, y, ai, verbose, sw=None):
        for name, m in zip(self.names, self.members):
            t = time.time()
            m.fit(peptides, y, ai, verbose=verbose, sw=sw)
            if verbose:
                print(f"    {name:14s} fitted in {time.time() - t:5.1f}s",
                      flush=True)

    def fit(self, peptides, y, alleles, verbose=True, blend_frac=0.15, sw=None):
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
        sw = None if sw is None else np.asarray(sw, dtype=np.float64)
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
            self._fit_members(peptides[~bm], y[~bm], ai[~bm], verbose=False,
                              sw=None if sw is None else sw[~bm])
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
        self._fit_members(peptides, y, ai, verbose, sw=sw)
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


def load_pan_model(mhc_class, models_dir=None):
    """
    Load a trained pan ensemble.

    Handles the case where the pickle was written by `python3 model_pan.py`,
    in which case the classes were recorded as `__main__.PanEnsemble` rather
    than `model_pan.PanEnsemble` and a plain `pickle.load` from another script
    raises AttributeError. The `__main__` block below avoids creating new
    pickles like that; this keeps older ones readable.
    """
    import __main__
    d = Path(models_dir or MODELS) / f"pan_{mhc_class}.pkl"
    if not d.exists():
        return None
    for name, cls in (("PanEnsemble", PanEnsemble), ("CoreScan", CoreScan)):
        if not hasattr(__main__, name):
            setattr(__main__, name, cls)
    with open(d, "rb") as fh:
        return pickle.load(fh)["ensemble"]


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


# Eluted ligands as a SECOND training signal.
#
# Class II affinity data is thin: 21,394 measurements over 4 alleles, against
# 97,746 eluted-ligand rows over 11 alleles that this project previously used
# only for validation. That imbalance is the most likely reason class II
# generalises so poorly, and training on both assays is what NetMHCIIpan 4.x
# does.
#
# The two assays measure different things - what binds in a tube, and what a
# cell actually presented - so eluted rows are mapped onto the same [0,1]
# scale at values that say "a presented peptide is at least a good binder"
# and "a peptide never detected is probably not", and carry less weight than
# a real measurement. They inform the model; they do not overrule it.
EL_POS_Y, EL_NEG_Y = 0.85, 0.05
EL_WEIGHT = 0.35
EL_DECOYS_PER_POS = 5          # the harvest ships 19; more only adds bulk


def load_eluted(mhc_class, exclude_peptides=None, seed=0):
    """Eluted ligands for one class, as (peptide, allele, y, weight) rows."""
    try:
        el = dataio.read("eluted_ligands")
    except FileNotFoundError:
        return None
    el = el[(el.mhc_class == mhc_class)
            & (el.length >= MIN_LEN[mhc_class])
            & (el.length <= MAX_LEN[mhc_class])].copy()
    if exclude_peptides:
        el = el[~el.peptide.isin(set(exclude_peptides))]
    if el.empty:
        return None
    pos = el[el.label == 1]
    neg = el[el.label == 0]
    want = len(pos) * EL_DECOYS_PER_POS
    if len(neg) > want:
        neg = neg.sample(n=want, random_state=seed)
    out = pd.concat([pos, neg], ignore_index=True)
    out["y"] = np.where(out.label == 1, EL_POS_Y, EL_NEG_Y)
    out["w"] = EL_WEIGHT
    return out[["peptide", "allele", "y", "w"]]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="small models, for checking the plumbing")
    ap.add_argument("--realign", type=int, default=4)
    ap.add_argument("--classes", default="I,II")
    ap.add_argument("--eluted", default="none",
                    help="classes to ALSO train on eluted ligands "
                         "('II', 'I,II', or 'none'). Default none: whether "
                         "this helps is being measured by "
                         "external_validation_pan.py --eluted, and the "
                         "shipped model should not change on a hunch before "
                         "that answers")
    args = ap.parse_args()

    from sklearn.metrics import roc_auc_score
    from scipy.stats import spearmanr
    from stats import auc_ci_delong, delong_test, fmt_p

    suffix = "_quick" if args.quick else ""
    if args.quick:
        print("  --quick: small models, one realignment pass. Writing to "
              "pan_*_quick.pkl so the shipped model is left alone. These "
              "numbers are NOT the real ones.")
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

        # Eluted ligands, if asked for. Any peptide in the affinity TEST set
        # is removed from them first: the two sources overlap, and a test
        # peptide arriving through the other assay would be leakage.
        peps = tr.peptide.values
        ys = tr.y.values
        als = tr.allele.values
        sw = np.ones(len(tr))
        if mclass in [c.strip() for c in args.eluted.split(",")]:
            el = load_eluted(mclass, exclude_peptides=set(te.peptide), seed=SEED)
            if el is not None and len(el):
                peps = np.concatenate([peps, el.peptide.values])
                ys = np.concatenate([ys, el.y.values])
                als = np.concatenate([als, el.allele.values])
                sw = np.concatenate([sw, el.w.values])
                alleles = sorted(set(alleles) | set(el.allele.unique()))
                print(f"  + {len(el):,} eluted-ligand rows over "
                      f"{el.allele.nunique()} alleles "
                      f"(weight {EL_WEIGHT}), test peptides excluded")

        ens = PanEnsemble(alleles, realign=args.realign, quick=args.quick)
        ens.fit(peps, ys, als, sw=sw)

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
        with open(MODELS / f"pan_{mclass}{suffix}.pkl", "wb") as f:
            pickle.dump({"ensemble": ens, "mhc_class": mclass,
                         "alleles": alleles}, f)

    # A --quick run writes pan_*_quick.pkl and pan_report_quick.json, never
    # the shipped files. Twice now a smoke test has silently replaced the
    # real model with a debug one - small models, one realignment pass - and
    # the second time it also dropped class I from the report, because a run
    # with --classes II rewrites the whole report from what that run did.
    # A debug run must not be able to touch what the application loads.
    if suffix:
        (MODELS / f"pan_report{suffix}.json").write_text(
            json.dumps(report, indent=2))
    else:
        # merge, so --classes II keeps the class I entry a previous run wrote
        f = MODELS / "pan_report.json"
        merged = {}
        if f.exists():
            try:
                merged = json.loads(f.read_text())
            except ValueError:
                merged = {}
        merged.update(report)
        f.write_text(json.dumps(merged, indent=2))
    print(f"\n  wrote {MODELS}/pan_report{suffix}.json and "
          f"pan_*{suffix}.pkl")
    print("\n  NOTE: these are random-split numbers and are optimistic. "
          "Run external_validation_pan.py for leave-one-pathogen-out.")


if __name__ == "__main__":
    # Run through the imported module namespace, so pickled models record
    # their classes as `model_pan.PanEnsemble` and load from other scripts.
    import model_pan
    model_pan.main()
