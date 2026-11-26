#!/usr/bin/env python3
"""
Peptide-MHC binding model trained on the harvested IEDB data.

Target:  log50k = 1 - log(IC50_nM) / log(50000)   [the standard IEDB transform]
         1.0 = very strong binder, 0.0 = non-binder at 50 uM.
Encoding: BLOSUM62 row per residue, fixed-length peptides, per allele.
Model:   gradient-boosted trees, one per allele.

Evaluation is on a peptide-disjoint held-out split, so no peptide seen in
training appears in test (prevents leakage from duplicate measurements).
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score

import dataio

HERE = Path(__file__).parent
DATA = HERE / "data"
MODELS = HERE / "models"
MODELS.mkdir(exist_ok=True)

CANON_LEN = {"I": 9, "II": 15}   # dominant length per MHC class
BINDER_NM = 500.0                # standard binder/non-binder cut-off
MIN_ROWS = 400                   # don't fit an allele with too little data
SEED = 0

AAS = "ACDEFGHIKLMNPQRSTVWY"

# BLOSUM62, rows ordered as AAS above.
BLOSUM62 = np.array([
    [4,0,-2,-1,-2,0,-2,-1,-1,-1,-1,-2,-1,-1,-1,1,0,0,-3,-2],
    [0,9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2],
    [-2,-3,6,2,-3,-1,-1,-3,-1,-4,-3,1,-1,0,-2,0,-1,-3,-4,-3],
    [-1,-4,2,5,-3,-2,0,-3,1,-3,-2,0,-1,2,0,0,-1,-2,-3,-2],
    [-2,-2,-3,-3,6,-3,-1,0,-3,0,0,-3,-4,-3,-3,-2,-2,-1,1,3],
    [0,-3,-1,-2,-3,6,-2,-4,-2,-4,-3,0,-2,-2,-2,0,-2,-3,-2,-3],
    [-2,-3,-1,0,-1,-2,8,-3,-1,-3,-2,1,-2,0,0,-1,-2,-3,-2,2],
    [-1,-1,-3,-3,0,-4,-3,4,-3,2,1,-3,-3,-3,-3,-2,-1,3,-3,-1],
    [-1,-3,-1,1,-3,-2,-1,-3,5,-2,-1,0,-1,1,2,0,-1,-2,-3,-2],
    [-1,-1,-4,-3,0,-4,-3,2,-2,4,2,-3,-3,-2,-2,-2,-1,1,-2,-1],
    [-1,-1,-3,-2,0,-3,-2,1,-1,2,5,-2,-2,0,-1,-1,-1,1,-1,-1],
    [-2,-3,1,0,-3,0,1,-3,0,-3,-2,6,-2,0,0,1,0,-3,-4,-2],
    [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2,7,-1,-2,-1,-1,-2,-4,-3],
    [-1,-3,0,2,-3,-2,0,-3,1,-2,0,0,-1,5,1,0,-1,-2,-2,-1],
    [-1,-3,-2,0,-3,-2,0,-3,2,-2,-1,0,-2,1,5,-1,-1,-3,-3,-2],
    [1,-1,0,0,-2,0,-1,-2,0,-2,-1,1,-1,0,-1,4,1,-2,-3,-2],
    [0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1,0,-1,-1,-1,1,5,0,-2,-2],
    [0,-1,-3,-2,-1,-3,-3,3,-2,1,1,-3,-2,-2,-3,-2,0,4,-3,-1],
    [-3,-2,-4,-3,1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11,2],
    [-2,-2,-3,-2,3,-3,2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1,2,7],
], dtype=np.float32) / 4.0

AA_IDX = {a: i for i, a in enumerate(AAS)}


def encode(peptides, length):
    """BLOSUM62 encoding -> (n, length*20) float matrix."""
    X = np.zeros((len(peptides), length * 20), dtype=np.float32)
    for i, p in enumerate(peptides):
        for j, aa in enumerate(p[:length]):
            idx = AA_IDX.get(aa)
            if idx is not None:
                X[i, j * 20:(j + 1) * 20] = BLOSUM62[idx]
    return X


def log50k(nm):
    """IC50 in nM -> [0,1] affinity score. Higher = stronger binder."""
    return np.clip(1.0 - np.log(np.clip(nm, 1e-3, 5e4)) / np.log(5e4), 0.0, 1.0)


def inv_log50k(score):
    """Inverse transform: [0,1] score -> IC50 nM."""
    return np.exp((1.0 - np.clip(score, 0, 1)) * np.log(5e4))


def load():
    d = dataio.read("binding_affinity")
    d = d.dropna(subset=["peptide", "allele", "affinity_nM", "mhc_class"])
    # One measurement per peptide/allele: geometric mean across studies.
    d["log_nm"] = np.log(d.affinity_nM)
    g = (d.groupby(["allele", "mhc_class", "peptide", "length"], as_index=False)
           .agg(log_nm=("log_nm", "mean"), n_meas=("log_nm", "size")))
    g["affinity_nM"] = np.exp(g.log_nm)
    g["y"] = log50k(g.affinity_nM)
    return g


def split_by_peptide(df, frac=0.2):
    """Peptide-disjoint split so the same peptide never spans train and test."""
    peps = np.asarray(df.peptide.unique(), dtype=object)
    rng = np.random.default_rng(SEED)
    rng.shuffle(peps)
    test = set(peps[:max(1, int(len(peps) * frac))])
    mask = df.peptide.isin(test)
    return df[~mask], df[mask]


def train_all():
    df = load()
    report = {}

    for (allele, mclass), sub in df.groupby(["allele", "mhc_class"]):
        L = CANON_LEN.get(mclass)
        sub = sub[sub.length == L]
        if len(sub) < MIN_ROWS:
            print(f"  {allele:16s} skipped ({len(sub)} rows of length {L})")
            continue

        tr, te = split_by_peptide(sub)
        Xtr, ytr = encode(tr.peptide.values, L), tr.y.values
        Xte, yte = encode(te.peptide.values, L), te.y.values

        m = HistGradientBoostingRegressor(
            max_iter=400, learning_rate=0.08, max_leaf_nodes=31,
            l2_regularization=1.0, early_stopping=True,
            validation_fraction=0.1, random_state=SEED)
        m.fit(Xtr, ytr)

        pred = m.predict(Xte)
        rho = spearmanr(yte, pred).statistic
        true_bind = (inv_log50k(yte) < BINDER_NM).astype(int)
        auc = (roc_auc_score(true_bind, pred)
               if 0 < true_bind.sum() < len(true_bind) else float("nan"))

        report[allele] = {
            "mhc_class": mclass, "peptide_length": L,
            "n_train": int(len(tr)), "n_test": int(len(te)),
            "spearman": round(float(rho), 3), "auc": round(float(auc), 3),
            "test_binder_frac": round(float(true_bind.mean()), 3),
        }
        with open(MODELS / f"{allele.replace('*','_').replace(':','_')}.pkl", "wb") as f:
            pickle.dump({"model": m, "length": L, "allele": allele,
                         "mhc_class": mclass}, f)
        print(f"  {allele:16s} class {mclass:2s} L={L:2d}  "
              f"train={len(tr):6d} test={len(te):5d}  "
              f"Spearman={rho:.3f}  AUC={auc:.3f}")

    with open(MODELS / "report.json", "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n  models + report.json -> {MODELS}")
    return report


if __name__ == "__main__":
    print("Training per-allele peptide-MHC binding models")
    train_all()
