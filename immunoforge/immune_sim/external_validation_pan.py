#!/usr/bin/env python3
"""
Leave-one-pathogen-out comparison: `model.py` versus `model_pan.py`.

The only question that matters for a vaccine design tool is whether a model
ranks the epitopes of a pathogen it has never seen. So for each organism:

    remove it entirely -> train both models on the rest -> predict it

and test the DIFFERENCE with a paired DeLong test on the same cases, with
Benjamini-Hochberg correction across pathogens. A pan-allele core-scanning
model has more parameters and more data; this is where it has to earn them.

    python3 external_validation_pan.py
    python3 external_validation_pan.py --quick
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import roc_auc_score

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
warnings.filterwarnings("ignore")

from external_validation import (MIN_TEST, MIN_TRAIN, N_PATHOGENS,   # noqa: E402
                                 aa_composition_baseline, organism)
from model import (BINDER_NM, CANON_LEN, SEED, encode,               # noqa: E402
                   inv_log50k, load)
from model_pan import MAX_LEN, MIN_LEN, PanEnsemble, load_eluted                  # noqa: E402
from stats import (auc_ci_delong, bh_fdr, calibration,               # noqa: E402
                   delong_test, fmt_p)
import dataio                                            # noqa: E402

RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)


def tag_organism(df):
    raw = dataio.read("binding_affinity")
    raw["org"] = organism(raw)
    pep_org = (raw.groupby("peptide").org
                  .agg(lambda s: s.value_counts().index[0]))
    df = df.copy()
    df["org"] = df.peptide.map(pep_org)
    return df[df.org.notna() & (df.org != "unknown")]


def baseline_model(train, test, L):
    """`model.py` exactly: canonical length only, flat BLOSUM, per class."""
    tr = train[train.length == L]
    te = test[test.length == L]
    if len(tr) < MIN_TRAIN or len(te) < 30:
        return None, None
    m = HistGradientBoostingRegressor(
        max_iter=400, learning_rate=0.08, max_leaf_nodes=31,
        l2_regularization=1.0, early_stopping=True,
        validation_fraction=0.1, random_state=SEED)
    m.fit(encode(tr.peptide.values, L), tr.y.values)
    return m.predict(encode(te.peptide.values, L)), te


def run(quick=False, realign=4, eluted=False, out_name=None,
        classes="I,II"):
    df = tag_organism(load())
    rows = []

    for mclass in [c.strip() for c in classes.split(",")]:
        L = CANON_LEN[mclass]
        sub = df[(df.mhc_class == mclass)
                 & (df.length >= MIN_LEN[mclass])
                 & (df.length <= MAX_LEN[mclass])]
        counts = sub.groupby("org").size().sort_values(ascending=False)
        pathogens = [o for o in counts.index
                     if counts[o] >= MIN_TEST][:N_PATHOGENS]
        print(f"\nclass {mclass}: {len(pathogens)} pathogens, "
              f"{len(sub)} usable rows "
              f"({(sub.length == L).sum()} at the canonical length)")

        for org in pathogens:
            test = sub[sub.org == org]
            train = sub[sub.org != org]
            test = test[~test.peptide.isin(set(train.peptide))]
            if len(test) < MIN_TEST or len(train) < MIN_TRAIN:
                continue

            # --- pan model, on every length -----------------------------
            peps, ys, als = (train.peptide.values, train.y.values,
                             train.allele.values)
            sw = np.ones(len(train))
            alleles = sorted(train.allele.unique())
            if eluted:
                # Leakage control, and the whole comparison depends on it:
                # the eluted set is keyed by peptide only, with no organism
                # column, so the held-out pathogen cannot be excluded by
                # name. Exclude by SEQUENCE instead - every test peptide is
                # dropped from the eluted rows before they are trained on.
                el = load_eluted(mclass, exclude_peptides=set(test.peptide),
                                 seed=SEED)
                if el is not None and len(el):
                    peps = np.concatenate([peps, el.peptide.values])
                    ys = np.concatenate([ys, el.y.values])
                    als = np.concatenate([als, el.allele.values])
                    sw = np.concatenate([sw, el.w.values])
                    alleles = sorted(set(alleles) | set(el.allele.unique()))
            ens = PanEnsemble(alleles, realign=realign, quick=quick)
            ens.fit(peps, ys, als, verbose=False, sw=sw)
            p_pan = ens.predict(test.peptide.values, test.allele.values)
            y_pan = (inv_log50k(test.y.values) < BINDER_NM).astype(int)
            if len(np.unique(y_pan)) < 2:
                continue

            # --- old model, canonical length only -----------------------
            p_old, te_old = baseline_model(train, test, L)

            auc_pan, (lo, hi) = auc_ci_delong(y_pan, p_pan)
            rec = {"pathogen": org, "mhc_class": mclass,
                   "n_test_pan": int(len(test)),
                   "n_train_pan": int(len(train)),
                   "auc_pan": round(float(auc_pan), 3),
                   "pan_lo": round(lo, 3), "pan_hi": round(hi, 3),
                   "ece_pan": None, "auc_old": None, "auc_diff": None,
                   "p_pan_vs_old": None, "n_shared": None}

            # compare ONLY on the rows both models can score, or the
            # difference would partly be the extra lengths, not the model
            if p_old is not None:
                m = test.length.values == L
                y_s = y_pan[m]
                if len(np.unique(y_s)) > 1 and m.sum() >= 30:
                    _, _, diff, pv = delong_test(y_s, p_pan[m], p_old)
                    rec.update({
                        "n_shared": int(m.sum()),
                        "auc_old": round(float(roc_auc_score(y_s, p_old)), 3),
                        "auc_pan_shared": round(
                            float(roc_auc_score(y_s, p_pan[m])), 3),
                        "auc_diff": round(float(diff), 3),
                        "p_pan_vs_old": float(pv)})

            # composition baseline, on the pan model's own test rows
            try:
                b = aa_composition_baseline(train, test, L)
                rec["auc_comp"] = round(float(roc_auc_score(y_pan, b)), 3)
                _, _, d2, p2 = delong_test(y_pan, p_pan, b)
                rec["p_vs_comp"] = float(p2)
                rec["gain_vs_comp"] = round(float(d2), 3)
            except Exception:
                rec["auc_comp"] = rec["p_vs_comp"] = rec["gain_vs_comp"] = None

            rows.append(rec)
            d = rec.get("auc_diff")
            print(f"  {org[:30]:30s} {mclass:2s} n={len(test):5d}  "
                  f"pan {auc_pan:.3f} [{lo:.3f}-{hi:.3f}]"
                  + (f"   old {rec['auc_old']:.3f}  diff {d:+.3f} "
                     f"(p={fmt_p(rec['p_pan_vs_old'])})"
                     if d is not None else "   old n/a"), flush=True)

    if not rows:
        sys.exit("no pathogen had enough data")

    tab = pd.DataFrame(rows)
    for col, newcol in (("p_pan_vs_old", "q_pan_vs_old"),
                        ("p_vs_comp", "q_vs_comp")):
        ok = tab[col].notna()
        if ok.sum():
            adj, sig = bh_fdr(tab.loc[ok, col].values)
            tab.loc[ok, newcol] = adj.round(4)
            tab.loc[ok, newcol.replace("q_", "sig_")] = sig
    # Honour --out. This line was NOT patched when the flag was added - only
    # the message that claims where it wrote was - so a run with --out
    # printed the name it was asked for and silently overwrote the shipped
    # results file instead. That dropped class I from a file the Evidence
    # tab reads. Third time a run has clobbered a shipped artefact; the
    # lesson each time is the same, so: one name, used for both.
    dest = RESULTS / (out_name or "external_validation_pan.csv")
    tab.to_csv(dest, index=False)

    print("\n" + "=" * 78)
    print("PAN-ALLELE CORE-SCAN vs model.py, on pathogens neither has seen")
    print("=" * 78)
    for mclass in [c.strip() for c in classes.split(",")]:
        s = tab[tab.mhc_class == mclass]
        if s.empty:
            continue
        c = s[s.auc_diff.notna()]
        print(f"\n  class {mclass}: {len(s)} pathogens")
        print(f"    pan  median AUC {s.auc_pan.median():.3f}  "
              f"range {s.auc_pan.min():.3f}-{s.auc_pan.max():.3f}")
        if not c.empty:
            print(f"    old  median AUC {c.auc_old.median():.3f}  "
                  f"(on the shared canonical-length rows)")
            print(f"    diff median {c.auc_diff.median():+.3f}  "
                  f"pan better in {int((c.auc_diff > 0).sum())}/{len(c)}  "
                  f"significant after BH in "
                  f"{int(c.get('sig_pan_vs_old', pd.Series(dtype=bool)).sum())}"
                  f"/{len(c)}")
        g = s[s.gain_vs_comp.notna()]
        if not g.empty:
            print(f"    beats composition baseline in "
                  f"{int(g.get('sig_vs_comp', pd.Series(dtype=bool)).sum())}"
                  f"/{len(g)} after BH")
    print(f"\n  wrote {dest}")
    return tab


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--realign", type=int, default=4)
    ap.add_argument("--eluted", action="store_true",
                    help="also train on eluted ligands; every test peptide is "
                         "excluded from them first, or the comparison would "
                         "be measuring leakage")
    ap.add_argument("--classes", default="I,II",
                    help="which MHC classes to fold over")
    ap.add_argument("--out", dest="out_name", default=None,
                    help="result CSV name under results/")
    run(**vars(ap.parse_args()))
