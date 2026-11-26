#!/usr/bin/env python3
"""
Read immunoforge output and turn it into simulator input.

WHY THIS IS THE PREFERRED PATH
------------------------------
immune_sim ships its own IEDB-trained binding model (model.py), which reaches
a mean held-out ROC-AUC of 0.916 for class I. That is respectable but it is
NOT state of the art. The immunoforge in this repository runs eleven
published predictors - NetMHCpan 4.2, MixMHCpred, MHCflurry 2, BigMHC, PRIME,
NetMHCIIpan 4.3, MixMHC2pred and others - and fuses them with population-
weighted HLA coverage. Where an engine run exists, its consensus is a better
input than our own model, so this adapter is the default and model.py is the
fallback for when the engine has not been run.

The engine reports PERCENTILE RANK against random human peptides rather than
nM. That is a better scale for comparing across alleles, so we keep it and
convert to an immunodominance weight directly.

Usage
-----
    python3 engine_adapter.py <run_dir_or_csv>            # inspect
    python3 simulate.py --engine <run_dir_or_csv> ...     # simulate from it
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

CLASS_II_PREFIX = ("DRB", "DRA", "DQ", "DP", "H2-I", "H-2I")

# Engine percentile-rank conventions (consensus.py defaults).
STRONG_RANK = {"I": 0.5, "II": 2.0}
WEAK_RANK = {"I": 2.0, "II": 5.0}


def mhc_class(allele):
    a = str(allele).upper().replace("HLA-", "").replace("*", "").replace("_", "")
    return "II" if a.startswith(CLASS_II_PREFIX) else "I"


def find_tables(path):
    """Accept a run directory or a CSV directly. Prefer the pair table."""
    p = Path(path)
    if p.is_file():
        return p, None
    if not p.is_dir():
        raise SystemExit(f"not found: {p}")
    pairs = sorted(p.glob("*_pairs.csv")) or sorted(p.rglob("*_pairs.csv"))
    peps = (sorted(p.glob("*_peptides.csv")) or sorted(p.glob("*_candidates.csv"))
            or sorted(p.rglob("*_peptides.csv")))
    if not pairs and not peps:
        raise SystemExit(
            f"no engine tables under {p}\n"
            f"expected *_pairs.csv or *_peptides.csv from an engine run")
    return (pairs[-1] if pairs else None), (peps[-1] if peps else None)


def _col(df, *names):
    for n in names:
        if n in df.columns:
            return n
    return None


def from_pairs(df):
    """One row per peptide x allele - the richest form."""
    pep = _col(df, "Peptide")
    alle = _col(df, "Allele", "Best_Allele")
    rank = _col(df, "Bind_Consensus_%Rank", "Bind_Consensus_PopW",
                "Bind_Consensus_Best")
    if not (pep and alle and rank):
        return None
    out = pd.DataFrame({
        "peptide": df[pep].astype(str).str.strip().str.upper(),
        "allele": df[alle].astype(str).str.strip(),
        "bind_rank": pd.to_numeric(df[rank], errors="coerce"),
    })
    for src, dst in (("Imm_TopPct_Consensus", "imm_rank"),
                     ("Imm_TopPct_Best", "imm_rank"),
                     ("Agree_Fraction", "agreement"),
                     ("Agree_Fraction_PopW", "agreement"),
                     ("Coverage_Binder", "coverage"),
                     ("Priority_Score", "priority"),
                     ("Start", "start")):
        if src in df.columns and dst not in out.columns:
            out[dst] = pd.to_numeric(df[src], errors="coerce")
    return out


def from_peptides(df):
    """One row per peptide; allele comes from Best_Allele."""
    pep = _col(df, "Peptide")
    alle = _col(df, "Best_Allele", "Allele")
    rank = _col(df, "Bind_Consensus_PopW", "Bind_Consensus_Best",
                "Bind_Consensus_%Rank")
    if not (pep and alle and rank):
        return None
    out = pd.DataFrame({
        "peptide": df[pep].astype(str).str.strip().str.upper(),
        "allele": df[alle].astype(str).str.strip(),
        "bind_rank": pd.to_numeric(df[rank], errors="coerce"),
    })
    for src, dst in (("Imm_TopPct_Best", "imm_rank"),
                     ("Agree_Fraction_PopW", "agreement"),
                     ("Coverage_Binder", "coverage"),
                     ("Priority_Score", "priority"),
                     ("Start", "start")):
        if src in df.columns:
            out[dst] = pd.to_numeric(df[src], errors="coerce")
    return out


def load(path, top_i=4, top_ii=3, use_priority=True):
    """
    Return the epitope frame the simulator expects:
    allele, mhc_class, peptide, start, pred_IC50_nM, immuno_p, weight.
    """
    pairs_f, peps_f = find_tables(path)
    tab = None
    if pairs_f is not None:
        tab = from_pairs(pd.read_csv(pairs_f))
    if tab is None and peps_f is not None:
        tab = from_peptides(pd.read_csv(peps_f))
    if tab is None:
        raise SystemExit("could not find the expected columns in the engine output")

    tab = tab.dropna(subset=["peptide", "allele", "bind_rank"])
    tab = tab[tab.peptide.str.len().between(6, 30)]
    if tab.empty:
        raise SystemExit("no usable rows in the engine output")

    tab["mhc_class"] = tab.allele.map(mhc_class)
    tab["length"] = tab.peptide.str.len()
    if "start" not in tab.columns:
        tab["start"] = np.nan

    # Percentile rank -> binding strength on [0, 1]. Rank 0.01 % is a top
    # binder, the weak-binder line scores about 0.3, rank 100 scores 0.
    def bind_w(row):
        r = float(np.clip(row.bind_rank, 1e-3, 100.0))
        return float(np.clip(1.0 - np.log(r / 1e-3) / np.log(100.0 / 1e-3), 0, 1))
    tab["bind_w"] = tab.apply(bind_w, axis=1)

    # Immunogenicity: the engine reports a top-percentile score, lower better.
    if "imm_rank" in tab.columns and tab.imm_rank.notna().any():
        ir = tab.imm_rank.clip(lower=0.1, upper=100.0)
        tab["immuno_p"] = np.clip(1.0 - np.log(ir / 0.1) / np.log(100.0 / 0.1), 0, 1)
    else:
        tab["immuno_p"] = np.nan

    # Priority_Score already fuses coverage, agreement, binding and
    # immunogenicity, so prefer it when present.
    if use_priority and "priority" in tab.columns and tab.priority.notna().any():
        w = tab.priority.clip(lower=0) / 100.0
    else:
        w = tab.bind_w * np.where(tab.immuno_p.notna(), tab.immuno_p, 1.0)
    tab["weight"] = w

    # A nominal nM, for display and for the C-ImmSim bit-string mapping only.
    tab["pred_IC50_nM"] = np.exp((1.0 - tab.bind_w) * np.log(5e4)).round(1)
    tab["call"] = np.where(
        tab.bind_rank <= tab.mhc_class.map(STRONG_RANK), "strong",
        np.where(tab.bind_rank <= tab.mhc_class.map(WEAK_RANK), "weak", "non-binder"))

    tab = tab[tab.call != "non-binder"]
    if tab.empty:
        raise SystemExit("engine run contains no binders at its own cut-offs")

    keep = []
    for cls, n in (("I", top_i), ("II", top_ii)):
        sub = tab[tab.mhc_class == cls]
        if not sub.empty:
            keep.append(sub.sort_values("weight", ascending=False)
                           .groupby("allele").head(n))
    tab = pd.concat(keep, ignore_index=True)
    tab["weight"] = tab.weight / tab.weight.max()
    return tab.sort_values(["mhc_class", "weight"], ascending=[True, False])


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    t = load(sys.argv[1])
    print(f"{len(t)} epitope clones from engine output\n")
    cols = [c for c in ["allele", "mhc_class", "peptide", "start", "bind_rank",
                        "immuno_p", "coverage", "priority", "weight", "call"]
            if c in t.columns]
    print(t[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))
