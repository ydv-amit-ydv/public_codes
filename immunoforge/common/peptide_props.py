"""Per-peptide details that help when ordering or handling a peptide. These are plain sequence calculations, not
predictions of anything biological: average mass, Kyte-Doolittle hydropathy (GRAVY), net charge, isoelectric point
and a few flags that matter for synthesis and solubility.
"""
import pandas as pd

MASS = {"G": 57.0519, "A": 71.0788, "S": 87.0782, "P": 97.1167, "V": 99.1326, "T": 101.1051, "C": 103.1388,
        "L": 113.1594, "I": 113.1594, "N": 114.1038, "D": 115.0886, "Q": 128.1307, "K": 128.1741, "E": 129.1155,
        "M": 131.1926, "H": 137.1411, "F": 147.1766, "R": 156.1875, "Y": 163.1760, "W": 186.2132}
WATER = 18.0153
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
      "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
# side-chain / terminal pKa (Lehninger-style set)
PKA_POS = {"K": 10.5, "R": 12.4, "H": 6.0}
PKA_NEG = {"D": 3.9, "E": 4.1, "C": 8.3, "Y": 10.1}
PKA_NTERM, PKA_CTERM = 9.69, 2.34
HYDROPHOBIC = set("AVILMFWC")


def charge(seq, ph=7.4):
    pos = 1 / (1 + 10 ** (ph - PKA_NTERM)) + sum(1 / (1 + 10 ** (ph - PKA_POS[c])) for c in seq if c in PKA_POS)
    neg = 1 / (1 + 10 ** (PKA_CTERM - ph)) + sum(1 / (1 + 10 ** (PKA_NEG[c] - ph)) for c in seq if c in PKA_NEG)
    return pos - neg


def isoelectric_point(seq):
    lo, hi = 0.0, 14.0
    for _ in range(40):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if charge(seq, mid) > 0 else (lo, mid)
    return (lo + hi) / 2


def longest_hydrophobic_run(seq):
    best = run = 0
    for c in seq:
        run = run + 1 if c in HYDROPHOBIC else 0
        best = max(best, run)
    return best


def notes(seq):
    n = []
    if "C" in seq: n.append("Cys: oxidises / dimerises")
    if "M" in seq: n.append("Met: oxidation-prone")
    if "W" in seq: n.append("Trp: oxidation-prone")
    if seq[0] in "QE": n.append("N-terminal Gln/Glu: pyroglutamate formation")
    if "DP" in seq or "DG" in seq: n.append("Asp-Pro / Asp-Gly: acid-labile or aspartimide-prone")
    if longest_hydrophobic_run(seq) >= 5: n.append("hydrophobic run >= 5: aggregation / low solubility")
    if sum(KD[c] for c in seq) / len(seq) > 1.0: n.append("GRAVY > 1: poorly soluble")
    return "; ".join(n)


def describe(peptides):
    """DataFrame indexed like `peptides` (unique sequences) with the detail columns."""
    rows = []
    for p in peptides:
        if set(p) - set(MASS):                       # X, U, B ... : no honest mass or charge
            rows.append({"Peptide": p, "Synthesis_Notes": "contains a non-standard residue: properties not calculated"})
            continue
        rows.append({"Peptide": p, "Mass_Da": round(sum(MASS[c] for c in p) + WATER, 2),
                     "GRAVY": round(sum(KD[c] for c in p) / len(p), 3), "Net_Charge_pH7.4": round(charge(p), 2),
                     "pI": round(isoelectric_point(p), 2), "Aromatic_Frac": round(sum(c in "FWY" for c in p) / len(p), 3),
                     "Hydrophobic_Frac": round(sum(c in HYDROPHOBIC for c in p) / len(p), 3),
                     "N_Cys": p.count("C"), "N_Met": p.count("M"), "N_Trp": p.count("W"), "Synthesis_Notes": notes(p)})
    return pd.DataFrame(rows)
