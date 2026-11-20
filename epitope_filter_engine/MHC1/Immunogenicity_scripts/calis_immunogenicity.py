"""IEDB Class I pMHC immunogenicity (Calis et al. 2013), implemented natively.

The model is small: each residue has a fixed immunogenicity value (a 20-entry scale), each of the
nine positions has a weight (longer peptides repeat the weight 0.30 for the extra middle positions),
and the score is the weighted sum over the positions that are not masked. Masked positions are the
MHC anchor residues, whose effect is already captured by binding. Higher = more immunogenic.

Reference: Calis JJ, Maybeno M, Greenbaum JA, Weiskopf D, De Silva AD, Sette A, Kesmir C, Peters B.
"Properties of MHC class I presented peptides that enhance immunogenicity." PLoS Comput Biol 2013;
9(10):e1003266. https://doi.org/10.1371/journal.pcbi.1003266
The constants below are those of IEDB's next-generation T Cell Class I package (v0.2-beta,
immunogenicity_predictor.py, Non-Profit Open Software License 3.0, (c) La Jolla Institute for
Immunology); this module reproduces its scores exactly (see test_engine.py).

Mask choices (as on the IEDB site):
  default    positions 1, 2 and the C-terminal residue
  by_allele  fixed positions per allele (table below); an allele not in the table falls back to default
  custom     positions you give, 1-based, e.g. "2,5,9"
"""
import numpy as np
import pandas as pd

SCALE = {"A": 0.127, "C": -0.175, "D": 0.072, "E": 0.325, "F": 0.380, "G": 0.110, "H": 0.105, "I": 0.432,
         "K": -0.700, "L": -0.036, "M": -0.570, "N": -0.021, "P": -0.036, "Q": -0.376, "R": 0.168, "S": -0.537,
         "T": 0.126, "V": 0.134, "W": 0.719, "Y": -0.012}
WEIGHT = [0.00, 0.00, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0.00]

BY_ALLELE = {"H-2-Db": "2,5,9", "H-2-Dd": "2,3,5", "H-2-Kb": "2,3,9", "H-2-Kd": "2,5,9", "H-2-Kk": "2,8,9",
             "H-2-Ld": "2,5,9", "HLA-A0101": "2,3,9", "HLA-A0201": "1,2,9", "HLA-A0202": "1,2,9",
             "HLA-A0203": "1,2,9", "HLA-A0206": "1,2,9", "HLA-A0211": "1,2,9", "HLA-A0301": "1,2,9",
             "HLA-A1101": "1,2,9", "HLA-A2301": "2,7,9", "HLA-A2402": "2,7,9", "HLA-A2601": "1,2,9",
             "HLA-A2902": "2,7,9", "HLA-A3001": "1,3,9", "HLA-A3002": "2,7,9", "HLA-A3101": "1,2,9",
             "HLA-A3201": "1,2,9", "HLA-A3301": "1,2,9", "HLA-A6801": "1,2,9", "HLA-A6802": "1,2,9",
             "HLA-A6901": "1,2,9", "HLA-B0702": "1,2,9", "HLA-B0801": "2,5,9", "HLA-B1501": "1,2,9",
             "HLA-B1502": "1,2,9", "HLA-B1801": "1,2,9", "HLA-B2705": "2,3,9", "HLA-B3501": "1,2,9",
             "HLA-B3901": "1,2,9", "HLA-B4001": "1,2,9", "HLA-B4002": "1,2,9", "HLA-B4402": "2,3,9",
             "HLA-B4403": "2,3,9", "HLA-B4501": "1,2,9", "HLA-B4601": "1,2,9", "HLA-B5101": "1,2,9",
             "HLA-B5301": "1,2,9", "HLA-B5401": "1,2,9", "HLA-B5701": "1,2,9", "HLA-B5801": "1,2,9"}

MASK_CHOICES = ("default", "by_allele", "custom")


def _weights(length, masked):
    """Position weights for a peptide of this length with the given 0-based masked positions zeroed."""
    w = np.array(WEIGHT[:5] + [0.30] * (length - 9) + WEIGHT[5:] if length > 9 else WEIGHT, dtype=float)
    # peptides shorter than 9 only have as many positions as residues (the original indexes the same way)
    w = w[:length]
    for m in masked:
        if 0 <= m < length:
            w[m] = 0.0
    return w


def allele_key(allele):
    """HLA-A*02:01 -> HLA-A0201, H2-Kb -> H-2-Kb (the table's spelling)."""
    return allele.replace("*", "").replace(":", "").replace("H2", "H-2") if allele else None


def parse_positions(text):
    try:
        pos = [int(x) for x in str(text).replace(" ", "").split(",") if x != ""]
    except ValueError:
        raise ValueError("Positions to mask must be numbers like 2,5,9")
    if not pos or any(p < 1 for p in pos):
        raise ValueError("Positions to mask must be 1 or higher, e.g. 2,5,9")
    return [p - 1 for p in pos]


def score_peptides(peptides, masked=None):
    """Scores for peptides (strings) with one mask. masked=None means default (1, 2, C-terminal);
    otherwise a list of 0-based positions. Returns a float array in input order."""
    out = np.empty(len(peptides))
    for length in sorted({len(p) for p in peptides}):
        idx = [i for i, p in enumerate(peptides) if len(p) == length]
        m = [0, 1, length - 1] if masked is None else masked
        w = _weights(length, m)
        vals = np.array([[SCALE[c] for c in peptides[i]] for i in idx])
        out[idx] = np.round((vals * w).sum(axis=1), 5)
    return out


def predict(peptides, alleles, mask_choice="default", positions=""):
    """-> (DataFrame[Peptide, Allele, IEDB_Immunogenicity], notes). Alleles are canonical (HLA-A*02:01)."""
    if mask_choice not in MASK_CHOICES:
        raise ValueError(f"mask_choice must be one of {MASK_CHOICES}")
    peps = sorted({p.upper() for p in peptides})
    bad = [p for p in peps if set(p) - set(SCALE)]
    if bad:
        raise ValueError(f"invalid residues in {bad[0]}")
    notes, frames = [], []
    if mask_choice == "by_allele":
        missing = []
        for a in alleles:
            pos = BY_ALLELE.get(allele_key(a))
            if pos is None:
                missing.append(a)
            masked = None if pos is None else parse_positions(pos)
            frames.append(pd.DataFrame({"Peptide": peps, "Allele": a, "IEDB_Immunogenicity": score_peptides(peps, masked)}))
        if missing:
            notes.append("IEDB immunogenicity has no allele-specific mask for " + ", ".join(missing)
                         + "; the default mask (1, 2, C-terminal) was used for them")
    else:
        masked = parse_positions(positions) if mask_choice == "custom" else None
        s = score_peptides(peps, masked)                       # allele-independent: same score for every allele
        for a in alleles:
            frames.append(pd.DataFrame({"Peptide": peps, "Allele": a, "IEDB_Immunogenicity": s}))
    return pd.concat(frames, ignore_index=True), notes
