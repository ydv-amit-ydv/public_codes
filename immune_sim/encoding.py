#!/usr/bin/env python3
"""
Peptide encoding for the advanced binding models.

Two ideas the simple flat encoding in `model.py` misses:

1. **Binding register.** An MHC class II groove is open at both ends, so a
   13-25mer is bound by a 9-residue CORE whose position in the peptide is not
   known in advance. Encoding a fixed 15-mer left-aligned therefore smears
   nine anchor positions across fifteen columns and the motif never sharpens.
   The fix (NN-align / NetMHCIIpan) is to encode a 9-mer core at a given
   OFFSET, plus summary features of the peptide flanking regions, and to let
   the model choose the offset.

2. **Variable length is data, not noise.** `model.py` keeps only the canonical
   length and discards 24% of class I rows (10- and 11-mers) and 23% of
   class II rows. A core at an offset handles any length >= 9.

Feature block per (peptide, offset):

    core BLOSUM62            9 x 20 = 180
    core physicochemical     9 x  5 =  45
    N-flank BLOSUM mean           20
    C-flank BLOSUM mean           20
    flank lengths + peptide len    3
                                 ---
                                 268
"""
import numpy as np

AAS = "ACDEFGHIKLMNPQRSTVWY"
AA_IDX = {a: i for i, a in enumerate(AAS)}
CORE = 9
FLANK = 3          # residues of peptide flanking region summarised each side

# BLOSUM62, rows ordered as AAS, scaled to ~unit range.
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

# Five orthogonal physicochemical scales, each standardised across the 20
# residues. Trees can only split on what is in the columns; hydrophobicity and
# volume are what MHC pockets actually select on, and BLOSUM only encodes them
# implicitly through substitution frequency.
_HYDRO = dict(zip("AVLIPFMWGSTCYNQDEKRH",      # Kyte-Doolittle
    [1.8,4.2,3.8,4.5,-1.6,2.8,1.9,-0.9,-0.4,-0.8,
     -0.7,2.5,-1.3,-3.5,-3.5,-3.5,-3.5,-3.9,-4.5,-3.2]))
_VOL = dict(zip("GASCDPNTEVQHMILKRFYW",         # residue volume A^3
    [60.1,88.6,89.0,108.5,111.1,112.7,114.1,116.1,138.4,140.0,
     143.8,153.2,162.9,166.7,166.7,168.6,173.4,189.9,193.6,227.8]))
_CHARGE = {a: 0.0 for a in AAS}
_CHARGE.update({"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.5})
_POLAR = dict(zip("ARNDCQEGHILKMFPSTWYV",       # Grantham polarity
    [8.1,10.5,11.6,13.0,5.5,10.5,12.3,9.0,10.4,5.2,
     4.9,11.3,5.7,5.2,8.0,9.2,8.6,5.4,6.2,5.9]))
_HELIX = dict(zip("AELMQKRHVIWFYCTGNPSD",       # Chou-Fasman helix propensity
    [1.42,1.51,1.21,1.45,1.11,1.16,0.98,1.00,1.06,1.08,
     1.08,1.13,0.69,0.70,0.83,0.57,0.67,0.57,0.77,1.01]))


def _phys_matrix():
    M = np.zeros((20, 5), dtype=np.float32)
    for i, a in enumerate(AAS):
        M[i] = (_HYDRO.get(a, 0.0), _VOL.get(a, 140.0),
                _CHARGE.get(a, 0.0), _POLAR.get(a, 9.0), _HELIX.get(a, 1.0))
    M -= M.mean(0)
    sd = M.std(0); sd[sd == 0] = 1.0
    return M / sd


PHYS = _phys_matrix()
NPHYS = PHYS.shape[1]
FEAT = CORE * 20 + CORE * NPHYS + 20 + 20 + 3


def n_offsets(length, core=CORE):
    """How many 9-mer registers a peptide of this length admits."""
    return max(1, length - core + 1)


def encode_cores(peptides, offsets, core=CORE):
    """
    Encode one (peptide, offset) pair per row.

    `offsets` is the 0-based start of the core in each peptide. A peptide
    shorter than `core` is left-padded with zeros, which the model reads as
    "no residue" rather than as any particular amino acid.
    """
    n = len(peptides)
    X = np.zeros((n, FEAT), dtype=np.float32)
    b_end = core * 20
    p_end = b_end + core * NPHYS
    for i, (pep, o) in enumerate(zip(peptides, offsets)):
        L = len(pep)
        o = int(min(max(o, 0), max(L - core, 0)))
        for j in range(core):
            k = o + j
            if k >= L:
                break
            idx = AA_IDX.get(pep[k])
            if idx is None:
                continue
            X[i, j * 20:(j + 1) * 20] = BLOSUM62[idx]
            X[i, b_end + j * NPHYS:b_end + (j + 1) * NPHYS] = PHYS[idx]

        # peptide flanking regions: mean BLOSUM row over up to FLANK residues
        for side, (lo, hi) in enumerate(((max(0, o - FLANK), o),
                                         (o + core, min(L, o + core + FLANK)))):
            acc = np.zeros(20, dtype=np.float32)
            cnt = 0
            for k in range(lo, hi):
                idx = AA_IDX.get(pep[k])
                if idx is not None:
                    acc += BLOSUM62[idx]
                    cnt += 1
            if cnt:
                X[i, p_end + side * 20:p_end + (side + 1) * 20] = acc / cnt

        X[i, p_end + 40] = min(o, 9) / 9.0                  # N-flank length
        X[i, p_end + 41] = min(max(L - o - core, 0), 9) / 9.0  # C-flank length
        X[i, p_end + 42] = min(L, 25) / 25.0                # peptide length
    return X


def expand_all_offsets(peptides, core=CORE):
    """
    Every register of every peptide, for scanning at prediction time.

    Returns (row_index, offset, X) where row_index maps each encoded row back
    to its peptide, so a max/soft-max pool over registers is a grouped
    reduction on `row_index`.
    """
    idx, offs = [], []
    for i, pep in enumerate(peptides):
        for o in range(n_offsets(len(pep), core)):
            idx.append(i)
            offs.append(o)
    idx = np.asarray(idx, dtype=np.int64)
    offs = np.asarray(offs, dtype=np.int64)
    X = encode_cores([peptides[i] for i in idx], offs, core)
    return idx, offs, X


def pool_max(idx, scores, n):
    """Per-peptide best register, and which register won."""
    best = np.full(n, -np.inf, dtype=np.float64)
    arg = np.zeros(n, dtype=np.int64)
    for r, (i, s) in enumerate(zip(idx, scores)):
        if s > best[i]:
            best[i] = s
            arg[i] = r
    return best, arg


# --------------------------------------------------------------------------
# allele representation
# --------------------------------------------------------------------------
import re as _re


class AlleleEncoder:
    """
    Hierarchical allele features, replacing a flat one-hot.

    A one-hot over training alleles cannot transfer to an allele it has never
    seen - there is no column to put it in - and a fallback to index 0 is
    worse than useless: it scores the new allele as if it were whichever
    allele happens to be first. Measured on mass-spectrometry eluted ligands,
    that is exactly what happened. On the nine alleles in the training set the
    model reached ROC-AUC 0.90-0.98; on alleles outside it, 0.54-0.76.

    So the allele is encoded at three levels of specificity:

        locus        HLA-A, HLA-B, HLA-C, HLA-DRB1, H2-K, ...
        group        the first field, e.g. B*44 - shared by B*44:02 and B*44:03
        allele       the exact allele

    A trained allele lights up all three. An unseen allele lights up the two
    coarser ones and degrades to its group's or locus's average motif instead
    of to an arbitrary allele's. Groups matter because MHC alleles within a
    first-field group share most of their binding pocket; this is the cheap,
    offline approximation to a pseudo-sequence encoding.
    """

    def __init__(self, alleles):
        self.alleles = list(dict.fromkeys(alleles))
        self.loci = sorted({self.locus(a) for a in self.alleles})
        self.groups = sorted({self.group(a) for a in self.alleles})
        self._li = {k: i for i, k in enumerate(self.loci)}
        self._gi = {k: i for i, k in enumerate(self.groups)}
        self._ai = {k: i for i, k in enumerate(self.alleles)}
        self.dim = len(self.loci) + len(self.groups) + len(self.alleles)

    @staticmethod
    def locus(a):
        a = str(a)
        m = _re.match(r"^([A-Za-z0-9]+-[A-Za-z0-9]+)", a.replace("*", "-"))
        if m:
            return m.group(1).split("-")[0] + "-" + m.group(1).split("-")[1]
        return a.split("*")[0] or "unknown"

    @staticmethod
    def group(a):
        a = str(a)
        if "*" in a:
            head, rest = a.split("*", 1)
            return f"{head}*{rest.split(':')[0]}"
        return a

    def transform(self, alleles):
        A = np.zeros((len(alleles), self.dim), dtype=np.float32)
        off_g = len(self.loci)
        off_a = off_g + len(self.groups)
        for i, a in enumerate(alleles):
            li = self._li.get(self.locus(a))
            if li is not None:
                A[i, li] = 1.0
            gi = self._gi.get(self.group(a))
            if gi is not None:
                A[i, off_g + gi] = 1.0
            ai = self._ai.get(a)
            if ai is not None:
                A[i, off_a + ai] = 1.0
        return A

    def coverage(self, allele):
        """Which levels an allele matches: for reporting, not for scoring."""
        return {"locus": self.locus(allele) in self._li,
                "group": self.group(allele) in self._gi,
                "allele": allele in self._ai}
