"""Codon usage models for human + mouse (Kazusa tables via python_codon_tables), CAI and optional tAI."""
from __future__ import annotations

import csv
import math

import python_codon_tables as pct

BASES = "UCAG"
AA_STRING = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
CODON_TO_AA = dict(zip(CODONS, AA_STRING))
AA_TO_CODONS: dict[str, list[str]] = {}
for _c, _a in CODON_TO_AA.items():
    AA_TO_CODONS.setdefault(_a, []).append(_c)

TABLE_NAMES = {"human": "h_sapiens_9606", "mouse": "m_musculus_10090"}
NO_SCORE_AA = {"M", "W", "*"}   # single-codon aa / stop are excluded from CAI


def translate(rna: str) -> str:
    return "".join(CODON_TO_AA.get(rna[i:i + 3], "X") for i in range(0, len(rna) - len(rna) % 3, 3))


def usage_table(species: str) -> dict[str, float]:
    """Per-codon frequency *within its amino-acid family* (RNA alphabet)."""
    t = pct.get_codons_table(TABLE_NAMES[species])
    out = {}
    for aa, d in t.items():
        for dna, f in d.items():
            out[dna.replace("T", "U")] = f
    return out


# dos Reis et al. 2004 wobble penalties: codon 3rd base -> [(anticodon pos34, s)]
_TAI_S = {"U": [("A", 0.0), ("G", 0.41)], "C": [("G", 0.0), ("A", 0.28)],
          "A": [("U", 0.0), ("A", 0.9999)], "G": [("C", 0.0), ("U", 0.68)]}
_COMP = {"A": "U", "U": "A", "G": "C", "C": "G"}


def _load_trna(path: str) -> dict[str, dict[str, float]]:
    out: dict[str, dict[str, float]] = {}
    with open(path) as fh:
        for row in csv.DictReader(fh):
            out.setdefault(row["species"].strip().lower(), {})[row["anticodon"].strip().upper().replace("T", "U")] = float(row["copies"])
    return out


class CodonModel:
    """Relative adaptiveness w per codon per species; log-weights for fast scoring."""

    def __init__(self, species=("human", "mouse"), metric="cai", trna_file=None, min_freq=0.0):
        self.species = tuple(species)
        self.metric = metric
        self.w: dict[str, dict[str, float]] = {}
        if metric == "tai" and not trna_file:
            raise ValueError("codon_metric='tai' needs trna_copy_file (CSV: species,anticodon,copies)")
        trna = _load_trna(trna_file) if metric == "tai" else None
        for sp in self.species:
            self.w[sp] = self._cai_w(sp) if metric == "cai" else self._tai_w(sp, trna[sp])
        self.min_freq = min_freq
        self.freq = {sp: usage_table(sp) for sp in self.species}

    @staticmethod
    def _cai_w(sp):
        f = usage_table(sp)
        w = {}
        for aa, cods in AA_TO_CODONS.items():
            m = max(f.get(c, 0.0) for c in cods) or 1.0
            for c in cods:
                w[c] = max(f.get(c, 0.0) / m, 0.01)
        return w

    @staticmethod
    def _tai_w(sp, copies):
        raw = {}
        for c in CODONS:
            if CODON_TO_AA[c] == "*":
                continue
            core = "".join(_COMP[b] for b in reversed(c[:2]))  # anticodon pos35-36 pairs codon pos2,1
            tot = 0.0
            for n34, s in _TAI_S[c[2]]:
                tot += (1 - s) * copies.get(n34 + core, 0.0)
            raw[c] = tot
        mx = max(raw.values()) or 1.0
        pos = [v / mx for v in raw.values() if v > 0]
        gm = math.exp(sum(math.log(v) for v in pos) / len(pos))
        w = {c: (v / mx if v > 0 else gm) for c, v in raw.items()}
        w.update({c: 1.0 for c in CODONS if CODON_TO_AA[c] == "*"})
        return w

    def options(self, aa: str) -> list[str]:
        """Synonymous codons allowed for aa (drops codons rarer than min_freq in any species)."""
        cods = AA_TO_CODONS[aa]
        if self.min_freq <= 0 or len(cods) == 1:
            return list(cods)
        keep = [c for c in cods if all(self.freq[sp].get(c, 0) >= self.min_freq for sp in self.species)]
        return keep or list(cods)

    def logw(self, species: str, codon: str) -> float:
        return math.log(self.w[species][codon])

    def index(self, codons: list[str], species: str) -> float:
        """CAI/tAI = geometric mean of w over scored codons."""
        vals = [math.log(self.w[species][c]) for c in codons if CODON_TO_AA[c] not in NO_SCORE_AA]
        return math.exp(sum(vals) / len(vals)) if vals else 1.0
