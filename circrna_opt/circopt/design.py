"""Design space: which codons may change, and how a genome (codon choices) becomes a sequence."""
from __future__ import annotations

from .codon import CODON_TO_AA, CodonModel, translate
from .config import Params


class Design:
    """Holds the fixed scaffold (IRES/UTRs/spacers) and the free synonymous codon slots.

    A *genome* is a list of ints: one index into `self.options[k]` per free codon slot.
    """

    def __init__(self, p: Params, model: CodonModel | None = None):
        self.p = p
        rna = p.rna.upper().replace("T", "U")
        if not rna:
            raise ValueError("params.rna (the full circular molecule) is required")
        bad = set(rna) - set("ACGU")
        if bad:
            raise ValueError(f"rna contains non-ACGU characters: {sorted(bad)}")
        self.rna = rna
        self.n = len(rna)
        self.model = model or CodonModel(
            species=(p.primary_species, p.secondary_species),
            metric=p.codon_metric, trna_file=p.trna_copy_file, min_freq=p.min_codon_freq)

        self.cds_start, self.cds_end = self._locate_cds()
        if (self.cds_end - self.cds_start) % 3:
            raise ValueError("CDS length is not a multiple of 3")
        self.cds = (self.cds_start, self.cds_end)
        self.protein = translate(rna[self.cds_start:self.cds_end])
        if p.protein:
            want = p.protein.upper().strip()
            got = self.protein[:-1] if self.protein.endswith("*") and not want.endswith("*") else self.protein
            if got != want:
                raise ValueError("params.protein does not match the CDS translation of params.rna")

        # positions locked by the user (IRES, UTRs, spacers) or outside the CDS
        self.locked = [True] * self.n
        for i in range(self.cds_start, self.cds_end):
            self.locked[i] = False
        for r in p.protected:
            for i in range(max(0, r["start"]), min(self.n, r["end"])):
                self.locked[i] = True

        self.slots: list[int] = []        # codon start positions that may change
        self.options: list[list[str]] = []
        for pos in range(self.cds_start, self.cds_end, 3):
            codon = rna[pos:pos + 3]
            aa = CODON_TO_AA.get(codon, "X")
            if aa == "X":
                raise ValueError(f"bad codon {codon} at {pos}")
            if any(self.locked[pos + k] for k in range(3)):
                continue
            opts = [c for c in self.model.options(aa)]
            if codon not in opts:
                opts.append(codon)        # never lose the sequence the user gave us
            if len(opts) > 1:
                self.slots.append(pos)
                self.options.append(sorted(opts))
        # designable mask for motif / repeat accounting
        self.free_mask = [False] * self.n
        for pos in self.slots:
            for k in range(3):
                self.free_mask[pos + k] = True

    def _locate_cds(self) -> tuple[int, int]:
        p = self.p
        if p.cds_start is not None and p.cds_end is not None:
            return int(p.cds_start), int(p.cds_end)
        if not p.protein:
            raise ValueError("give either cds_start+cds_end or protein (to locate the CDS)")
        want = p.protein.upper().strip().rstrip("*")
        for frame_start in range(self.n - 3 * len(want) + 1):
            if self.rna[frame_start:frame_start + 3] != "AUG":
                continue
            end = frame_start + 3 * len(want)
            if translate(self.rna[frame_start:end]) == want:
                stop = self.rna[end:end + 3]
                if CODON_TO_AA.get(stop) == "*":
                    end += 3
                return frame_start, end
        raise ValueError("could not locate the protein's CDS in rna; pass cds_start/cds_end")

    # ---- genome <-> sequence -------------------------------------------------
    def genome_of(self, rna: str | None = None) -> list[int]:
        rna = (rna or self.rna).upper().replace("T", "U")
        g = []
        for k, pos in enumerate(self.slots):
            codon = rna[pos:pos + 3]
            g.append(self.options[k].index(codon) if codon in self.options[k] else 0)
        return g

    def build(self, genome: list[int]) -> str:
        chars = list(self.rna)
        for k, pos in enumerate(self.slots):
            codon = self.options[k][genome[k]]
            chars[pos:pos + 3] = codon
        return "".join(chars)

    def cds_codons(self, seq: str) -> list[str]:
        return [seq[i:i + 3] for i in range(self.cds_start, self.cds_end, 3)]

    def summary(self) -> dict:
        return {"length": self.n, "cds": [self.cds_start, self.cds_end],
                "protein_len": len(self.protein.rstrip("*")),
                "free_codons": len(self.slots),
                "total_cds_codons": (self.cds_end - self.cds_start) // 3,
                "protected": [dict(r) for r in self.p.protected],
                "choices": sum(len(o) for o in self.options)}
