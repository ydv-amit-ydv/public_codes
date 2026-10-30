"""Accept whatever the user pastes: protein, DNA or RNA, payload-only or a whole circle.

Everything here is about turning one text box into an unambiguous design problem, and saying
out loud what was assumed. Nothing guesses silently: every inference returns a note that the
caller shows the user.
"""
from __future__ import annotations

import re

from .codon import AA_TO_CODONS, CODON_TO_AA, translate

AA20 = set("ACDEFGHIKLMNPQRSTVWY")
# letters that can ONLY be amino acids (never a nucleotide code)
AA_ONLY = AA20 - set("ACGTUN")
NT = set("ACGU")
IUPAC_NT = set("ACGTURYSWKMBDHVN")


def clean(text: str) -> str:
    """Strip FASTA headers, numbers, whitespace and punctuation; upper-case the rest."""
    lines = [ln for ln in (text or "").splitlines() if not ln.strip().startswith((">", ";"))]
    s = "".join(lines).upper()
    return re.sub(r"[^A-Z*]", "", s)


def detect(seq: str) -> tuple[str, str]:
    """Return (kind, note) where kind is 'protein' | 'dna' | 'rna'.

    The ambiguity that matters: A, C, G and T are all valid amino-acid letters too, so a short
    peptide of only those letters is indistinguishable from DNA. Nucleotide wins, because that
    is overwhelmingly the commoner input -- and the caller can override it.
    """
    if not seq:
        raise ValueError("empty sequence")
    s = set(seq.rstrip("*"))
    if s & AA_ONLY:
        bad = s - AA20 - {"*"}
        if bad:
            raise ValueError(f"not a valid protein, DNA or RNA sequence: unexpected {sorted(bad)}")
        return "protein", ""
    if s <= {"A", "C", "G", "T", "U"}:
        has_t, has_u = "T" in s, "U" in s
        if has_t and has_u:
            raise ValueError("sequence mixes T and U -- is it DNA or RNA?")
        kind = "rna" if has_u else "dna"
        note = ""
        if s <= AA20 and len(seq) < 60:
            note = (f"read as {kind.upper()}; these letters are also valid amino acids, so set "
                    f"the sequence type explicitly if you meant a protein")
        return kind, note
    if s <= AA20 | {"*"}:
        return "protein", ""
    bad = s - AA20 - IUPAC_NT - {"*"}
    if bad:
        raise ValueError(f"unexpected characters in sequence: {sorted(bad)}")
    raise ValueError("ambiguous sequence codes (N/R/Y/...): give plain ACGT/ACGU or a protein")


def to_rna(seq: str) -> str:
    return seq.replace("T", "U")


def reverse_translate(protein: str, model, primary: str, secondary: str, w: float) -> str:
    """Build a starting CDS for a protein, using the best dual-species codon per residue.

    This is only a seed -- the optimiser rewrites it. Starting from the codon-optimal sequence
    rather than a random one means cycle 1 begins inside the region of interest.
    """
    import math
    out = []
    for aa in protein:
        if aa == "*":
            out.append("UAA")
            continue
        opts = model.options(aa)
        out.append(max(opts, key=lambda c: w * math.log(model.w[primary][c])
                       + (1 - w) * math.log(model.w[secondary][c])))
    return "".join(out)


def find_orfs(rna: str) -> list[tuple[int, int]]:
    """All AUG..stop ORFs as (start, end_exclusive_including_stop), longest first."""
    out = []
    for i in range(len(rna) - 2):
        if rna[i:i + 3] != "AUG":
            continue
        for j in range(i + 3, len(rna) - 2, 3):
            if CODON_TO_AA.get(rna[j:j + 3]) == "*":
                out.append((i, j + 3))
                break
    out.sort(key=lambda t: t[1] - t[0], reverse=True)
    return out


def looks_like_payload(rna: str) -> bool:
    """True if the sequence looks like a bare CDS: in frame, starts AUG, no internal stop."""
    if len(rna) < 6 or len(rna) % 3:
        return False
    if rna[:3] != "AUG":
        return False
    prot = translate(rna)
    return "*" not in prot[:-1]


def resolve(p, model) -> dict:
    """Turn `p.sequence` (+ optional scaffolds) into rna / protein / protected / notes.

    Returns a dict of fields to apply to the Params. `p.rna` set explicitly always wins, so
    the existing full-molecule workflow is untouched.
    """
    notes: list[str] = []
    seq = clean(p.sequence)
    if not seq:
        raise ValueError("no sequence given")

    kind = p.sequence_type
    if kind == "auto":
        kind, note = detect(seq)
        if note:
            notes.append(note)
        notes.append(f"sequence read as {kind.upper()} ({len(seq)} "
                     f"{'residues' if kind == 'protein' else 'nt'})")

    if kind == "protein":
        protein = seq.rstrip("*")
        bad = set(protein) - AA20
        if bad:
            raise ValueError(f"protein contains non-standard residues: {sorted(bad)}")
        cds = reverse_translate(protein, model, p.primary_species, p.secondary_species,
                                p.primary_weight)
        notes.append("no nucleotide sequence given, so the CDS was built by reverse translation "
                     "with the best dual-species codon per residue; the optimiser takes it from "
                     "there, and the 'original' column is that seed, not a real construct")
        role = "payload"
    else:
        rna = to_rna(seq)
        protein = None
        role = p.sequence_role
        if role == "auto":
            role = "payload" if looks_like_payload(rna) else "molecule"
            notes.append(f"treated as the {'CDS payload' if role == 'payload' else 'whole molecule'}")
        if role == "payload":
            if len(rna) % 3:
                raise ValueError(f"payload length {len(rna)} is not a multiple of 3 -- "
                                 "is this the whole molecule rather than just the CDS?")
            cds = rna
        else:
            orfs = find_orfs(rna)
            if not orfs:
                raise ValueError("no AUG..stop ORF found in the molecule; give the protein, "
                                 "or cds_start/cds_end")
            s, e = orfs[0]
            prot = translate(rna[s:e]).rstrip("*")
            notes.append(f"CDS taken as the longest ORF, {s}..{e} ({len(prot)} aa); "
                         "everything outside it is held fixed")
            return {"rna": rna, "protein": prot, "cds_start": s, "cds_end": e, "notes": notes}

    # payload (or reverse-translated protein): add a stop, then wrap in the scaffolds
    if p.add_stop and CODON_TO_AA.get(cds[-3:]) != "*":
        cds += "UAA"
        notes.append("added a UAA stop codon (set add_stop=false for an infinite-ORF design)")
    prot = translate(cds).rstrip("*")
    if "*" in translate(cds)[:-1]:
        raise ValueError("the payload has an internal stop codon -- check the reading frame")

    five = to_rna(clean(p.scaffold_5))
    three = to_rna(clean(p.scaffold_3))
    for name, sc in (("scaffold_5", five), ("scaffold_3", three)):
        if sc and set(sc) - NT:
            raise ValueError(f"{name} must be plain DNA/RNA (ACGT/ACGU)")
    rna = five + cds + three
    protected = []
    if five:
        protected.append({"name": "scaffold_5", "start": 0, "end": len(five)})
    if three:
        protected.append({"name": "scaffold_3", "start": len(rna) - len(three), "end": len(rna)})
    if not five and not three:
        notes.append("NOTE: only the payload was given, so the circle folded here is the CDS "
                     "alone. The real construct also contains the IRES and spacers, which will "
                     "change the fold -- re-run with the full molecule before you commit to a "
                     "sequence")
    return {"rna": rna, "protein": prot, "cds_start": len(five), "cds_end": len(five) + len(cds),
            "protected": protected + list(p.protected), "notes": notes}
