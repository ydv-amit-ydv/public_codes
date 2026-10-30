"""DERNA integration (Gu, Qi & El-Kebir 2024, J Comput Biol; bioconda package `derna`).

DERNA designs a CDS from a protein by balancing MFE against CAI with dynamic programming, and
returns Pareto-optimal sequences. Two things to know before reading this module:

  * It optimises toward *minimal* MFE, i.e. the most structured sequence -- the opposite
    direction to this package's default `structure_goal="minimize_pairing"`. So it is used
    as (a) a high-CAI seed at low lambda and (b) an explicit comparison arm showing how much
    structure a stability-first design would carry.
  * `-c` (custom codon usage table) aborts on every documented-looking CSV layout, so DERNA
    runs with its built-in human table. Mouse/human codon usage differ marginally; the
    dual-species balancing is done by this package's optimiser, not by DERNA.

Only the codons at designable slots are adopted from DERNA output; protected regions
(IRES, UTRs, spacers) are never touched.
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import tempfile

from .design import Design

_RNA_LINE = re.compile(r"^(?:zuker cai rna|zuker rna|rna):\s*([ACGUX]+)", re.M)


def find_binary(p) -> str | None:
    cand = p.derna_bin or "derna"
    if os.path.sep in cand:
        return cand if os.path.exists(cand) else None
    return shutil.which(cand)


def _parse(text: str) -> list[str]:
    out = []
    for m in _RNA_LINE.finditer(text):
        s = m.group(1)
        if "X" in s:          # DERNA writes X where the MFE-only pass left a codon undetermined
            continue
        if s not in out:
            out.append(s)
    return out


def run_derna(design: Design, log=print) -> list[dict]:
    """Run DERNA at each configured lambda. Returns [{lambda, rna, full_seq, adopted}]."""
    p = design.p
    binary = find_binary(p)
    if not binary:
        log("DERNA: binary not found (conda install -c conda-forge -c bioconda derna) - skipping")
        return []
    protein = design.protein.rstrip("*")
    if len(protein) > p.derna_max_aa:
        log(f"DERNA: protein is {len(protein)} aa > derna_max_aa={p.derna_max_aa} - skipping "
            f"(DERNA is O(n^3); raise derna_max_aa to force it)")
        return []
    results = []
    with tempfile.TemporaryDirectory(prefix="derna_") as td:
        fa = os.path.join(td, "p.fasta")
        with open(fa, "w") as fh:
            fh.write(">protein\n" + protein + "\n")
        for lam in p.derna_lambdas:
            outp = os.path.join(td, f"o_{lam}.txt")
            cmd = [binary, "-i", fa, "-o", outp, "-m", "1", "-s", "2", "-l", str(lam)]
            log(f"DERNA: lambda={lam} ({len(protein)} aa, this is the slow step)")
            try:
                proc = subprocess.run(cmd, capture_output=True, text=True,
                                      timeout=p.derna_timeout_s, cwd=td)
            except subprocess.TimeoutExpired:
                log(f"DERNA: timed out after {p.derna_timeout_s}s at lambda={lam}")
                continue
            text = ""
            if os.path.exists(outp):
                with open(outp) as fh:
                    text = fh.read()
            text += "\n" + (proc.stdout or "")
            seqs = _parse(text)
            if not seqs:
                log(f"DERNA: no sequence parsed at lambda={lam} (rc={proc.returncode}) "
                    f"{(proc.stderr or '').strip()[:200]}")
                continue
            for rna in seqs:
                full, adopted = graft(design, rna)
                if full:
                    results.append({"lambda": lam, "cds_rna": rna, "full_seq": full,
                                    "adopted_codons": adopted})
    log(f"DERNA: produced {len(results)} seed(s)")
    return results


def graft(design: Design, cds_rna: str) -> tuple[str | None, int]:
    """Put a DERNA CDS into the molecule at designable slots only.

    DERNA emits the CDS without the stop codon, in the same reading frame as `design.protein`.
    """
    prot = design.protein.rstrip("*")
    if len(cds_rna) != 3 * len(prot):
        return None, 0
    chars = list(design.rna)
    adopted = 0
    for k, pos in enumerate(design.slots):
        ci = (pos - design.cds_start) // 3
        codon = cds_rna[3 * ci: 3 * ci + 3]
        if codon in design.options[k]:
            chars[pos:pos + 3] = codon
            adopted += 1
    return "".join(chars), adopted
