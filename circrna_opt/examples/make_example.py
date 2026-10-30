"""Build examples/example_params.json: a synthetic circRNA vaccine construct.

The construct is deliberately *bad* so a run has something to fix:
  * the CDS is written with low-usage synonymous codons (poor CAI in both species),
  * a perfect ~34 bp inverted repeat is planted inside the CDS using only synonymous
    codon choices, which folds into a long duplex (above the default 30 bp dsRNA threshold),
  * the region just after the start codon is made self-complementary, so the start
    context starts out paired.

The IRES block is a PLACEHOLDER, not a real IRES. Replace it with your own (CVB3, EMCV, ...)
and keep it in `protected` so it is never edited.
"""
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from circopt.codon import AA_TO_CODONS, CODON_TO_AA, translate  # noqa: E402

# A benign standard test protein (EGFP N-terminal fragment) standing in for an antigen.
PROTEIN = ("MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQ"
           "HDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNS")

rng = random.Random(7)
FREQ_RARE = {}   # pick a low-usage codon per aa to leave optimisation headroom
import python_codon_tables as pct
hum = {c.replace("T", "U"): f for aa, d in pct.get_codons_table("h_sapiens_9606").items()
       for c, f in d.items()}
for aa, cods in AA_TO_CODONS.items():
    FREQ_RARE[aa] = min(cods, key=lambda c: hum.get(c, 0))

cds = "".join(FREQ_RARE[a] for a in PROTEIN)

# plant a synonymous inverted repeat: take a 36 nt window, and rewrite a downstream window
# to its reverse complement wherever synonymous codons allow it
COMP = {"A": "U", "U": "A", "G": "C", "C": "G"}


def revcomp(s):
    return "".join(COMP[b] for b in reversed(s))


def nudge(cds, start_codon_i, target):
    """Rewrite codons from start_codon_i toward `target` nt, synonymously, best-effort."""
    out = list(cds)
    for k in range(len(target) // 3):
        ci = (start_codon_i + k) * 3
        aa = CODON_TO_AA[cds[ci:ci + 3]]
        want = target[3 * k:3 * k + 3]
        best = min(AA_TO_CODONS[aa], key=lambda c: sum(1 for x, y in zip(c, want) if x != y))
        out[ci:ci + 3] = best
    return "".join(out)


src_i = 20                      # codon index of the repeat's left arm
arm = cds[src_i * 3:src_i * 3 + 36]
cds = nudge(cds, 60, revcomp(arm))      # right arm ~120 nt downstream

# make the start context pairable: codons 2..7 complementary to the 5' side of the CDS
head = cds[0:18]
cds = nudge(cds, 8, revcomp(head))

assert translate(cds) == PROTEIN, "the example must stay synonymous"
cds += "UAA"

ires = ("GGGAGACCCAAGCUGGCUAGCGUUUAAACUUAAGCUUGGUACCGAGCUCGGAUCCACUAGUCCAGUGUGGUGGAAUUCUGCAG"
        "AUAUCCAGCACAGUGGCGGCCGCUCGAGUCUAGAGGGCCCGUUUAAACCCGCUGAUCAGCCU")   # PLACEHOLDER
spacer5, spacer3 = "AAACAAACAAA", "AAAGAAACAAAUAAA"
rna = ires + spacer5 + cds + spacer3

params = {
    "name": "circ_demo",
    "protein": PROTEIN,
    "rna": rna,
    "protected": [
        {"name": "IRES_PLACEHOLDER", "start": 0, "end": len(ires)},
        {"name": "spacer5", "start": len(ires), "end": len(ires) + len(spacer5)},
        {"name": "spacer3", "start": len(rna) - len(spacer3), "end": len(rna)},
    ],
    "primary_species": "mouse",
    "secondary_species": "human",
    "primary_weight": 0.6,
    "cai_floor": 0.78,
    "structure_goal": "minimize_pairing",
    "dsrna_min_bp": 30,
    "gc_min": 0.40, "gc_max": 0.58,
    "gc_window": 60, "gc_window_min": 0.30, "gc_window_max": 0.65,
    "max_poly_u": 4,
    "avoid_slippery": True, "avoid_g4": True, "avoid_cryptic_splice": True, "avoid_re_sites": True,
    "min_codon_freq": 0.08,
    "population": 48, "generations_per_round": 10, "rounds": 4,
    "elite_k": 6, "top_n": 6, "polish_codons": 6,
    "inner_span": 200, "stress_test": True, "linear_compare": True,
    "use_derna": False,
    "derna_lambdas": [0.0, 0.3],
    "seed": 1, "workers": 0,
}
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "example_params.json")
with open(out, "w") as fh:
    json.dump(params, fh, indent=2)
print(f"wrote {out}: {len(rna)} nt molecule, CDS {len(cds)} nt, protein {len(PROTEIN)} aa")
