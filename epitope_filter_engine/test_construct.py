"""Tests for construct.py: assembly, properties, junction screening and the VaxiJen parser."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import construct as cn                                                           # noqa: E402

ok = 0


def eq(got, want, what):
    global ok
    assert got == want, f"{what}: got {got!r}, want {want!r}"
    ok += 1


def close(got, want, what, tol=0.05):
    global ok
    assert abs(got - want) <= tol, f"{what}: got {got!r}, want {want!r}"
    ok += 1


# ---------------------------------------------------------------- physicochemistry
# Checked against ExPASy ProtParam for these sequences.

# Gly-Gly: 2x75.07 - 18.02 = 132.12
close(cn.mw("GG"), 132.12, "MW of GG", 0.01)
close(cn.mw("A"), 89.09, "MW of a single residue is the free amino acid", 0.01)
eq(cn.mw(""), 0.0, "MW of nothing")
# insulin A chain, ProtParam: 2383.8 Da
close(cn.mw("GIVEQCCTSICSLYQLENYCN"), 2383.8, "MW of the insulin A chain", 1.5)

# pI: a peptide of only acidic residues must be acidic, only basic must be basic
assert cn.pi("DDDD") < 4.5 and cn.pi("KKKK") > 10, (cn.pi("DDDD"), cn.pi("KKKK"))
ok += 1
close(cn.pi("GGGG"), 5.55, "pI with no ionisable side chain sits between the termini", 0.6)
# adding a lysine must raise the pI, adding a glutamate must lower it
assert cn.pi("GGGGK") > cn.pi("GGGG") > cn.pi("GGGGE"), (cn.pi("GGGGK"), cn.pi("GGGG"), cn.pi("GGGGE"))
ok += 1

# GRAVY: all-Ile is the top of the Kyte-Doolittle scale, all-Arg the bottom
close(cn.gravy("IIII"), 4.5, "GRAVY of poly-Ile", 0.001)
close(cn.gravy("RRRR"), -4.5, "GRAVY of poly-Arg", 0.001)
# by hand: (-0.4+4.5+3.8-0.4+2.8+4.2+2.8-0.7+3.8)/9 = 20.4/9
close(cn.gravy("GILGFVFTL"), 20.4 / 9, "GRAVY of GILGFVFTL", 0.001)
eq(cn.gravy(""), 0.0, "GRAVY of nothing")

# aliphatic index: poly-Ala is 100 by definition, poly-Ile 390, glycine 0
close(cn.aliphatic_index("AAAA"), 100.0, "aliphatic index of poly-Ala", 0.01)
close(cn.aliphatic_index("IIII"), 390.0, "aliphatic index of poly-Ile", 0.01)
close(cn.aliphatic_index("GGGG"), 0.0, "aliphatic index with no aliphatic residues", 0.01)

# instability needs the table; a tiny one proves the formula: (10/N) * sum of weights
eq(cn.instability("AA", {"AA": 10.0}), 50.0, "instability of a 2-mer = (10/2)*10")
eq(cn.instability("A", {}), 0.0, "a single residue has no dipeptide")
eq(cn.instability("AAA", {}), round((10.0 / 3) * 2, 2), "unlisted pairs default to 1.0")

p = cn.properties("GILGFVFTLAAYSIINFEKL", {"AA": 1.0})
eq(sorted(p), ["aliphatic_index", "gravy", "instability", "length", "mw", "pi"], "property keys")
eq(p["length"], 20, "length")

# the checklist reads the numbers the right way round
v = dict((m.split(":")[0], okv) for okv, m in cn.verdicts(cn.properties("D" * 60, {})))
assert any("GRAVY" in k for k in v), v
good = cn.verdicts({"length": 400, "mw": 44000, "pi": 6.5, "gravy": -0.3, "aliphatic_index": 85, "instability": 30})
eq([o for o, _ in good], [True] * 5, "a well-behaved construct passes every check")
bad = cn.verdicts({"length": 10, "mw": 1000, "pi": 12.0, "gravy": 1.2, "aliphatic_index": 20, "instability": 60})
eq([o for o, _ in bad], [False] * 5, "a badly-behaved one fails every check")

# ---------------------------------------------------------------- assembly

# CTL ordering puts favoured proteasomal C-termini first and is stable otherwise
eq(cn.order_ctl(["SIINFEKQ", "GILGFVFTL", "AAAAAAAW", "KLGGALQAV"]),
   ["GILGFVFTL", "KLGGALQAV", "SIINFEKQ", "AAAAAAAW"], "CTL order")
eq(cn.order_ctl([]), [], "ordering nothing")

seq, parts = cn.assemble(ctl=["GILGFVFTL", "KLGGALQAK"], htl=["PKYVKQNTLKLATGM"],
                         bcell=["DDDDDDDD"], adjuvant="APPHALS", ctag="HHHHHH")
eq(seq,
   "APPHALS" "EAAAK" "GILGFVFTL" "AAY" "KLGGALQAK" "GGGS" "PKYVKQNTLKLATGM" "GGGS" "DDDDDDDD" "GGGS" "HHHHHH",
   "the assembled construct")
eq([x["kind"] for x in parts],
   ["adjuvant", "linker", "epitope", "linker", "epitope", "linker", "epitope", "linker", "epitope", "linker", "tag"],
   "part kinds in order")
# coordinates are contiguous, 1-based and inclusive
eq(parts[0]["start"], 1, "first part starts at 1")
eq(parts[-1]["end"], len(seq), "last part ends at the end")
assert all(parts[i]["end"] + 1 == parts[i + 1]["start"] for i in range(len(parts) - 1)), "parts must tile the sequence"
assert "".join(x["seq"] for x in parts) == seq, "parts must reconstruct the sequence"
ok += 2

# a single epitope needs no intra-block linker
s1, _ = cn.assemble(ctl=["GILGFVFTL"])
eq(s1, "GILGFVFTL", "one epitope, nothing added")
# an adjuvant alone is not given a dangling linker
s2, _ = cn.assemble(adjuvant="APPHALS")
eq(s2, "APPHALS", "adjuvant alone")
# only B-cell epitopes: KK between them, no block linkers
s3, _ = cn.assemble(bcell=["AAAA", "CCCC"])
eq(s3, "AAAAKKCCCC", "B-cell block")
# ordering can be switched off
s4, _ = cn.assemble(ctl=["SIINFEKQ", "GILGFVFTL"], optimise_ctl=False)
eq(s4, "SIINFEKQAAYGILGFVFTL", "CTL order kept as given")
# linkers can be overridden
s5, _ = cn.assemble(ctl=["AAAA", "CCCC"], linkers={"CTL_intra": "GGGGS"})
eq(s5, "AAAAGGGGSCCCC", "custom linker")
eq(cn.assemble()[0], "", "nothing in, nothing out")

# ---------------------------------------------------------------- junction screening

_, parts2 = cn.assemble(ctl=["AAAAAAAAA", "CCCCCCCCC"])      # AAAAAAAAA AAY CCCCCCCCC
j = cn.junctions(parts2, lengths=(9,))
# every 9-mer that crosses a seam, and none that lies inside one epitope
assert all(x["peptide"] not in ("AAAAAAAAA", "CCCCCCCCC") for x in j), j
assert any(x["peptide"] == "AAAAAAAAY" for x in j), [x["peptide"] for x in j]   # spans epitope 1 + linker
assert any("AAY" in x["peptide"] and x["peptide"].endswith("C") for x in j), [x["peptide"] for x in j]
ok += 3
# coordinates point back into the construct
full = "".join(p["seq"] for p in parts2)
assert all(full[x["start"] - 1:x["end"]] == x["peptide"] for x in j), "junction coordinates must be right"
ok += 1
# each junction says which parts it crosses
assert all(len(x["spans"]) >= 2 for x in j), j
ok += 1
# a construct with no seam has no junctions
eq(cn.junctions(cn.assemble(ctl=["GILGFVFTL"])[1], lengths=(9,)), [], "a single part has no junctions")
# windows longer than the construct produce nothing, and do not crash
eq(cn.junctions(parts2, lengths=(500,)), [], "window longer than the construct")
# two lengths give both sets
j2 = cn.junctions(parts2, lengths=(9, 15))
eq(sorted({x["length"] for x in j2}), [9, 15], "both window lengths")
s = cn.junction_summary(parts2, lengths=(9, 15))
eq(s["n"], len(j2), "summary counts match")

# ---------------------------------------------------------------- FASTA out

eq(cn.to_fasta(["GILGFVFTL", "SIINFEKL"]), ">seq_1\nGILGFVFTL\n>seq_2\nSIINFEKL\n", "FASTA from a list")
eq(cn.to_fasta([("epA", "GILGFVFTL")]), ">epA\nGILGFVFTL\n", "FASTA keeps given names")
eq(cn.to_fasta([]), "", "FASTA from nothing")
eq(cn.to_fasta(["A" * 130]).count("\n"), 4, "long sequences wrap at 60 columns")

# ---------------------------------------------------------------- VaxiJen parser

# the classic 2.0 wording
rows, w = cn.parse_vaxijen(">seq_1\nOverall Prediction for the Protective Antigen = 0.4521 ( Probable ANTIGEN ).\n"
                           ">seq_2\nOverall Prediction for the Protective Antigen = 0.2113 ( Probable NON-ANTIGEN ).")
eq([r["name"] for r in rows], ["seq_1", "seq_2"], "names")
eq([r["antigen"] for r in rows], [True, False], "antigen calls")
eq([r["score"] for r in rows], [0.4521, 0.2113], "scores")
eq(w, [], "clean parse warns about nothing")

# the wording the user's own parser was written for, with a percentage
rows, _ = cn.parse_vaxijen("seq_1  Probable IMMUNOGEN with a probability of 87%\n"
                           "seq_2  Probable NON-IMMUNOGEN with a probability of 12%")
eq([r["name"] for r in rows], ["seq_1", "seq_2"], "names from bare seq_N")
eq([r["antigen"] for r in rows], [True, False], "immunogen calls")
eq([r["score"] for r in rows], [0.87, 0.12], "percentages become 0-1 scores")

# a verdict with no score still parses, and says so
rows, w = cn.parse_vaxijen(">p1\nProbable ANTIGEN")
eq(rows[0]["antigen"], True, "verdict without a score")
eq(rows[0]["score"], None, "no score available")
assert any("no score" in x for x in w), w
ok += 1

# the threshold is applied when only a score is present
rows, _ = cn.parse_vaxijen(">p1\n= 0.5000 ( Probable ANTIGEN )", threshold=0.4)
eq(rows[0]["antigen"], True, "score above the threshold")
# nonsense input is reported, not silently empty
rows, w = cn.parse_vaxijen("the quick brown fox")
eq(rows, [], "nothing parsed")
assert w and "Nothing recognisable" in w[0], w
ok += 1
eq(cn.parse_vaxijen("")[0], [], "empty input")

# a round trip: build the FASTA the user would paste, then parse what comes back
fa = cn.to_fasta(["GILGFVFTL", "SIINFEKL"])
eq(fa.count(">"), 2, "two records to submit")
back, _ = cn.parse_vaxijen(">seq_1\n= 0.61 ( Probable ANTIGEN )\n>seq_2\n= 0.30 ( Probable NON-ANTIGEN )")
eq([r["name"] for r in back], ["seq_1", "seq_2"], "names line up with the submitted FASTA")

# ---------------------------------------------------------------- epitope + allele input

m, order, w = cn.parse_epitope_alleles("""
# a comment
GILGFVFTL  HLA-A*02:01, HLA-A*02:06
KLGGALQAK	HLA-A*03:01
SIINFEKL
GILGFVFTL  HLA-A*68:01
BADsequence!
""")
eq(order, ["GILGFVFTL", "KLGGALQAK", "SIINFEKL"], "order kept, duplicates merged not repeated")
eq(m["GILGFVFTL"], ["HLA-A*02:01", "HLA-A*02:06", "HLA-A*68:01"], "a repeated peptide merges its alleles")
eq(m["KLGGALQAK"], ["HLA-A*03:01"], "tab separated works")
eq(m["SIINFEKL"], [], "a peptide with no alleles is still kept")
eq("NOALLELES" in cn.parse_epitope_alleles("NOALLELES")[0], False, "O is not an amino acid")
assert any("non-standard" in x for x in w), w
ok += 1
eq(cn.parse_epitope_alleles("")[0], {}, "empty input")

# ---------------------------------------------------------------- coverage report
class _FakePop:
    """Two regions, one allele with data and one without, so the arithmetic is checkable."""
    label = "Test"
    def resolve(self, a):
        class R: pass
        r = R(); r.source = "none" if a.startswith("HLA-DP") else "exact"; return r
    def per_region(self, alleles):
        return {"North": 0.5 if alleles else 0.0, "South": 0.25 if alleles else 0.0,
                "national": 0.4 if alleles else 0.0}

rep = cn.coverage_report({"AAA": ["HLA-A*02:01"], "BBB": ["HLA-DPA1*01:03/DPB1*04:01"], "CCC": []}, _FakePop())
eq(rep["national"], 0.4, "national coverage")
eq([r["region"] for r in rep["regions"]], ["North", "South"], "regions listed, national not among them")
eq(rep["alleles_without_data"], ["HLA-DPA1*01:03/DPB1*04:01"], "untyped alleles reported, not counted")
eq(rep["n_epitopes_without_alleles"], 1, "epitopes with no alleles are counted and reported")
eq(rep["weakest_region"]["region"], "South", "weakest region identified")
eq(isinstance(rep["national"], float), True, "plain floats, so it serialises to json")
# an epitope whose only allele has no data contributes nothing rather than crashing
eq([e for e in rep["per_epitope"] if e["peptide"] == "BBB"][0]["national"], 0.0, "untyped allele gives no coverage")
eq(cn.coverage_report({"AAA": ["HLA-DPA1*01:03/DPB1*04:01"]}, _FakePop()).get("error") is not None, True,
   "nothing typed at all is an explained zero, not a silent one")
eq(cn.coverage_report({"AAA": ["HLA-A*02:01"]}, None).get("error") is not None, True, "no population file says so")

# the instability table is the published one, all 400 pairs
eq(len(cn.DIWV), 400, "DIWV is complete")

# ---------------------------------------------------------------- screening wrappers
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "Screening_scripts"))
import screen_runner as sr                                                       # noqa: E402

# "Non-Allergen" contains "allergen": a substring test would call everything positive, which
# is what the first version did
for call, token, want in (("Toxin", "toxin", True), ("Non-Toxin", "toxin", False),
                          ("Allergen", "allergen", True), ("Non-Allergen", "allergen", False),
                          ("IFN-γ inducer", "inducer", True), ("Non-inducer", "inducer", False),
                          ("Hemolytic", "hemolytic", True), ("Non-Hemolytic", "hemolytic", False),
                          ("", "toxin", False), (None, "toxin", False)):
    eq(sr.is_positive(call, token), want, f"is_positive({call!r})")

# every tool declares whether a positive call is bad news or good news
eq(sorted(k for k, v in sr.TOOLS.items() if v["adverse"]),
   ["algpred2", "hemopi2", "ntxpred2", "toxinpred3"], "adverse screens")
# ToxinPred2 is the older, protein-trained version of ToxinPred3 and must not come back:
# two versions of one model, and on epitopes it calls benign ones toxic
eq("toxinpred2" in sr.TOOLS, False, "ToxinPred2 is not shipped")
# the two heavy screens are optional, so a missing one is not a failure
eq(sorted(k for k, v in sr.TOOLS.items() if v.get("optional")), ["hemopi2", "ntxpred2"], "optional screens")
eq(sorted(k for k, v in sr.TOOLS.items() if not v.get("optional")),
   ["algpred2", "ifnepitope2", "toxinpred3"], "screens installed by default")
# allergen homology is a second, independent allergenicity method, not another model
assert hasattr(sr, "allergen_homology"), "homology check missing"
_h = sr.allergen_homology([], {"python": ""})
eq(list(_h.columns), ["Peptide", "Allergen_Homology", "Allergen_Identity", "Allergen_Match",
                      "Allergen_Align_Len", "Allergen_Evalue"], "homology columns")
try:
    sr.allergen_homology(["GILGFVFTL"], {"python": "/nonexistent/python"})
    raise AssertionError("homology must refuse without AlgPred2's blastp")
except sr.ScreenError:
    ok += 1
eq([k for k, v in sr.TOOLS.items() if not v["adverse"]], ["ifnepitope2"], "IFN-gamma is a benefit, not a flag")

import pandas as _pd                                                             # noqa: E402
_t = _pd.DataFrame({"Peptide": ["A", "B"], "Toxicity": ["Toxin", "Non-Toxin"],
                    "Allergenicity": ["Non-Allergen", "Non-Allergen"],
                    "IFN_gamma": ["IFN-γ inducer", "Non-inducer"]})
eq(sr.flags(_t), {"A": ["Toxicity"]}, "only adverse calls are flagged")
eq(sr.benefits(_t), {"A": ["IFN_gamma"]}, "IFN-gamma reported separately as a benefit")
eq(sr.flags(_pd.DataFrame({"Peptide": []})), {}, "an empty table flags nothing")

# a screen that is not installed is a note, never an exception that sinks the others
_tbl, _notes = sr.screen(["GILGFVFTL"], {"python": "/nonexistent/python"}, tools=["toxinpred3"])
eq(list(_tbl.columns), ["Peptide"], "a missing tool contributes no columns")
assert _notes and "not installed" in _notes[0], _notes
ok += 1
eq(sr.screen([], {"python": ""})[0].empty, True, "no peptides, no rows")
# ToxinPred3/2 and AlgPred2 reject a batch of one ("axis 1 is out of bounds", "Reshape your
# data... single sample"), so a filler is sent alongside and dropped. It must never surface,
# including when the peptide being screened IS the filler.
eq(sr.FILLER != sr.FILLER2, True, "two fillers, so a peptide equal to one still gets padding")
eq(len(sr.FILLER) >= 8 and set(sr.FILLER) <= set("ACDEFGHIKLMNPQRSTVWY"), True, "the filler is a real peptide")
try:
    sr.predict("nosuchtool", ["AAA"], {"python": "x"})
    raise AssertionError("an unknown screen must be refused")
except sr.ScreenError:
    ok += 1

print(f"construct + screening tests passed ({ok} checks)")
