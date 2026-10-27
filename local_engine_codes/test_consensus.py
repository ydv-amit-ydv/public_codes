"""Known-answer tests for consensus.py and common/population.py (no prediction tool needed)."""
import glob, os, sys
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import consensus as C
from common.population import PopulationModel

# a tiny population: two regions weighted 0.6 / 0.4
tab = pd.DataFrame({"Allele": ["A*02:01", "A*11:01", "B*07:02"], "Locus": ["A", "A", "B"],
                    "f_X": [0.20, 0.10, 0.30], "f_Y": [0.10, 0.40, 0.10]})
tab["PopWeightedFreq"] = 0.6 * tab.f_X + 0.4 * tab.f_Y
pop = PopulationModel(tab)
assert abs(pop.region_weights["X"] - 0.6) < 1e-9
# one allele: 1-(1-f)^2 per region, then weighted
exp = 0.6 * (1 - 0.8 ** 2) + 0.4 * (1 - 0.9 ** 2)
assert abs(pop.coverage(["HLA-A*02:01"]) - exp) < 1e-12, pop.coverage(["HLA-A*02:01"])
# two loci independent: 1 - (1-cA)(1-cB)
cA = 0.6 * (1 - .8 ** 2) + 0.4 * (1 - .9 ** 2); cB = 0.6 * (1 - .7 ** 2) + 0.4 * (1 - .9 ** 2)
assert abs(pop.coverage(["HLA-A*02:01", "HLA-B*07:02"]) - (0.6 * (1 - .64 * .49) + 0.4 * (1 - .81 * .81))) < 1e-12
# same locus adds frequencies
assert abs(pop.coverage(["HLA-A*02:01", "HLA-A*11:01"]) - (0.6 * (1 - .7 ** 2) + 0.4 * (1 - .5 ** 2))) < 1e-12


def pairs(rows):
    return pd.DataFrame(rows, columns=["Protein", "Start", "Peptide", "Allele", "NMP_%Rank_EL", "MixMHC_%Rank", "MHCF_Presentation_%Rank",
                                       "IEDB_Immunogenicity"])

# peptide AAAAAAAAA: A*02:01 ranks 0.1 / 1 / 5 -> 2 of 3 vote -> majority binder, median 1; A*11:01 all fail
rows = [("P", 1, "AAAAAAAAA", "HLA-A*02:01", 0.1, 1.0, 5.0, 0.5), ("P", 1, "AAAAAAAAA", "HLA-A*11:01", 9, 8, 7, 0.5)]
R = C.resolve_rules(1, "iedb"); R.update(bind_min_agree="majority", bind_rank_max=2.0)
res = C.analyse(pairs(rows), 1, R, pop)
pr = res["pairs"]
a = pr[pr.Allele == "HLA-A*02:01"].iloc[0]
assert (a.Bind_Votes, a.Bind_Tools, bool(a.Binder), a["Bind_Consensus_%Rank"]) == (2, 3, True, 1.0), a
assert not pr[pr.Allele == "HLA-A*11:01"].iloc[0].Binder
pe = res["peptides"].iloc[0]
assert abs(pe.Coverage_Binder - pop.coverage(["HLA-A*02:01"])) < 1e-9 and pe.N_Alleles_Binder == 1
assert pe.Best_Allele == "HLA-A*02:01"
# "all" rule: 2 of 3 fails; "any": passes
for rule, want in (("all", False), ("any", True), ("majority", True)):
    R2 = dict(R, bind_min_agree=rule)
    assert bool(C.analyse(pairs(rows), 1, R2, pop)["pairs"].iloc[0].Binder) is want, rule
# 1 of 3 votes is not a majority
r1 = pairs([("P", 1, "AAAAAAAAA", "HLA-A*02:01", 0.1, 5.0, 5.0, 0.5)])
assert not C.analyse(r1, 1, dict(R, bind_min_agree="majority"), pop)["pairs"].iloc[0].Binder
# unscored tool does not count against: only 2 tools scored, both vote -> binder under "all"
r2 = pairs([("P", 1, "AAAAAAAAA", "HLA-A*02:01", 0.1, 1.0, np.nan, 0.5)])
assert bool(C.analyse(r2, 1, dict(R, bind_min_agree="all"), pop)["pairs"].iloc[0].Binder)

# relative immunogenicity: 8 binder peptides, top 25 % -> exactly the 2 best vote
rows = [("P", i + 1, "ACDEFGHIK"[i % 9] + "L" * 8, "HLA-A*02:01", 0.5, 0.5, 0.5, float(i)) for i in range(8)]
rows = [(p, s, "ACDEFGHIK"[i % 9] + "LLLLLLL" + "ACDEFGHIK"[s % 9], al, *rest) for (p, s, pep, al, *rest), i in zip(rows, range(8))]
R3 = dict(C.resolve_rules(1, "balanced"), imm_top_q=0.25, imm_min_agree="any")
d = C.analyse(pairs(rows), 1, R3, pop)["pairs"]
assert d.Vote_IEDB_IMM.sum() == 2 and set(d[d.Vote_IEDB_IMM == 1].Start) == {7, 8}, d[["Start", "TopPct_IEDB_IMM", "Vote_IEDB_IMM"]]
assert d.Pair_Candidate.sum() == 2

# overlap collapse: windows 1-9, 2-10 (overlap 8/9) merge; 30-38 does not; different protein does not
cand = pd.DataFrame({"Protein": ["P", "P", "P", "Q"], "Start": [1, 2, 30, 1], "End": [9, 10, 38, 9]})
cl = C.collapse_overlaps(cand, 0.7)
assert [c["members"] for c in cl] == [[0, 1], [2], [3]], cl
# anchored on the representative: 1-9, 4-12 (6/9=0.67<0.7) is separate even though chained via 2-10
cand = pd.DataFrame({"Protein": ["P"] * 3, "Start": [1, 2, 4], "End": [9, 10, 12]})
assert [c["members"] for c in C.collapse_overlaps(cand, 0.7)] == [[0, 1], [2]]

# greedy panel: peptide1 covers A*02:01, peptide2 covers B*07:02 -> panel needs both; third (A*02:01 again) adds nothing
rows = [("P", 1, "AAAAAAAAA", "HLA-A*02:01", .1, .1, .1, 1.0), ("P", 1, "AAAAAAAAA", "HLA-B*07:02", 50, 50, 50, 1.0),
        ("P", 40, "CCCCCCCCC", "HLA-A*02:01", 50, 50, 50, 1.0), ("P", 40, "CCCCCCCCC", "HLA-B*07:02", .1, .1, .1, 1.0),
        ("P", 80, "DDDDDDDDD", "HLA-A*02:01", .2, .2, .2, 1.0), ("P", 80, "DDDDDDDDD", "HLA-B*07:02", 50, 50, 50, 1.0)]
res = C.analyse(pairs(rows), 1, dict(C.resolve_rules(1, "iedb"), panel_target=0.99), pop)
assert [x["peptide"] for x in res["panel"]][:2] and len(res["panel"]) == 2, res["panel"]
assert abs(res["panel"][-1]["cumulative_coverage"] - pop.coverage(["HLA-A*02:01", "HLA-B*07:02"])) < 1e-9

# every preset runs and states every rule field
for cls in (1, 2):
    for name, r in C.preset_rules(cls).items():
        assert set(C.RULE_FIELDS) <= set(r), (cls, name)
# bad input is a readable error
try:
    C.resolve_rules(1, "balanced", {"rank_by": "nonsense"}); raise SystemExit("no error")
except ValueError:
    pass

# real Indian table, when available: coverage is monotone and bounded
for path in glob.glob(os.path.join(HERE, "..", "..", "HEV", "Outputs", "india_hla_allele_weightage.csv")):
    P = PopulationModel.from_csv(path)
    top = P.alleles_for_coverage("A", 0.9)
    cov = [P.coverage(["HLA-" + a for a in top[:k]]) for k in range(1, len(top) + 1)]
    assert all(b >= a - 1e-12 for a, b in zip(cov, cov[1:])) and 0 < cov[-1] <= 1
    print("real population table ok:", round(cov[-1], 3))
print("consensus tests passed")
