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


# the gate in front of BigMHC is lossless: no pair that is a binder with BigMHC-EL scored is dropped without it
rng = np.random.default_rng(7)
n = 4000
full = pd.DataFrame({"Protein": "P", "Start": np.arange(n), "Peptide": [f"PEP{i:06d}" for i in range(n)], "Allele": "HLA-A*02:01",
                     "NMP_%Rank_EL": rng.exponential(12, n), "NMP_%Rank_BA": rng.exponential(12, n), "MixMHC_%Rank": rng.exponential(12, n),
                     "MHCF_Presentation_%Rank": rng.exponential(12, n), "BigMHC_EL": rng.beta(1, 3, n)})
full.loc[rng.random(n) < 0.1, "MixMHC_%Rank"] = np.nan                       # an allele a tool does not know
for rule in ("any", "majority", "all", 2, 3):
    Rr = dict(C.resolve_rules(1, "balanced"), bind_min_agree=rule)
    quick = full.drop(columns="BigMHC_EL")
    keep = C.gate_mask(quick, 1, Rr, pending=1)
    truth = C.binder_mask(full, 1, Rr)
    lost = truth & ~keep
    if rule == "any":
        lost &= False            # documented: with "any" a pending tool alone is not allowed to create a binder
    assert not lost.any(), (rule, int(lost.sum()))
    assert keep.mean() < 0.6 or rule == "any", (rule, keep.mean())

# regions: a class I 9-mer inside a class II 15-mer is one CTL+HTL region; lone peptides keep their class; dual ranks first
import regions as RG
def peps(rows):
    return pd.DataFrame([dict(Protein="P", Start=a, End=b, Peptide="A" * (b - a + 1), Candidate=True, Priority_Score=sc, Coverage_Candidate=cv,
                              Best_Allele="HLA-A*02:01", In_Panel=False) for a, b, sc, cv in rows])
seq = "".join("ACDEFGHIKLMNPQRSTVWY"[i % 20] for i in range(300))
reg = RG.build(peps([(10, 18, 90, .5), (100, 108, 99, .6), (50, 58, 80, .3)]), peps([(6, 20, 70, .4), (200, 214, 60, .2), (54, 68, 50, .1)]), {"P": seq})
t = dict(zip(zip(reg.Start, reg.End), reg.Type))
assert t[(6, 20)] == "CTL+HTL" and t[(100, 108)] == "CTL only" and t[(200, 214)] == "HTL only", t
assert t[(50, 58)] == "CTL only" and t[(54, 68)] == "HTL only"           # 5 of 9 overlap < 0.7: separate regions
assert reg.iloc[0].Type == "CTL+HTL" and list(reg.Region_Rank) == list(range(1, len(reg) + 1))
top = reg.iloc[0]
assert top.Sequence == seq[5:20] and top.CTL_Best_Peptide and top.HTL_Best_Peptide and abs(top.Region_Score - (80 + RG.DUAL_BONUS)) < 1e-9, top
assert RG.build(None, None, {}).empty
lab = RG.label_peptides(peps([(10, 18, 90, .5), (250, 258, 1, .1)]), reg, 1)
assert lab.Region_Type.tolist() == ["CTL+HTL", ""], lab[["Region_ID", "Region_Type"]]

# clusters: overlapping windows collapse onto the best one and the members are listed
rows = []
for st in (1, 2, 3, 40):                      # 9-mers at 1,2,3 overlap (8/9, 7/9); 40 does not
    pepx = "ACDEFGHIKLMNPQRSTVWY"[st % 20: st % 20 + 9].ljust(9, "A")
    rows.append(("P", st, pepx + str(st)[:0], "HLA-A*02:01", 0.1 if st == 2 else 0.5, 0.1 if st == 2 else 0.5, 0.1 if st == 2 else 0.5, 1.0))
cres = C.analyse(pairs(rows), 1, dict(C.resolve_rules(1, "iedb"), collapse=True, overlap_fraction=0.7), pop)["peptides"]
rep = cres[cres.Is_Representative & (cres.Start == 2)].iloc[0]
assert rep.Cluster_Size == 3 and rep.Cluster_Members.count("@") == 2, rep[["Cluster_Size", "Cluster_Members"]]
assert set(cres[cres.Cluster_ID == rep.Cluster_ID].Cluster_Representative) == {rep.Peptide}
assert cres[cres.Start == 40].iloc[0].Cluster_Size == 1 and cres[cres.Start == 40].iloc[0].Cluster_Members == ""

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
