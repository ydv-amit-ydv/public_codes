"""Tests for comparative.py: the panel file, and the statistics behind a confidence.

The statistics are checked against values worked out by hand, because a plausible-looking
but wrong likelihood ratio would be invisible in the page and would mislead a wet-lab
decision.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comparative as cp                                                         # noqa: E402

ok = 0


def eq(got, want, what):
    global ok
    assert got == want, f"{what}: got {got!r}, want {want!r}"
    ok += 1


def close(got, want, what, tol=1e-6):
    global ok
    assert got is not None and abs(got - want) <= tol, f"{what}: got {got!r}, want {want!r}"
    ok += 1


# ---------------------------------------------------------------- the panel file

# our own layout, no header: MHC class, epitope, alleles
df, w = cp.parse_panel("1, GILGFVFTL, HLA-A*02:01\n2, PKYVKQNTLKLATGM, HLA-DRB1*01:01\n1, LLFGYPVYV,")
eq(len(df), 3, "three rows")
eq(list(df["Class"]), [1, 2, 1], "classes")
eq(list(df["Peptide"]), ["GILGFVFTL", "PKYVKQNTLKLATGM", "LLFGYPVYV"], "peptides")
eq(df.iloc[0]["Alleles"], ["HLA-A*02:01"], "allele parsed")
eq(df.iloc[2]["Alleles"], [], "empty allele cell is allowed")
eq(set(df["Outcome"]), {"positive"}, "no outcome column means every row is a tested positive")
eq(w, [], "a clean panel warns about nothing")

# several alleles in one cell, any of the usual separators
df, _ = cp.parse_panel("1, GILGFVFTL, HLA-A*02:01; HLA-A*02:06 | HLA-B*07:02")
eq(df.iloc[0]["Alleles"], ["HLA-A*02:01", "HLA-A*02:06", "HLA-B*07:02"], "multi-allele cell")
# a class II allele is ONE name containing a slash: splitting on "/" would break every DQ/DP pair
df, _ = cp.parse_panel("2, PKYVKQNTLKLATGM, HLA-DQA1*05:01/DQB1*02:01, positive")
eq(df.iloc[0]["Alleles"], ["HLA-DQA1*05:01/DQB1*02:01"], "DQ pair kept whole")
df, _ = cp.parse_panel("2, PKYVKQNTLKLATGM, HLA-DQA1*05:01/DQB1*02:01; HLA-DPA1*01:03/DPB1*04:01")
eq(len(df.iloc[0]["Alleles"]), 2, "two class II pairs, not four chains")
eq(cp._parse_alleles("HLA-A*02:01, n/a, -"), ["HLA-A*02:01"], "placeholders dropped")

# class written in words
for text, want in (("I, GILGFVFTL", 1), ("II, PKYVKQNTLKLATGM", 2), ("MHC-I, GILGFVFTL", 1),
                   ("class II, PKYVKQNTLKLATGM", 2), ("mhc class i, GILGFVFTL", 1)):
    d, _ = cp.parse_panel(text)
    eq(int(d.iloc[0]["Class"]), want, f"class from {text.split(',')[0]!r}")

# the optional 4th column: this is what makes a real comparison possible
df, _ = cp.parse_panel("1, GILGFVFTL, HLA-A*02:01, positive\n1, AAAAAAAAA,, negative\n1, CCCCCCCCC,, Negative")
eq(list(df["Outcome"]), ["positive", "negative", "negative"], "outcome column")
for word in ("pos", "Positive-High", "yes", "1", "immunogenic"):
    d, _ = cp.parse_panel(f"1, GILGFVFTL,, {word}")
    eq(d.iloc[0]["Outcome"], "positive", f"{word!r} is positive")
for word in ("neg", "Negative", "no", "0", "non-immunogenic"):
    d, _ = cp.parse_panel(f"1, GILGFVFTL,, {word}")
    eq(d.iloc[0]["Outcome"], "negative", f"{word!r} is negative")
d, w = cp.parse_panel("1, GILGFVFTL,, maybe")
eq(d.iloc[0]["Outcome"], "positive", "an unreadable outcome is counted as positive")
assert any("not understood" in x for x in w), "and it says so"
ok += 1

# an IEDB-style header row is understood, so an export can be pasted unchanged
df, w = cp.parse_panel('"Description","MHC - Allele Names","Qualitative Measure","Assay - Method"\n'
                       '"GILGFVFTL","HLA-A*02:01","Positive-High","ELISPOT"\n'
                       '"AAAAAAAAA","HLA-A*02:01","Negative","ELISPOT"')
eq(list(df["Peptide"]), ["GILGFVFTL", "AAAAAAAAA"], "IEDB header: peptides")
eq(list(df["Outcome"]), ["positive", "negative"], "IEDB header: outcomes")
eq(df.iloc[0]["Assay"], "ELISPOT", "IEDB header: assay")
eq(df.iloc[0]["Alleles"], ["HLA-A*02:01"], "IEDB header: allele")

# our own header wording, tab separated
df, _ = cp.parse_panel("MHC class\tepitope\talleles\toutcome\n1\tGILGFVFTL\tHLA-A*02:01\tpositive")
eq(len(df), 1, "tab separated with our header")
eq(df.iloc[0]["Peptide"], "GILGFVFTL", "tab separated peptide")

# a RECOGNISED header with no peptide column is refused with an explanation, not guessed at
df, w = cp.parse_panel("class,alleles,outcome\n1,HLA-A*02:01,positive")
eq(len(df), 0, "recognised header without a peptide column -> nothing parsed")
assert w and "No peptide column" in w[0], w
ok += 1
# an unrecognised header is not a header: it falls through to the positional layout and the
# junk row is skipped by name, so nothing is parsed wrongly either way
df, w = cp.parse_panel("foo,bar\n1,GILGFVFTL")
eq(list(df["Peptide"]), ["GILGFVFTL"], "unrecognised header falls through to positions")
assert any("BAR" in x for x in w), w
ok += 1

# bad rows are skipped by name, never silently
df, w = cp.parse_panel("1, GILGFVFTL\n1, GILG1VFTL\n1,\n2, PKYVKQNTLKLATGM")
eq(list(df["Peptide"]), ["GILGFVFTL", "PKYVKQNTLKLATGM"], "bad rows dropped")
eq(len(w), 2, "one warning per bad row")
assert any("non-standard" in x for x in w) and any("no peptide" in x for x in w), w
ok += 1

# a missing class is guessed from length, and the guess is declared
df, w = cp.parse_panel("GILGFVFTL\nPKYVKQNTLKLATGM")
eq(list(df["Class"]), [1, 2], "class guessed from length")
eq(len(w), 2, "both guesses declared")

# duplicates collapse
df, w = cp.parse_panel("1, GILGFVFTL, HLA-A*02:01\n1, GILGFVFTL, HLA-A*02:06")
eq(len(df), 1, "duplicate peptide collapsed")
assert any("duplicate" in x for x in w), w
ok += 1

# comments and blank lines
df, _ = cp.parse_panel("# my panel\n\n1, GILGFVFTL\n\n")
eq(len(df), 1, "comments and blanks ignored")
eq(len(cp.parse_panel("")[0]), 0, "empty input")

# ---------------------------------------------------------------- the statistics

# AUC: worked out by hand. positives all above negatives -> 1.0, reversed -> 0.0
close(cp.auc([3, 4, 5], [0, 1, 2]), 1.0, "auc perfect")
close(cp.auc([0, 1, 2], [3, 4, 5]), 0.0, "auc inverted")
close(cp.auc([1, 2], [1, 2]), 0.5, "auc identical distributions = chance")
# one positive, one negative, tied: a tie counts a half
close(cp.auc([1], [1]), 0.5, "auc single tie")
# pos 2,4 vs neg 1,3: pairs (2>1),(2<3),(4>1),(4>3) -> 3/4
close(cp.auc([2, 4], [1, 3]), 0.75, "auc worked by hand")
eq(cp.auc([], [1]), None, "auc needs both groups")
eq(cp.auc([1], []), None, "auc needs both groups (other way)")

# percentile: share of the reference at or below the value
close(cp.percentile_of(5, [1, 2, 3, 4]), 100.0, "percentile above all")
close(cp.percentile_of(0, [1, 2, 3, 4]), 0.0, "percentile below all")
close(cp.percentile_of(2, [1, 2, 3, 4]), 50.0, "percentile in the middle")
eq(cp.percentile_of(1, []), None, "percentile needs a reference")

# the likelihood ratio is a DENSITY ratio, so it can argue both ways.
# A tail ratio cannot: with the positives at the top of the range, a low score still has
# nearly every positive above it, so the ratio saturates near 1 and junk reads as merely
# inconclusive. These cases pin that down.
pos, neg = [10, 11, 12, 13, 14, 15, 16, 17, 18], [0, 1, 2, 3, 4, 5, 6, 7, 8]
assert cp.likelihood_ratio(14, pos, neg) > 10, "a score in the positives' bulk must support"
assert cp.likelihood_ratio(4, pos, neg) < 1 / 3, "a score in the negatives' bulk must argue AGAINST"
close(cp.likelihood_ratio(9, pos, neg), 1.0, "a score between the two groups is uninformative", tol=0.15)
ok += 2
# identical distributions carry no information whichever score you ask about
for v in (1, 5, 9):
    close(cp.likelihood_ratio(v, [1, 3, 5, 7, 9], [1, 3, 5, 7, 9]), 1.0, f"no information at {v}", tol=0.02)
# beyond the panel's range the ratio stays finite and keeps its direction, rather than
# collapsing to 1 (which would call a peptide better than every known epitope "inconclusive")
assert cp.likelihood_ratio(99, pos, neg) > 10, "better than every positive must still support"
assert cp.likelihood_ratio(-99, pos, neg) < 1 / 3, "worse than every negative must still argue against"
ok += 2
# the ratio rises monotonically with the score, which is what makes the verdict readable
seq = [cp.likelihood_ratio(v, pos, neg) for v in (0, 4, 9, 14, 18)]
assert all(a < b for a, b in zip(seq, seq[1:])), seq
ok += 1
eq(cp.likelihood_ratio(5, [], [1]), None, "LR needs positives")
eq(cp.likelihood_ratio(5, [1], []), None, "LR needs negatives")
eq(cp.likelihood_ratio(None, pos, neg), None, "LR needs a value")
eq(cp.likelihood_ratio(5, [5, 5, 5], [5, 5, 5]), None, "a panel with one score everywhere gives no baseline")

# the tail share is reported alongside, as a plain "how many validated epitopes do this well"
close(cp.tail_share(10, pos), 9.5 / 10, "tail share, all positives at or above")
close(cp.tail_share(99, pos), 0.5 / 10, "tail share, none")
eq(cp.tail_share(1, []), None, "tail share needs a reference")

# verdict bands
eq(cp.verdict(19.0)[0], "strong support", "band >= 10")
eq(cp.verdict(5.0)[0], "moderate support", "band 3-10")
eq(cp.verdict(1.0)[0], "inconclusive", "band around 1")
eq(cp.verdict(0.1)[0], "evidence against", "band <= 1/3")
eq(cp.verdict(None)[0], "no baseline", "no baseline")

# score direction: a %Rank is better when lower, so it must be flipped before comparing
s = pd.Series([0.1, 5.0])
eq(list(cp._oriented(s, "NMP_%Rank_EL")), [-0.1, -5.0], "a %Rank is inverted")
eq(list(cp._oriented(s, "BigMHC_EL")), [0.1, 5.0], "a score is not inverted")
# and that flip must change the verdict: rank 0.1 is good, rank 50 is bad
_pos = [-0.2, -0.3, -0.1, -0.4, -0.5, -0.15, -0.25, -0.35, -0.45]        # oriented good ranks
_neg = [-40.0, -50.0, -60.0, -45.0, -55.0, -42.0, -48.0, -52.0, -58.0]   # oriented bad ranks
assert cp.likelihood_ratio(-0.1, _pos, _neg) > 10 > cp.likelihood_ratio(-50.0, _pos, _neg), "direction must drive the verdict"
ok += 1

# ---------------------------------------------------------------- analyse end to end

# a panel of 9 positives and 9 negatives, plus two queries: one epitope-like, one not
panel_text = "\n".join(
    [f"1, {p},, positive" for p in ("GILGFVFTL", "NLVPMVATV", "GLCTLVAML", "CLGGLLTMV", "YLLEMLWRL",
                                    "KLGGALQAK", "ILRGSVAHK", "RMFPNAPYL", "VLFSSDFRI")] +
    [f"1, {p},, negative" for p in ("AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD", "EEEEEEEEE", "FFFFFFFFF",
                                    "GGGGGGGGG", "HHHHHHHHH", "IIIIIIIII", "KKKKKKKKK")])
panel, pw = cp.parse_panel(panel_text)
eq(len(panel), 18, "panel parsed")
eq(pw, [], "panel clean")

rows = []
for p in panel["Peptide"]:
    good = panel[panel["Peptide"] == p].iloc[0]["Outcome"] == "positive"
    rows.append({"Peptide": p, "Class": 1, "Priority_Score": 0.8 if good else 0.2,
                 "NMP_%Rank_EL": 0.3 if good else 40.0})
rows.append({"Peptide": "QQQQQQQQQ", "Class": 1, "Priority_Score": 0.82, "NMP_%Rank_EL": 0.25})   # looks like an epitope
rows.append({"Peptide": "WWWWWWWWW", "Class": 1, "Priority_Score": 0.18, "NMP_%Rank_EL": 45.0})   # looks like junk
scored = pd.DataFrame(rows)

res = cp.analyse(scored, panel, ["QQQQQQQQQ", "WWWWWWWWW", "ZZZZZZZZZ"])
eq(res["primary"], "Priority_Score", "primary axis is the combined score")
eq(res["panel"]["n_positive"], 9, "positives counted")
eq(res["panel"]["n_negative"], 9, "negatives counted")
close(res["panel"]["separation"]["Priority_Score"]["class_1"]["auc"], 1.0, "panel separates perfectly here")
close(res["panel"]["separation"]["NMP_%Rank_EL"]["class_1"]["auc"], 1.0, "and so does the rank, once oriented")

q = {x["peptide"]: x for x in res["queries"]}
eq(q["QQQQQQQQQ"]["verdict"], "strong support", "epitope-like query")
eq(q["WWWWWWWWW"]["verdict"], "evidence against", "junk-like query")
assert q["QQQQQQQQQ"]["likelihood_ratio"] > 10, q["QQQQQQQQQ"]
assert q["WWWWWWWWW"]["likelihood_ratio"] < 1, q["WWWWWWWWW"]
ok += 2
eq(q["ZZZZZZZZZ"].get("error") is not None, True, "a peptide that was not scored is reported, not invented")
# the evidence is per column, with both percentiles and the counts behind them
ev = q["QQQQQQQQQ"]["evidence"]["Priority_Score"]
eq(ev["n_pos"], 9, "evidence carries the group sizes")
close(ev["pct_of_positives"], 100.0, "query above every positive")
close(ev["pct_of_negatives"], 100.0, "query above every negative")
eq(ev["score"], 0.82, "and the raw score")

# positives only: no negatives means no confidence, and it must say so rather than invent one
pos_only, _ = cp.parse_panel("\n".join(f"1, {p},, positive" for p in panel[panel["Outcome"] == "positive"]["Peptide"]))
res2 = cp.analyse(scored, pos_only, ["QQQQQQQQQ"])
eq(res2["panel"]["n_negative"], 0, "no negatives")
eq(res2["queries"][0]["verdict"], "no baseline", "positives alone cannot give a confidence")
assert any("small" in c for c in res2["caveats"]), res2["caveats"]
ok += 1

# the caveats are always present and name the real problems
assert any("because a predictor flagged them" in c for c in res["caveats"]), "circularity must be stated"
assert any("assay" in c.lower() for c in res["caveats"]), "assay heterogeneity must be stated"
ok += 2
scored.attrs["decoys"] = True
assert any("decoys" in c for c in cp.analyse(scored, panel, ["QQQQQQQQQ"])["caveats"]), "decoy fallback must be labelled"
ok += 1
scored.attrs["decoys"] = False

# length matching: a 15-mer query is not compared with 9-mer references when enough 15-mers exist
mixed = pd.concat([scored, pd.DataFrame([
    {"Peptide": "P" * 15, "Class": 1, "Priority_Score": 0.5, "NMP_%Rank_EL": 1.0}])], ignore_index=True)
pos15 = [f"{'A' * i}{'Y' * (15 - i)}" for i in range(1, 10)]
neg15 = [f"{'C' * i}{'W' * (15 - i)}" for i in range(1, 10)]
big_panel, _ = cp.parse_panel("\n".join([f"1, {p},, positive" for p in pos15] + [f"1, {p},, negative" for p in neg15]))
rows15 = [{"Peptide": p, "Class": 1, "Priority_Score": 0.9, "NMP_%Rank_EL": 0.2} for p in pos15] + \
         [{"Peptide": p, "Class": 1, "Priority_Score": 0.1, "NMP_%Rank_EL": 50.0} for p in neg15] + \
         [{"Peptide": "P" * 15, "Class": 1, "Priority_Score": 0.95, "NMP_%Rank_EL": 0.15}]
r15 = cp.analyse(pd.DataFrame(rows15), big_panel, ["P" * 15])
eq(r15["queries"][0]["evidence"]["Priority_Score"]["length_matched"], True, "15-mers compared with 15-mers")
eq(r15["queries"][0]["verdict"], "strong support", "15-mer query judged on its own length")

print(f"comparative tests passed ({ok} checks)")
