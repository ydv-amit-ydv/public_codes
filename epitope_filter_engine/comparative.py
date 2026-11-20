"""Comparative prediction: score a query peptide against peptides with a known wet-lab outcome.

Why this exists. None of the predictors reports a probability that a peptide is immunogenic.
A NetMHCpan %Rank of 0.4 means "in the top 0.4 % of *random* peptides for binding": it is
calibrated against noise, not against real epitopes, and not against immunogenicity at all.
So an absolute cut-off answers a question nobody asked.

What this does instead. You give it peptides whose wet-lab outcome is known. They go through
the same pipeline as the query, and the query's score is then read against two empirical
distributions:

    positives   peptides that were tested and WORKED
    negatives   peptides that were tested and did NOT work (or, failing that, decoys)

The answer is a **tail likelihood ratio**: how many times more often a validated epitope
reaches this score than a non-epitope does. That is a real statement about evidence, and it
is the same logic as a positive and a negative control on a plate.

Positives alone cannot do this. Known epitopes score high, but so does a lot of junk, so a
percentile among positives only ever flatters the query: without a comparison group there is
no specificity and the number is not a confidence. Hence the negatives, and hence the honest
label when they are decoys rather than tested peptides.

What it deliberately does not hide (see `caveats`): published epitopes were often found
*because* a predictor flagged them and are skewed towards well-studied alleles, so recovering
them is partly circular; "positive" spans assays that ask different questions, and binding is
not immunogenicity; and a small panel gives a noisy baseline.

Everything here is pure: it works on the tables a finished run already saved, so it never
re-runs a predictor and never changes a stored result.
"""
import math
import re

import pandas as pd

MIN_GROUP = 8                  # below this a distribution is too small to read a percentile from
AA = set("ACDEFGHIKLMNPQRSTVWY")

# how a score is read: True = lower is better (a %Rank), False = higher is better
LOWER_IS_BETTER = {
    "Priority_Score": False, "Agree_Fraction_Best": False, "Agree_Fraction": False,
    "Coverage_Binder": False, "Coverage_Candidate": False,
    "Bind_Consensus_Best": True, "Bind_Consensus_%Rank": True, "Bind_Consensus_PopW": True,
    "NMP_%Rank_EL": True, "NMP_%Rank_BA": True, "NMP_Aff_nM": True,
    "NMP_Score_EL": False, "NMP_Score_BA": False,
    "MixMHC_%Rank": True, "PRIME_%Rank": True, "MHCF_Presentation_%Rank": True,
    "MHCF_Affinity_%Rank": True, "MHCF_Processing": False,
    "CTLpan_%Rank": True, "NetChop_Cterm": False,
    "BigMHC_EL": False, "BigMHC_IM": False, "IEDB_Immunogenicity": False,
    "NII_%Rank_EL": True, "NII_%Rank_BA": True, "MIX2_%Rank": True,
    "CD4_Episcore": True, "ImmuScope_IM": False, "TLimmuno2": False,
}
# the axis a verdict is based on, best first: the engine's combined score, then agreement, then binding
PRIMARY = ("Priority_Score", "Agree_Fraction_Best", "Bind_Consensus_Best")


# ---------------------------------------------------------------- the panel file

_CLASS_NAMES = {"1": 1, "i": 1, "mhc1": 1, "mhci": 1, "mhc-i": 1, "class1": 1, "classi": 1, "mhc class i": 1,
                "2": 2, "ii": 2, "mhc2": 2, "mhcii": 2, "mhc-ii": 2, "class2": 2, "classii": 2, "mhc class ii": 2}
_POSITIVE = re.compile(r"^(pos|positive|yes|y|true|1|immunogenic|epitope)", re.I)
_NEGATIVE = re.compile(r"^(neg|negative|no|n|false|0|non-?immunogenic|non-?epitope)", re.I)

# column name -> field. Our own layout first; the usual IEDB export headers are accepted too,
# so an export can be pasted without being rewritten by hand.
_HEADERS = {
    "class": "cls", "mhc": "cls", "mhc class": "cls", "mhc_class": "cls", "mhcclass": "cls",
    "mhc - class": "cls", "mhc allele class": "cls", "class i/ii": "cls",
    "epitope": "peptide", "peptide": "peptide", "sequence": "peptide", "description": "peptide",
    "epitope - name": "peptide", "epitope name": "peptide", "linear sequence": "peptide",
    "alleles": "alleles", "allele": "allele_one", "allele name": "allele_one",
    "mhc - allele names": "alleles", "mhc allele names": "alleles", "mhc restriction": "alleles",
    "hla": "alleles", "restriction": "alleles",
    "outcome": "outcome", "qualitative measure": "outcome", "measure": "outcome",
    "assay - qualitative measurement": "outcome", "qualitative measurement": "outcome",
    "result": "outcome", "immunogenicity": "outcome",
    "assay": "assay", "assay - method": "assay", "method": "assay", "method/technique": "assay",
    "assay method": "assay",
    "name": "name", "id": "name", "epitope id": "name", "reference": "name", "source": "name",
}


def _norm_header(h):
    return re.sub(r"\s+", " ", str(h).strip().lower().strip('"')).rstrip(":")


def _split_row(line):
    """One row of a pasted table: tab, comma or semicolon separated, quotes honoured."""
    if "\t" in line:
        parts = line.split("\t")
    elif ";" in line and line.count(";") >= line.count(","):
        parts = line.split(";")
    else:
        parts = next(__import__("csv").reader([line])) if "," in line else [line]
    return [p.strip().strip('"').strip() for p in parts]


def _looks_peptide(s):
    s = str(s).strip().upper()
    return len(s) >= 7 and set(s) <= AA


def _parse_class(v):
    return _CLASS_NAMES.get(_norm_header(v).replace(" ", "") if v else "", _CLASS_NAMES.get(_norm_header(v), None))


_PLACEHOLDER = {"na", "n/a", "-", "", "none", "unknown", "not determined", "nd"}


def _parse_alleles(v):
    """Split an allele cell. Not on '/': a class II allele is one name with a slash in it
    (HLA-DQA1*05:01/DQB1*02:01), so splitting there would break every DQ and DP pair."""
    if not v:
        return []
    return [a for a in (x.strip() for x in re.split(r"[,;|]| and ", str(v)))
            if a and a.lower() not in _PLACEHOLDER]


def parse_panel(text, default_class=None):
    """Panel text -> (DataFrame[Class, Peptide, Alleles, Outcome, Assay, Name], warnings).

    Our layout, one peptide per line, with or without a header row:

        MHC class , epitope , alleles (comma separated, optional) [, outcome] [, assay] [, name]
        1, GILGFVFTL, HLA-A*02:01
        2, PKYVKQNTLKLATGM, HLA-DRB1*01:01, positive
        1, LLFGYPVYV,, negative

    Alleles may be left empty; the comparison is then allele-agnostic for that peptide.
    Outcome may be left out entirely, in which case every row counts as a tested positive.
    A header row in IEDB's wording is also understood, so an export can be pasted as it is.
    """
    warns, rows = [], []
    lines = [l for l in (text or "").splitlines() if l.strip() and not l.lstrip().startswith("#")]
    if not lines:
        return pd.DataFrame(columns=["Class", "Peptide", "Alleles", "Outcome", "Assay", "Name"]), ["The panel is empty."]

    # a header row is one whose cells are known names and which holds no peptide-looking cell
    first = _split_row(lines[0])
    mapped = [_HEADERS.get(_norm_header(c)) for c in first]
    has_header = sum(m is not None for m in mapped) >= 2 and not any(
        len(c) >= 8 and set(c.upper()) <= AA for c in first)
    if has_header:
        cols, body = mapped, lines[1:]
        unknown = [c for c, m in zip(first, mapped) if m is None]
        if unknown:
            warns.append("Columns ignored: " + ", ".join(unknown[:8]))
        if "peptide" not in cols:
            return (pd.DataFrame(columns=["Class", "Peptide", "Alleles", "Outcome", "Assay", "Name"]),
                    ["No peptide column found. Expected a header such as 'epitope' or 'peptide', "
                     "or no header at all and the order: MHC class, epitope, alleles."])
    else:                                       # positional: class, peptide, alleles, [outcome], [assay], [name]
        cols, body = ["cls", "peptide", "alleles", "outcome", "assay", "name"], lines

    for n, line in enumerate(body, 1):
        cells = _split_row(line)
        # Without a header the first cell is the MHC class, but a bare peptide list is a
        # perfectly reasonable thing to paste. If the first cell is a peptide rather than a
        # class, shift the row along and let the class be worked out from the length below.
        if not has_header and cells and _parse_class(cells[0]) is None and _looks_peptide(cells[0]):
            cells = [""] + cells
        rec = {}
        for i, cell in enumerate(cells):
            if i < len(cols) and cols[i]:
                key = cols[i]
                if key == "allele_one":
                    rec.setdefault("alleles", "")
                    rec["alleles"] = (rec["alleles"] + "," + cell).strip(",") if cell else rec["alleles"]
                else:
                    rec[key] = cell
        pep = (rec.get("peptide") or "").upper().strip()
        if not pep:
            warns.append(f"row {n}: no peptide; skipped")
            continue
        bad = set(pep) - AA
        if bad:
            warns.append(f"row {n}: '{pep[:24]}' has non-standard residues {''.join(sorted(bad))}; skipped")
            continue
        cls = _parse_class(rec.get("cls")) or default_class
        if cls is None:                          # fall back to length, and say so
            cls = 1 if len(pep) <= 11 else 2
            warns.append(f"row {n}: '{pep[:18]}' has no MHC class; guessed class "
                         f"{'I' if cls == 1 else 'II'} from its length ({len(pep)} aa)")
        out = rec.get("outcome") or ""
        if not out:
            outcome = "positive"
        elif _POSITIVE.match(out):
            outcome = "positive"
        elif _NEGATIVE.match(out):
            outcome = "negative"
        else:
            outcome = "positive"
            warns.append(f"row {n}: outcome '{out[:20]}' not understood; counted as positive")
        rows.append({"Class": int(cls), "Peptide": pep, "Alleles": _parse_alleles(rec.get("alleles")),
                     "Outcome": outcome, "Assay": rec.get("assay") or "", "Name": rec.get("name") or ""})

    df = pd.DataFrame(rows, columns=["Class", "Peptide", "Alleles", "Outcome", "Assay", "Name"])
    if len(df):
        dup = df.duplicated(subset=["Peptide", "Class"], keep="first")
        if dup.any():
            warns.append(f"{int(dup.sum())} duplicate peptide(s) in the panel; the first row of each is used")
            df = df[~dup].reset_index(drop=True)
    return df, warns


# ---------------------------------------------------------------- riding along in a normal run

PANEL_PROTEIN = "__panel__"        # reserved 'Protein' marking a row that came from the panel, not the scan


def panel_rows(panel, cls, bounds):
    """Panel peptides of one class as pipeline rows [(Protein, Start, Peptide)].

    They are scored by the same models, in the same run, as the scan they will be compared
    with, which is what makes the two comparable. The reserved Protein name keeps them out of
    the scan's own ranking and out of the candidate list.
    """
    out, seen = [], set()
    for p in panel[panel["Class"] == cls]["Peptide"]:
        if bounds[0] <= len(p) <= bounds[1] and p not in seen:
            seen.add(p)
            out.append((PANEL_PROTEIN, len(out) + 1, p))
    return out


CORE = {1: 8, 2: 9}            # shared stretch that counts as "the same epitope region"

VALIDATED_MODES = ("flag", "prefer", "boost", "top")
NEGATIVE_MODES = ("flag", "demote", "remove")


def annotate(pep, panel, cls):
    """Mark each scan peptide against the wet-lab panel. Returns the table with:

    Validated           this exact peptide was tested and worked
    Validated_Negative  this exact peptide was tested and did NOT work
    Overlaps_Validated  shares a core-length stretch with one that worked, so it is probably
                        the same epitope seen at another offset or length
    Validated_Match     which panel peptide it matched
    Epitope_Status      validated / overlaps validated / tested negative / novel

    The overlap uses a shared 8-mer (class I) or 9-mer (class II), which is the MHC-bound
    core: two peptides sharing one are presenting the same thing.
    """
    pep = pep.copy()
    pan = panel[panel["Class"] == cls] if len(panel) else panel
    pos = list(pan[pan["Outcome"] == "positive"]["Peptide"]) if len(pan) else []
    neg = set(pan[pan["Outcome"] == "negative"]["Peptide"]) if len(pan) else set()
    k = CORE.get(cls, 9)
    cores = {}
    for p in pos:
        for i in range(max(1, len(p) - k + 1)):
            cores.setdefault(p[i:i + k], p)
    pos_set = set(pos)

    def mark(s):
        s = str(s)
        if s in pos_set:
            return True, False, s
        if s in neg:
            return False, True, s
        for i in range(max(1, len(s) - k + 1)):
            hit = cores.get(s[i:i + k])
            if hit:
                return False, False, hit
        return False, False, ""

    marks = [mark(s) for s in pep["Peptide"]] if len(pep) else []
    pep["Validated"] = [m[0] for m in marks]
    pep["Validated_Negative"] = [m[1] for m in marks]
    pep["Validated_Match"] = [m[2] for m in marks]
    pep["Overlaps_Validated"] = [bool(m[2]) and not m[0] and not m[1] for m in marks]
    pep["Epitope_Status"] = ["validated" if m[0] else "tested negative" if m[1]
                             else "overlaps validated" if m[2] else "novel" for m in marks]
    return pep


def score_flags(pep, rules):
    """Let the flags change the ranking, by the rule the user chose. Returns (table, note).

    The likelihood ratio is deliberately NOT used here. It is computed from the same
    predictor scores that already produced Priority_Score, so feeding it back would count
    the same evidence twice. Whether a peptide was validated in a wet lab is independent
    evidence, and that is what these rules act on.
    """
    mode = str(rules.get("validated_mode", "flag"))
    neg_mode = str(rules.get("validated_negative", "flag"))
    bonus = float(rules.get("validated_bonus", 0.10) or 0)
    if "Validated" not in pep.columns or not len(pep):
        return pep, ""
    pep = pep.copy()
    notes = []
    if "Priority_Score" in pep.columns:
        pep["Priority_Score_Raw"] = pep["Priority_Score"]
    val, negv = pep["Validated"].fillna(False), pep["Validated_Negative"].fillna(False)

    if neg_mode == "remove" and negv.any() and "Candidate" in pep.columns:
        pep.loc[negv, "Candidate"] = False
        notes.append(f"{int(negv.sum())} candidate(s) dropped: tested in the wet lab and did not work")
    if mode == "boost" and "Priority_Score" in pep.columns:
        pep.loc[val, "Priority_Score"] = (pep.loc[val, "Priority_Score"].fillna(0) + bonus).clip(upper=1.0)
        notes.append(f"validated epitopes given +{bonus:g} priority")
    if neg_mode == "demote" and "Priority_Score" in pep.columns:
        pep.loc[negv, "Priority_Score"] = (pep.loc[negv, "Priority_Score"].fillna(0) - bonus).clip(lower=0.0)
        notes.append(f"peptides that failed in the wet lab given -{bonus:g} priority")

    if mode in ("prefer", "boost", "top") or neg_mode in ("demote", "remove"):
        order = pd.DataFrame({
            "top": (~val).astype(int) if mode == "top" else 0,          # validated first
            "neg": negv.astype(int),                                     # failures last
            "pref": (~val).astype(int) if mode == "prefer" else 0,       # tie-break only
            "score": -pep["Priority_Score"].fillna(-1) if "Priority_Score" in pep.columns else 0})
        keys = ["top", "neg"] + (["score", "pref"] if mode == "prefer" else ["pref", "score"])
        pep = pep.assign(_o0=order["top"], _o1=order["neg"], _o2=order[keys[2]], _o3=order[keys[3]])
        pep = pep.sort_values(["_o0", "_o1", "_o2", "_o3"], kind="stable").drop(columns=["_o0", "_o1", "_o2", "_o3"])
        if "Priority_Rank" in pep.columns:
            pep["Priority_Rank"] = range(1, len(pep) + 1)
        pep = pep.reset_index(drop=True)
    return pep, "; ".join(notes)


def split_panel(df):
    """(scan rows, panel rows) of a table that carries a Protein column."""
    if "Protein" not in df.columns:
        return df, df.iloc[0:0]
    m = df["Protein"].astype(str) == PANEL_PROTEIN
    return df[~m], df[m]


# ---------------------------------------------------------------- the statistics

def _oriented(series, col):
    """Scores turned so that larger always means more epitope-like."""
    s = pd.to_numeric(series, errors="coerce")
    return -s if LOWER_IS_BETTER.get(col, False) else s


def tail_share(value, ref):
    """Share of `ref` at or above `value`, with a +0.5/+1 correction so it is never 0 or 1."""
    if not len(ref):
        return None
    return (sum(1 for x in ref if x >= value) + 0.5) / (len(ref) + 1)


def _density(value, sample, h):
    """Gaussian kernel density of `sample` at `value`."""
    if not sample or h <= 0:
        return 0.0
    c = 1.0 / (len(sample) * h * math.sqrt(2 * math.pi))
    return c * sum(math.exp(-0.5 * ((value - x) / h) ** 2) for x in sample)


def likelihood_ratio(value, pos, neg):
    """p(score | validated epitope) / p(score | non-epitope) at this score.

    A *density* ratio, not a tail ratio. The tail version ("how many positives reach this
    score or better, over how many negatives do") cannot argue against a peptide: when the
    positives sit at the top of the range, every low score still has nearly all positives
    above it, so the ratio saturates near 1 and a junk-looking peptide reads as merely
    inconclusive. Comparing the densities *at* the score fixes that: a score sitting in the
    negatives' bulk gives a ratio below 1, which is what "evidence against" has to mean.

    Both densities use one Gaussian kernel with a shared Silverman bandwidth over the pooled
    panel, and both are floored by the density a point three bandwidths away would have, so a
    score beyond every negative gives a large but finite ratio instead of dividing by zero.
    """
    if not pos or not neg or value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    pooled = list(pos) + list(neg)
    sd = _pstdev(pooled)
    if sd <= 0:                                  # every panel peptide scored the same
        return None
    h = 0.9 * sd * len(pooled) ** -0.2
    floor = math.exp(-4.5) / (len(pooled) * h * math.sqrt(2 * math.pi))
    # Evaluate inside the observed range. Far outside it both densities hit the floor and the
    # ratio collapses to 1, so a peptide scoring better than every validated epitope would read
    # as "inconclusive". Clamping judges it at the best score the panel actually shows, which is
    # the strongest claim the data can support without extrapolating off the end of it.
    v = min(max(value, min(pooled)), max(pooled))
    return (_density(v, pos, h) + floor) / (_density(v, neg, h) + floor)


def _pstdev(xs):
    n = len(xs)
    if n < 2:
        return 0.0
    m = sum(xs) / n
    return math.sqrt(sum((x - m) ** 2 for x in xs) / n)


def auc(pos, neg):
    """Probability that a random positive outscores a random negative (ties count a half).

    Mann-Whitney U on average ranks, so no scipy and ties are handled properly.
    """
    npos, nneg = len(pos), len(neg)
    if not npos or not nneg:
        return None
    allv = sorted([(v, 0) for v in pos] + [(v, 1) for v in neg])
    ranks, i = {}, 0
    while i < len(allv):
        j = i
        while j + 1 < len(allv) and allv[j + 1][0] == allv[i][0]:
            j += 1
        r = (i + j) / 2.0 + 1                    # average rank of the tied block
        for k in range(i, j + 1):
            ranks.setdefault(k, r)
        i = j + 1
    rsum = sum(ranks[k] for k, (_, grp) in enumerate(allv) if grp == 0)
    return (rsum - npos * (npos + 1) / 2.0) / (npos * nneg)


def percentile_of(value, ref):
    """Share of `ref` at or below `value`, as a percentage (None if there is nothing to compare)."""
    if not len(ref) or value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    return 100.0 * sum(1 for x in ref if x <= value) / len(ref)


def verdict(lr):
    """Plain words for a likelihood ratio. Bands follow the usual reading of diagnostic evidence."""
    if lr is None:
        return "no baseline", "There was nothing to compare this score with."
    if lr >= 10:
        return "strong support", "This score is typical of validated epitopes and rare among non-epitopes."
    if lr >= 3:
        return "moderate support", "This score is more typical of validated epitopes than of non-epitopes."
    if lr > 1 / 3:
        return "inconclusive", "Validated epitopes and non-epitopes reach this score about equally often."
    return "evidence against", "This score is more typical of non-epitopes than of validated epitopes."


def _groups(scored, panel, cls, length, col):
    """Oriented positives and negatives for one class, matched on length where there are enough.

    Score distributions differ by peptide length (a 9-mer and a 15-mer are not comparable), so
    the reference is taken from peptides of the same length; if either group would be too small
    the whole class is used instead and the caller says so.
    """
    p = panel[panel["Class"] == cls]
    pos_pep = set(p[p["Outcome"] == "positive"]["Peptide"])
    neg_pep = set(p[p["Outcome"] == "negative"]["Peptide"])
    tbl = scored[scored["Peptide"].isin(pos_pep | neg_pep)]
    same = tbl[tbl["Peptide"].str.len() == length] if length else tbl
    for frame, matched in ((same, True), (tbl, False)):
        pos = _oriented(frame[frame["Peptide"].isin(pos_pep)][col], col).dropna().tolist()
        neg = _oriented(frame[frame["Peptide"].isin(neg_pep)][col], col).dropna().tolist()
        if (len(pos) >= MIN_GROUP and len(neg) >= MIN_GROUP) or not matched:
            return pos, neg, matched
    return [], [], False


def caveats(panel, n_pos, n_neg, decoys, allele_agnostic):
    """What the numbers cannot tell you. Always reported with the result, never buried."""
    out = []
    if decoys:
        out.append("The negatives are length-matched decoys from the host proteome, not peptides that were tested "
                   "and failed. They are a weaker comparison: a decoy is only probably not an epitope. Add rows with "
                   "outcome 'negative' for a real comparison.")
    if n_pos < MIN_GROUP or n_neg < MIN_GROUP:
        out.append(f"The baseline is small ({n_pos} positive, {n_neg} negative). Percentiles and ratios from groups "
                   f"this size are noisy; treat the verdict as a hint, not a measurement.")
    out.append("Published epitopes were often found because a predictor flagged them, and they are skewed towards "
               "well-studied alleles such as HLA-A*02:01. A pipeline that ranks them highly is therefore partly "
               "agreeing with the predictors that found them, not independently confirming them.")
    if len(panel) and (panel["Assay"] == "").all():
        out.append("No assay was given for the panel. 'Positive' covers tetramer staining, IFN-gamma ELISpot and "
                   "MHC binding, which ask different questions: a peptide that binds MHC is not necessarily "
                   "immunogenic. Mixing them widens the baseline.")
    elif len(panel):
        kinds = sorted({a for a in panel["Assay"] if a})
        if len(kinds) > 1:
            out.append("The panel mixes assays (" + ", ".join(kinds[:5]) + "). They ask different questions, so the "
                       "baseline is wider than any single assay would give.")
    if allele_agnostic:
        out.append("Some panel peptides carry no allele, so for those the comparison uses each peptide's best score "
                   "over the alleles of the run rather than a matched allele.")
    return out


def analyse(scored, panel, queries, columns=None):
    """Compare each query peptide with the panel.

    scored   DataFrame of peptide-level results for panel + query peptides: one row per
             peptide, with the score columns of the run (what a finished run saved).
    panel    the parsed panel (parse_panel).
    queries  peptides to judge.
    columns  score columns to use; by default the primary axis plus every known column present.

    Returns {queries: [...], panel: {...}, caveats: [...]}: for each query a verdict on the
    primary axis and the evidence per column; for the panel, how well each column separated
    positives from negatives on this very data (AUC), which is also the honest way to see
    which predictor is worth trusting here.
    """
    scored = scored.copy()
    scored["Peptide"] = scored["Peptide"].astype(str).str.upper()
    present = [c for c in (columns or LOWER_IS_BETTER) if c in scored.columns]
    primary = next((c for c in PRIMARY if c in present), present[0] if present else None)
    pos_all = set(panel[panel["Outcome"] == "positive"]["Peptide"])
    neg_all = set(panel[panel["Outcome"] == "negative"]["Peptide"])
    decoys = bool(scored.attrs.get("decoys"))
    agnostic = bool(len(panel)) and any(not a for a in panel["Alleles"])

    # how well each column separates the panel's own positives from its negatives
    panel_cls = sorted(set(panel["Class"])) or [1]
    sep = {}
    for col in present:
        per_cls = {}
        for cls in panel_cls:
            pos, neg, _ = _groups(scored, panel, cls, None, col)
            a = auc(pos, neg)
            if a is not None:
                per_cls[f"class_{cls}"] = {"auc": round(a, 3), "n_pos": len(pos), "n_neg": len(neg)}
        if per_cls:
            sep[col] = per_cls

    out = []
    for q in queries:
        q = str(q).upper().strip()
        row = scored[scored["Peptide"] == q]
        if row.empty:
            out.append({"peptide": q, "error": "not in the scored table: it was not part of the run"})
            continue
        cls = int(row.iloc[0]["Class"]) if "Class" in row.columns else (panel_cls[0] if panel_cls else 1)
        ev, best = {}, None
        for col in present:
            v = _oriented(row[col], col).dropna()
            if v.empty:
                continue
            v = float(v.iloc[0])
            pos, neg, matched = _groups(scored, panel, cls, len(q), col)
            lr = likelihood_ratio(v, pos, neg)
            share = tail_share(v, pos)
            ev[col] = {"score": None if pd.isna(row.iloc[0][col]) else float(row.iloc[0][col]),
                       "pct_of_positives": None if percentile_of(v, pos) is None else round(percentile_of(v, pos), 1),
                       "pct_of_negatives": None if percentile_of(v, neg) is None else round(percentile_of(v, neg), 1),
                       "positives_at_least_as_good": None if share is None else round(100 * share, 1),
                       "likelihood_ratio": None if lr is None else round(lr, 2),
                       "n_pos": len(pos), "n_neg": len(neg), "length_matched": matched}
            if col == primary:
                best = lr
        label, why = verdict(best)
        out.append({"peptide": q, "length": len(q), "class": cls, "primary": primary,
                    "likelihood_ratio": None if best is None else round(best, 2),
                    "verdict": label, "explanation": why, "evidence": ev})

    n_pos = len(pos_all & set(scored["Peptide"]))
    n_neg = len(neg_all & set(scored["Peptide"]))
    return {"queries": out, "primary": primary,
            "panel": {"n_positive": n_pos, "n_negative": n_neg, "decoys": decoys,
                      "separation": sep, "classes": panel_cls},
            "caveats": caveats(panel, n_pos, n_neg, decoys, agnostic)}
