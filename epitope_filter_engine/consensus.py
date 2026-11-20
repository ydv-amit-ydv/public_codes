"""Consensus across tools: votes, agreement, a population-weighted peptide ranking, a minimum peptide panel.

Input is the pair table the pipeline builds (one row per peptide x allele, one column per tool output); output is the
same table plus consensus columns, a one-row-per-peptide table, the shortlist, and a report (funnel, suggestions).
Everything is plain pandas on columns, so it is fast and testable without any prediction tool.

THE METHOD, so the numbers can be defended
------------------------------------------
Every tool answers one of three questions: is this a ligand of the allele (binding), is it immunogenic, is it
processed. Each tool casts one vote per peptide-allele pair:

  * Tools that report a PERCENTILE against random human peptides (NetMHCpan, MixMHCpred, MHCflurry, NetMHCIIpan,
    MixMHC2pred, PRIME, NetCTLpan) vote on that absolute scale: percentile <= a cut-off (binding 2 %, the field's
    usual weak-binder line, unless a preset says otherwise).
  * BigMHC EL has no percentile; it votes at a score cut-off (default 0.2, which sits at about the top 2 % of random
    human 9-mers for A*02:01 and B*07:02 when measured).
  * Tools with no absolute scale for immunogenicity (BigMHC IM, the IEDB/Calis score, CD4episcore, ImmuScope, TLimmuno2) vote RELATIVELY:
    top q % among the binder pairs of this run (as the Workbench does: immunogenicity matters only for peptides
    that are presented). PRIME also has a relative rank for the consensus mean, besides its absolute vote.

Agreement = votes / tools that scored the pair, per group and overall. A pair is a BINDER when enough binding tools
agree (strict majority by default). The binding consensus is the MEDIAN percentile of the percentile tools (IEDB's own
way of combining methods). Immunogenicity consensus is the mean relative rank of its tools.

Per peptide, alleles are weighted by how common they are in the chosen population (national allele frequency):
  * Coverage = share of the population carrying at least one allele the peptide is a binder for (Hardy-Weinberg,
    loci independent; see common/population.py)
  * population-weighted agreement and population-weighted binding percentile (geometric mean, so a few non-binding
    alleles at percentile 60 do not swamp the score)
The Priority score is the weighted mean of the peptide's percentile standing on coverage, agreement, binding and
immunogenicity (weights are rules, equal-ish by default; each component is in the CSV so any other ranking can be rebuilt).
"""
import math
from dataclasses import dataclass

import numpy as np
import pandas as pd

from common import hostself

from common import peptide_props

KEYS = ["Protein", "Start", "Peptide"]
MIN_RELATIVE = 8          # a "top q %" vote needs enough binder pairs to rank among


@dataclass(frozen=True)
class Tool:
    key: str
    label: str
    col: str
    group: str            # binding | immunogenicity | processing
    kind: str             # abs (percentile, lower better) | score (own score cut-off) | rel (relative rank)
    direction: str        # low | high  (which end is better)
    allele_dependent: bool
    cls: int


TOOLS = [
    Tool("NMP_EL", "NetMHCpan-EL", "NMP_%Rank_EL", "binding", "abs", "low", True, 1),
    Tool("NMP_BA", "NetMHCpan-BA", "NMP_%Rank_BA", "binding", "abs", "low", True, 1),
    Tool("MIXMHC", "MixMHCpred", "MixMHC_%Rank", "binding", "abs", "low", True, 1),
    Tool("MHCF", "MHCflurry", "MHCF_Presentation_%Rank", "binding", "abs", "low", True, 1),
    Tool("BIGMHC_EL", "BigMHC-EL", "BigMHC_EL", "binding", "score", "high", True, 1),
    Tool("PRIME", "PRIME", "PRIME_%Rank", "immunogenicity", "abs", "low", True, 1),
    Tool("BIGMHC_IM", "BigMHC-IM", "BigMHC_IM", "immunogenicity", "rel", "high", True, 1),
    Tool("IEDB_IMM", "IEDB-Calis", "IEDB_Immunogenicity", "immunogenicity", "rel", "high", True, 1),
    Tool("CTLPAN", "NetCTLpan", "CTLpan_%Rank", "processing", "abs", "low", True, 1),
    Tool("NETCHOP", "NetChop", "NetChop_Cterm", "processing", "score", "high", False, 1),
    Tool("NII_EL", "NetMHCIIpan-EL", "NII_%Rank_EL", "binding", "abs", "low", True, 2),
    Tool("NII_BA", "NetMHCIIpan-BA", "NII_%Rank_BA", "binding", "abs", "low", True, 2),
    Tool("MIX2", "MixMHC2pred", "MIX2_%Rank", "binding", "abs", "low", True, 2),
    Tool("CD4", "CD4episcore", "CD4_Episcore", "immunogenicity", "rel", "low", False, 2),
    Tool("IMMUSCOPE", "ImmuScope-IM", "ImmuScope_IM", "immunogenicity", "rel", "high", True, 2),
    Tool("TLIMMUNO2", "TLimmuno2", "TLimmuno2", "immunogenicity", "rel", "high", True, 2),
]
GROUPS = ("binding", "immunogenicity", "processing")
GROUP_PREFIX = {"binding": "Bind", "immunogenicity": "Imm", "processing": "Proc"}

# ---------------------------------------------------------------------------------------------- rules
RULE_FIELDS = ("bind_rank_max", "strong_rank_max", "bigmhc_el_min", "bind_min_agree", "imm_top_q", "imm_abs_rank_max",
               "imm_min_agree", "ctlpan_rank_max", "netchop_min", "require_processing", "min_alleles", "min_coverage",
               "collapse", "overlap_fraction", "cluster_by", "host_mode", "host_near_min", "panel_target", "panel_cap", "rank_by", "rank_weights")
RANK_BY = ("consensus", "coverage", "promiscuity", "binding", "immunogenicity")
DEFAULT_WEIGHTS = {"coverage": 0.35, "agreement": 0.25, "binding": 0.20, "immunogenicity": 0.20}


def _rules(cls, **over):
    base = {"bind_rank_max": 2.0 if cls == 1 else 5.0, "strong_rank_max": 0.5 if cls == 1 else 2.0, "bigmhc_el_min": 0.2,
            "bind_min_agree": "majority", "imm_top_q": 0.25, "imm_abs_rank_max": 5.0, "imm_min_agree": "any",
            "ctlpan_rank_max": 1.0, "netchop_min": None, "require_processing": False, "min_alleles": 1,
            "min_coverage": 0.0, "host_mode": "flag", "host_near_min": 8 / 9, "collapse": True, "overlap_fraction": 0.7, "cluster_by": "core" if cls == 2 else "peptide", "panel_target": 0.95, "panel_cap": 25,
            "rank_by": "consensus", "rank_weights": dict(DEFAULT_WEIGHTS)}
    base.update(over)
    return base


# Every preset states every field (nothing is inherited from whatever was set last), and says why.
PRESET_INFO = {
    "balanced": ("Balanced", "A sensible default: a peptide must be called a binder by a majority of the binding tools, at least one "
                             "immunogenicity tool must rank it in the top quarter of the binders, and overlapping windows collapse to one."),
    "strict": ("Strict", "A short list for synthesis: binding by every binding tool at a tight cut-off, immunogenicity by most tools "
                         "in the top tenth, and at least two alleles; peptides with a 9-mer identical to a host protein are removed."),
    "coverage": ("Max coverage", "Broadly presented peptides: a majority of binding tools, no immunogenicity requirement, at least four "
                                 "alleles, ranked by population coverage. For maximising how many people a panel reaches."),
    "iedb": ("IEDB cut-offs", "IEDB's own suggested starting point: percentile <= 1 for class I and <= 10 for class II, any single tool."),
    "wide": ("Wide net", "A first look: any one binding tool at a loose cut-off, no immunogenicity requirement, overlapping windows kept."),
}


def preset_rules(cls):
    if cls == 1:
        return {"balanced": _rules(1),
                "strict": _rules(1, bind_rank_max=0.5, strong_rank_max=0.5, bigmhc_el_min=0.5, bind_min_agree="all", imm_top_q=0.10,
                                 imm_abs_rank_max=1.0, imm_min_agree="majority", min_alleles=2, host_mode="exclude_identical"),
                "coverage": _rules(1, imm_min_agree="off", min_alleles=4, rank_by="coverage"),
                "iedb": _rules(1, bind_rank_max=1.0, bind_min_agree="any", imm_min_agree="off"),
                "wide": _rules(1, bind_rank_max=5.0, bigmhc_el_min=0.1, bind_min_agree="any", imm_min_agree="off", collapse=False)}
    return {"balanced": _rules(2),
            "strict": _rules(2, bind_rank_max=2.0, bind_min_agree="all", imm_top_q=0.10, imm_min_agree="majority", min_alleles=2, host_mode="exclude_identical"),
            "coverage": _rules(2, imm_min_agree="off", min_alleles=4, rank_by="coverage"),
            "iedb": _rules(2, bind_rank_max=10.0, strong_rank_max=2.0, bind_min_agree="any", imm_min_agree="off"),
            "wide": _rules(2, bind_rank_max=10.0, bind_min_agree="any", imm_min_agree="off", collapse=False)}


def resolve_rules(cls, preset="balanced", overrides=None):
    presets = preset_rules(cls)
    r = dict(presets.get(preset, presets["balanced"]))
    r["rank_weights"] = dict(r["rank_weights"])
    for k, v in (overrides or {}).items():
        if k not in RULE_FIELDS or v in (None, "") and k not in ("netchop_min",):
            continue
        r[k] = ({**r["rank_weights"], **v} if k == "rank_weights" else v)
    if r.get("host_mode") not in hostself.MODES:
        raise ValueError(f"host_mode must be one of {hostself.MODES}")
    if not 0 < float(r.get("host_near_min", 8 / 9)) <= 1:
        raise ValueError("host_near_min must be a fraction between 0 and 1")
    if r.get("cluster_by") not in ("peptide", "core"):
        raise ValueError("cluster_by must be peptide or core")
    if r["rank_by"] not in RANK_BY:
        raise ValueError(f"rank_by must be one of {RANK_BY}")
    for k in ("bind_min_agree", "imm_min_agree"):
        v = r[k]
        if not (v in ("any", "majority", "all", "off") or (isinstance(v, (int, float)) and v >= 1)):
            raise ValueError(f"{k} must be any, majority, all, off or a number of tools")
    return r


def _need(rule, n):
    """Votes required out of n tools. n is a Series of ints; returns floats (inf when it can never be met)."""
    n = n.astype(float)
    if rule == "off":
        return pd.Series(0.0, index=n.index)
    if rule == "any":
        need = pd.Series(1.0, index=n.index)
    elif rule == "majority":
        need = np.floor(n / 2) + 1
    elif rule == "all":
        need = n
    else:
        need = np.minimum(float(rule), n)
    return need.where(n > 0, np.inf)


def _vote(cond, series):
    v = cond.astype("Int8")
    v[series.isna()] = pd.NA
    return v


def active_tools(df, cls):
    return [t for t in TOOLS if t.cls == cls and t.col in df.columns and df[t.col].notna().any()]


def _binding_votes(d, tools, rules):
    """(votes, scored) per row over the binding tools that have a value right now."""
    votes = pd.Series(0, index=d.index)
    scored = pd.Series(0, index=d.index)
    for t in (t for t in tools if t.group == "binding"):
        c = d[t.col]
        yes = (c >= rules["bigmhc_el_min"]) if t.key == "BIGMHC_EL" else (c <= rules["bind_rank_max"])
        votes += yes.fillna(False).astype(int)
        scored += c.notna().astype(int)
    return votes, scored


def gate_mask(df, cls, rules, pending=0):
    """Pairs worth sending to a slow binding tool (`pending` of them still to run, e.g. BigMHC-EL = 1).

    Lossless for the agreement rule: a pair is kept when the quick tools' votes plus every pending vote could still reach
    the number of votes the rule needs. A pair that fails this can never become a binder, so skipping it changes nothing
    in the result. With the "any one tool" rule a pending tool alone could make a binder; it is a confirmation tool
    there, so at least one quick vote is still asked for.
    """
    tools = [t for t in active_tools(df, cls) if t.group == "binding" and t.kind == "abs"]
    if not tools:
        return pd.Series(True, index=df.index)
    votes, scored = _binding_votes(df, tools, rules)
    need = _need(rules["bind_min_agree"], scored + pending)
    return votes >= np.maximum(1, need - pending)


def binder_mask(df, cls, rules):
    """Binder call from the binding tools present now (used to decide which pairs the immunogenicity tools must score)."""
    tools = [t for t in active_tools(df, cls) if t.group == "binding"]
    if not tools:
        return pd.Series(True, index=df.index)
    votes, scored = _binding_votes(df, tools, rules)
    return (scored > 0) & (votes >= _need(rules["bind_min_agree"], scored))


# ---------------------------------------------------------------------------------------------- votes
def add_votes(d, cls, rules):
    """Adds Vote_*, TopPct_*, Plausible, group counts, Agree_*, consensus columns. Returns (d, tools)."""
    tools = active_tools(d, cls)
    for t in tools:
        c = d[t.col]
        if t.kind == "abs":
            thr = {"binding": rules["bind_rank_max"], "immunogenicity": rules["imm_abs_rank_max"], "processing": rules["ctlpan_rank_max"]}[t.group]
            d[f"Vote_{t.key}"] = _vote(c <= thr, c)
        elif t.key == "BIGMHC_EL":
            d[f"Vote_{t.key}"] = _vote(c >= rules["bigmhc_el_min"], c)
        elif t.key == "NETCHOP" and rules["netchop_min"] is not None:
            d[f"Vote_{t.key}"] = _vote(c >= float(rules["netchop_min"]), c)
    d["Plausible"] = binder_mask(d, cls, rules)         # immunogenicity is ranked among the pairs presented (binders)

    # relative ranks among the binders (lower = better, in percent)
    d.attrs["notes"] = []
    for t in (t for t in tools if t.group == "immunogenicity"):
        s = d[t.col].where(d["Plausible"])
        if t.kind == "rel" and s.notna().sum() < MIN_RELATIVE:
            d.attrs["notes"].append(f"{t.label} gives no vote: only {int(s.notna().sum())} binder pair(s), too few to rank among "
                                    f"(needs {MIN_RELATIVE}). Its raw values are still in the tables.")
            tools = [x for x in tools if x is not t]
            continue
        asc = t.direction == "low"
        if t.allele_dependent:
            pct = s.groupby(d["Allele"]).rank(pct=True, ascending=asc, method="average")
        else:
            u = d.loc[s.notna(), KEYS + [t.col]].drop_duplicates(KEYS)
            u["_p"] = u[t.col].rank(pct=True, ascending=asc, method="average")
            pct = d[KEYS].merge(u[KEYS + ["_p"]], on=KEYS, how="left")["_p"].set_axis(d.index)
        d[f"TopPct_{t.key}"] = (pct * 100).round(4)
        if t.kind == "rel":
            d[f"Vote_{t.key}"] = _vote(d[f"TopPct_{t.key}"] <= rules["imm_top_q"] * 100, d[f"TopPct_{t.key}"])
    top_cols = [f"TopPct_{t.key}" for t in tools if t.group == "immunogenicity"]
    d["Imm_TopPct_Consensus"] = d[top_cols].mean(axis=1, skipna=True).round(4) if top_cols else np.nan

    for g in GROUPS:
        cols = [f"Vote_{t.key}" for t in tools if t.group == g and f"Vote_{t.key}" in d]
        p = GROUP_PREFIX[g]
        d[f"{p}_Votes"] = (d[cols] == 1).sum(axis=1).astype(int) if cols else 0
        d[f"{p}_Tools"] = d[cols].notna().sum(axis=1).astype(int) if cols else 0
    d["Agree_Votes"] = d["Bind_Votes"] + d["Imm_Votes"] + d["Proc_Votes"]
    d["Agree_Tools"] = d["Bind_Tools"] + d["Imm_Tools"] + d["Proc_Tools"]
    d["Agree_Fraction"] = (d["Agree_Votes"] / d["Agree_Tools"].where(d["Agree_Tools"] > 0)).round(4)
    vcols = [(f"Vote_{t.key}", t.label) for t in tools if f"Vote_{t.key}" in d]
    mat = np.column_stack([(d[c] == 1).fillna(False).to_numpy(bool) for c, _ in vcols]) if vcols else np.zeros((len(d), 0), bool)
    labels = np.array([l for _, l in vcols], dtype=object)
    d["Agree_List"] = [";".join(labels[row]) for row in mat] if vcols else ""

    pc = [t.col for t in tools if t.group == "binding" and t.kind == "abs"]
    d["Bind_Consensus_%Rank"] = d[pc].median(axis=1, skipna=True) if pc else np.nan
    return d, tools


def add_calls(d, tools, rules, filters_ok=None):
    d["Binder"] = (d["Bind_Tools"] > 0) & (d["Bind_Votes"] >= _need(rules["bind_min_agree"], d["Bind_Tools"]))
    d["Strong"] = d["Binder"] & (d["Bind_Consensus_%Rank"] <= rules["strong_rank_max"])
    has_imm = d["Imm_Tools"].max() > 0
    need_i = _need(rules["imm_min_agree"], d["Imm_Tools"])
    d["Imm_OK"] = (d["Imm_Votes"] >= need_i) & (d["Imm_Tools"] > 0) if (has_imm and rules["imm_min_agree"] != "off") else True
    d["Proc_OK"] = ((d["Proc_Tools"] > 0) & (d["Proc_Votes"] >= 1)) if rules["require_processing"] and d["Proc_Tools"].max() > 0 else True
    d["Filters_OK"] = True if filters_ok is None else filters_ok
    d["Pair_Candidate"] = d["Binder"] & d["Imm_OK"] & d["Proc_OK"] & d["Filters_OK"]
    d["Pair_Call"] = np.where(d["Strong"], "strong binder", np.where(d["Binder"], "binder", ""))
    return d


# ---------------------------------------------------------------------------------------------- peptide level
def _wmean(values, weights):
    ok = values.notna() & (weights > 0)
    return float(np.average(values[ok], weights=weights[ok])) if ok.any() else np.nan


def peptide_table(d, tools, rules, pop, cls):
    """One row per peptide: aggregates over its alleles, population coverage, priority score."""
    n_alleles = d["Allele"].nunique()
    weights = {a: (pop.allele_info(a)[0] if pop else 0.0) for a in d["Allele"].unique()}
    d["Allele_Freq"] = d["Allele"].map(weights)
    d["Allele_Data"] = d["Allele"].map(lambda a: pop.allele_info(a)[1] if pop else "none")
    has_pop = bool(pop) and any(w > 0 for w in weights.values())

    g = d.groupby(KEYS, sort=False)
    out = g.agg(Length=("Peptide", lambda s: len(s.iloc[0])), Alleles_Tested=("Allele", "nunique"),
                N_Alleles_Binder=("Binder", "sum"), N_Alleles_Strong=("Strong", "sum"),
                N_Alleles_Candidate=("Pair_Candidate", "sum")).reset_index()
    out["End"] = out["Start"] + out["Length"] - 1

    def listed(mask, name):
        s = d[mask].groupby(KEYS, sort=False)["Allele"].agg(lambda a: ", ".join(sorted(a)))
        return s.rename(name).reset_index()

    out = out.merge(listed(d["Binder"], "Alleles_Binder"), on=KEYS, how="left")
    out = out.merge(listed(d["Pair_Candidate"], "Alleles_Candidate"), on=KEYS, how="left")
    for c in ("Alleles_Binder", "Alleles_Candidate"):
        out[c] = out[c].fillna("")

    # best allele = the binder allele with the lowest binding consensus percentile (or highest Agree_Fraction without one)
    b = d[d["Binder"]].copy()
    b["_k"] = b["Bind_Consensus_%Rank"].fillna(1e9)
    best = b.sort_values(["_k", "Agree_Fraction"], ascending=[True, False]).drop_duplicates(KEYS)
    core_col = [c for c in ("NII_Core", "MIX2_Core") if c in best.columns][:1]
    best = best[KEYS + ["Allele", "Bind_Consensus_%Rank", "Agree_Fraction", "Imm_Votes", "Imm_Tools", "Proc_Votes"] + core_col].rename(
        columns={"Allele": "Best_Allele", "Bind_Consensus_%Rank": "Bind_Consensus_Best", "Agree_Fraction": "Agree_Fraction_Best",
                 "Imm_Votes": "Imm_Votes_Best", "Imm_Tools": "Imm_Tools_Best", "Proc_Votes": "Proc_Votes_Best", **{c: "Best_Core" for c in core_col}})
    out = out.merge(best, on=KEYS, how="left")
    if "Best_Core" in out:                       # where the binding core sits in the protein (class II: the 9-aa core in a longer window)
        off = [pep.find(core) if isinstance(core, str) else -1 for pep, core in zip(out["Peptide"], out["Best_Core"])]
        out["Core_Start"] = [st + o if o >= 0 else np.nan for st, o in zip(out["Start"], off)]
        out["Core_End"] = [cs + len(core) - 1 if cs == cs else np.nan for cs, core in zip(out["Core_Start"], out["Best_Core"].fillna(""))]
    immb = b.groupby(KEYS, sort=False)["Imm_TopPct_Consensus"].min().rename("Imm_TopPct_Best").reset_index()
    out = out.merge(immb, on=KEYS, how="left")

    # population-weighted agreement / binding over every tested allele that has data
    w = d["Allele_Freq"] if has_pop else pd.Series(1.0, index=d.index)
    d["_w"] = w
    d["_wa"] = d["Agree_Fraction"] * w
    d["_lr"] = np.log10(d["Bind_Consensus_%Rank"].clip(lower=1e-3)) * w
    d["_wl"] = w.where(d["Bind_Consensus_%Rank"].notna(), 0.0)
    d["_wf"] = w.where(d["Agree_Fraction"].notna(), 0.0)
    sums = d.groupby(KEYS, sort=False)[["_wa", "_lr", "_wl", "_wf"]].sum().reset_index()
    out = out.merge(sums, on=KEYS, how="left")
    out["Agree_Fraction_PopW"] = (out["_wa"] / out["_wf"].where(out["_wf"] > 0)).round(4)
    out["Bind_Consensus_PopW"] = (10 ** (out["_lr"] / out["_wl"].where(out["_wl"] > 0))).round(4)
    out = out.drop(columns=["_wa", "_lr", "_wl", "_wf"])
    d.drop(columns=["_w", "_wa", "_lr", "_wl", "_wf"], inplace=True)

    # coverage of the binder alleles and of the candidate alleles
    def cov_frame(col, prefix):
        sets = d[d["Binder" if prefix == "Binder" else "Pair_Candidate"]].groupby(KEYS, sort=False)["Allele"].agg(lambda a: tuple(sorted(a)))
        cache, rows = {}, []
        for key, alleles in sets.items():
            if alleles not in cache:
                cache[alleles] = pop.per_region(list(alleles)) if has_pop else {"national": len(alleles) / max(n_alleles, 1)}
            rows.append((*key, *(cache[alleles].get(r, np.nan) for r in ["national"] + (pop.regions if has_pop else []))))
        cols = KEYS + [col] + ([f"{col}_{r}" for r in pop.regions] if has_pop else [])
        return pd.DataFrame(rows, columns=cols)
    out = out.merge(cov_frame("Coverage_Binder", "Binder"), on=KEYS, how="left")
    out = out.merge(cov_frame("Coverage_Candidate", "Candidate"), on=KEYS, how="left")
    for c in [c for c in out.columns if c.startswith("Coverage_")]:
        out[c] = out[c].fillna(0.0).round(5)
    out.attrs["coverage_metric"] = "population" if has_pop else "allele fraction"

    # candidate rule (peptide level): enough alleles, enough coverage
    out["Candidate"] = ((out["N_Alleles_Candidate"] >= max(1, int(rules["min_alleles"]))) &
                        (out["Coverage_Candidate"] >= float(rules["min_coverage"])))

    # host proteome: flag (and optionally drop) candidates that look like the host's own proteins
    notes = d.attrs.setdefault("notes", [])
    out = hostself.annotate(out, rules, notes)

    # priority: weighted mean of percentile standing on four components (1 = best of this run)
    def pct(s, higher_better):
        return s.rank(pct=True, ascending=higher_better, method="average")
    comp = {"coverage": pct(out["Coverage_Binder"], True),
            "agreement": pct(out["Agree_Fraction_PopW"], True) if out["Agree_Fraction_PopW"].notna().any() else None,
            "binding": pct(out["Bind_Consensus_PopW"], False) if out["Bind_Consensus_PopW"].notna().any() else None,
            "immunogenicity": pct(out["Imm_TopPct_Best"], False) if out["Imm_TopPct_Best"].notna().any() else None}
    wts = {k: float(rules["rank_weights"].get(k, 0)) for k, v in comp.items() if v is not None}
    tot = sum(wts.values()) or 1.0
    score = sum(comp[k].fillna(0) * wts[k] for k in wts) / tot
    out["Priority_Score"] = (score * 100).round(2)
    mode = rules["rank_by"]
    if mode == "consensus":
        key = [out["Priority_Score"]]
    elif mode == "coverage":
        key = [out["Coverage_Binder"], out["Priority_Score"]]
    elif mode == "promiscuity":
        key = [out["N_Alleles_Binder"], out["Coverage_Binder"], out["Priority_Score"]]
    elif mode == "binding":
        key = [-np.log10(out["Bind_Consensus_PopW"].clip(lower=1e-3)).fillna(-99), out["Priority_Score"]]
    else:
        key = [-(out["Imm_TopPct_Best"].fillna(101)), out["Priority_Score"]]
    order = pd.DataFrame({"c": out["Candidate"].astype(int), **{f"k{i}": k.to_numpy() for i, k in enumerate(key)}})
    order = order.sort_values(list(order.columns), ascending=False, kind="stable")
    out = out.loc[order.index].reset_index(drop=True)
    out.insert(0, "Priority_Rank", range(1, len(out) + 1))

    out = out.merge(peptide_props.describe(out["Peptide"].drop_duplicates().tolist()), on="Peptide", how="left")
    return out


# ---------------------------------------------------------------------------------------------- shortlist
def _spans(cands, by_core):
    """(start, end) used to decide what clusters together: the binding core where known and asked for, else the whole peptide."""
    s, e = cands["Start"].to_numpy(float), cands["End"].to_numpy(float)
    if by_core and "Core_Start" in cands:
        cs, ce = cands["Core_Start"].to_numpy(float), cands["Core_End"].to_numpy(float)
        ok = ~np.isnan(cs)
        s, e = np.where(ok, cs, s), np.where(ok, ce, e)
    return s, e


def collapse_overlaps(cands, min_overlap, by_core=False):
    """Group candidates (already best-first) whose spans overlap the representative by >= min_overlap of the
    shorter one. Anchored on the representative's own span (not the growing cluster), so a region cannot creep along
    the antigen. Same protein only. With by_core, class II windows are compared by their 9-aa binding cores, so the many
    15-25-mers built around one core fall into one cluster."""
    clusters = []
    cs, ce = _spans(cands, by_core)
    for i, r in enumerate(cands.itertuples(index=False)):
        placed = False
        for c in clusters:
            if c["protein"] != r.Protein:
                continue
            a0, a1 = c["anchor"]
            ov = min(ce[i], a1) - max(cs[i], a0) + 1
            shorter = min(ce[i] - cs[i] + 1, a1 - a0 + 1)
            if shorter > 0 and ov / shorter >= min_overlap:
                c["members"].append(i)
                c["lo"], c["hi"] = min(c["lo"], r.Start), max(c["hi"], r.End)
                placed = True
                break
        if not placed:
            clusters.append({"protein": r.Protein, "anchor": (cs[i], ce[i]), "lo": r.Start, "hi": r.End, "members": [i]})
    return clusters


def shortlist(pep, rules, pop, n_alleles_tested):
    """Marks cluster / representative / panel columns on the peptide table; returns (pep, panel list, curve)."""
    pep["Cluster_ID"] = pd.array([pd.NA] * len(pep), dtype="Int64")
    pep["Is_Representative"] = False
    pep["Cluster_Size"] = pd.array([pd.NA] * len(pep), dtype="Int64")
    pep["Region"] = ""
    pep["Cluster_Representative"] = ""        # the peptide that stands for this one (itself for a representative)
    pep["Cluster_Overlap"] = np.nan           # overlap with the representative, as a fraction of the shorter peptide
    pep["Cluster_Members"] = ""               # on a representative: the other peptides clustered into it, "SEQ@start"
    pep["In_Panel"] = False
    pep["Panel_Order"] = pd.array([pd.NA] * len(pep), dtype="Int64")
    pep["Panel_Cum_Coverage"] = np.nan
    cands = pep[pep["Candidate"]]
    if cands.empty:
        return pep, [], []
    if rules["collapse"]:
        clusters = collapse_overlaps(cands, float(rules["overlap_fraction"]), rules.get("cluster_by") == "core")
    else:
        clusters = [{"protein": r.Protein, "anchor": (r.Start, r.End), "lo": r.Start, "hi": r.End, "members": [i]}
                    for i, r in enumerate(cands.itertuples(index=False))]
    idx = list(cands.index)
    for cid, c in enumerate(clusters, 1):
        members = [idx[i] for i in c["members"]]
        pep.loc[members, "Cluster_ID"] = cid
        pep.loc[members, "Cluster_Size"] = len(members)
        pep.loc[members, "Region"] = f"{c['protein']}:{c['lo']}-{c['hi']}"
        pep.loc[members[0], "Is_Representative"] = True
        rep = pep.loc[members[0]]
        pep.loc[members, "Cluster_Representative"] = rep["Peptide"]
        for m in members:
            use_core = rules.get("cluster_by") == "core" and "Core_Start" in pep and pep.at[m, "Core_Start"] == pep.at[m, "Core_Start"] and rep["Core_Start"] == rep["Core_Start"]
            a0, a1, b0, b1 = ((pep.at[m, "Core_Start"], pep.at[m, "Core_End"], rep["Core_Start"], rep["Core_End"]) if use_core
                              else (pep.at[m, "Start"], pep.at[m, "End"], rep["Start"], rep["End"]))
            ov = min(a1, b1) - max(a0, b0) + 1
            pep.at[m, "Cluster_Overlap"] = round(ov / max(1, min(a1 - a0 + 1, b1 - b0 + 1)), 3)
        pep.at[members[0], "Cluster_Members"] = "; ".join(f"{pep.at[m, 'Peptide']}@{pep.at[m, 'Start']}" for m in members[1:])

    # greedy panel over the representatives: add the peptide that gains the most coverage, until the target is reached
    reps = pep[pep["Is_Representative"]]
    metric = (lambda S: pop.coverage(sorted(S))) if pop and pop_has_data(pop, pep) else (lambda S: len(S) / max(n_alleles_tested, 1))
    remaining = {i: set(a.strip() for a in r.Alleles_Candidate.split(",") if a.strip()) for i, r in reps.iterrows()}
    union, cov, panel, curve = set(), 0.0, [], []
    while remaining and cov < float(rules["panel_target"]) and len(panel) < int(rules["panel_cap"]):
        best_i, best_cov = None, cov
        for i, S in remaining.items():
            c2 = metric(union | S)
            if c2 > best_cov + 1e-12 or (best_i is None and c2 > cov):
                best_i, best_cov = i, c2
        if best_i is None:
            break
        added = remaining.pop(best_i) - union
        union |= added
        gain, cov = best_cov - cov, best_cov
        panel.append({"order": len(panel) + 1, "peptide": pep.at[best_i, "Peptide"], "protein": pep.at[best_i, "Protein"],
                      "start": int(pep.at[best_i, "Start"]), "length": int(pep.at[best_i, "Length"]),
                      "gain": round(float(gain), 5), "cumulative_coverage": round(float(cov), 5), "new_alleles": sorted(added),
                      "region": pep.at[best_i, "Region"]})
        pep.at[best_i, "In_Panel"] = True
        pep.at[best_i, "Panel_Order"] = len(panel)
        pep.at[best_i, "Panel_Cum_Coverage"] = round(float(cov), 5)
    # the coverage a long panel would reach (curve), for the chart: cumulative over the first panel_cap representatives in rank order
    u2, cc = set(), []
    for i in list(reps.index)[:int(rules["panel_cap"])]:
        u2 |= set(a.strip() for a in pep.at[i, "Alleles_Candidate"].split(",") if a.strip())
        cc.append({"n": len(cc) + 1, "peptide": pep.at[i, "Peptide"], "coverage": round(float(metric(u2)), 5)})
    return pep, panel, cc


def pop_has_data(pop, pep):
    return pop is not None


# ---------------------------------------------------------------------------------------------- the whole analysis
def analyse(pairs, cls, rules, pop=None, filters_ok=None, light=False):
    """-> dict(pairs=master pair table, peptides=peptide table, candidates=representatives in rank order, panel, curve, funnel, tools, rules)."""
    d = pairs.copy()
    if filters_ok is not None:
        d["Filters_OK"] = filters_ok.reindex(d.index).fillna(False)
        fok = d["Filters_OK"]
    else:
        fok = None
    d, tools = add_votes(d, cls, rules)
    d = add_calls(d, tools, rules, fok)
    pep = peptide_table(d, tools, rules, pop, cls)
    n_alleles = d["Allele"].nunique()
    pep, panel, curve = shortlist(pep, rules, pop, n_alleles)
    cand = pep[pep["Is_Representative"]].copy() if rules["collapse"] or True else pep[pep["Candidate"]]
    funnel = make_funnel(d, pep, rules, cls)
    return {"notes": d.attrs.get("notes", []), "pairs": d, "peptides": pep, "candidates": cand, "panel": panel, "curve": curve, "funnel": funnel,
            "tools": [{"key": t.key, "label": t.label, "group": t.group, "kind": t.kind, "col": t.col} for t in tools],
            "rules": rules, "coverage_metric": "population" if (pop is not None and (d["Allele_Freq"] > 0).any()) else "allele fraction"}


def make_funnel(d, pep, rules, cls):
    n = d[KEYS].drop_duplicates().shape[0]
    steps = [("Peptides scanned", n)]
    m1 = d["Binder"]
    steps.append(("Binder in at least one allele (tool agreement)", int(d[m1][KEYS].drop_duplicates().shape[0])))
    if d["Imm_Tools"].max() > 0 and rules["imm_min_agree"] != "off":
        m2 = m1 & d["Imm_OK"]
        steps.append(("... and immunogenic (tool agreement)", int(d[m2][KEYS].drop_duplicates().shape[0])))
    else:
        m2 = m1
    if rules["require_processing"] and d["Proc_Tools"].max() > 0:
        m3 = m2 & d["Proc_OK"]
        steps.append(("... and processed", int(d[m3][KEYS].drop_duplicates().shape[0])))
    else:
        m3 = m2
    m4 = m3 & d["Filters_OK"]
    if not bool(d["Filters_OK"].all()):
        steps.append(("... and your per-tool filters", int(d[m4][KEYS].drop_duplicates().shape[0])))
    steps.append((f"... in at least {int(rules['min_alleles'])} allele(s) and {float(rules['min_coverage']):.0%} coverage" if (rules["min_alleles"] > 1 or rules["min_coverage"] > 0) else "Candidates", int(pep["Candidate"].sum())))
    if rules["collapse"]:
        steps.append(("Regions (overlapping windows merged)", int(pep["Is_Representative"].sum())))
    return [{"label": a, "n": b} for a, b in steps]


# ---------------------------------------------------------------------------------------------- suggestions
def quick_counts(pairs, cls, rules, pop, filters_ok):
    r = analyse(pairs, cls, rules, pop, filters_ok)
    pep = r["peptides"]
    return {"candidates": int(pep["Candidate"].sum()), "regions": int(pep["Is_Representative"].sum()),
            "panel": len(r["panel"]), "panel_coverage": (r["panel"][-1]["cumulative_coverage"] if r["panel"] else 0.0)}


def suggestions(pairs, cls, rules, pop, filters_ok, current, max_pairs=250_000):
    """What the same data gives under each preset and under single-rule relaxations, plus plain-language advice."""
    out = {"presets": [], "relax": [], "advice": []}
    if len(pairs) > max_pairs:
        out["advice"].append("The what-if table is skipped for very large runs.")
        return out
    for name, rr in preset_rules(cls).items():
        c = current if rr == rules else quick_counts(pairs, cls, rr, pop, filters_ok)
        out["presets"].append({"preset": name, "label": PRESET_INFO[name][0], **c})
    tries = []
    if rules["bind_min_agree"] in ("all", "majority"):
        tries.append(("Accept a binding call from one tool fewer", {"bind_min_agree": "majority" if rules["bind_min_agree"] == "all" else "any"}))
    if rules["imm_min_agree"] not in ("off",):
        tries.append(("Do not require immunogenicity agreement", {"imm_min_agree": "off"}))
    if rules["imm_top_q"] < 0.5:
        tries.append((f"Widen the immunogenicity cut to the top {min(rules['imm_top_q'] * 2, 1):.0%}", {"imm_top_q": min(rules["imm_top_q"] * 2, 1.0)}))
    tries.append((f"Loosen the binding cut-off to {rules['bind_rank_max'] * 2:g} %", {"bind_rank_max": rules["bind_rank_max"] * 2}))
    if rules["min_alleles"] > 1:
        tries.append(("Require one allele fewer", {"min_alleles": int(rules["min_alleles"]) - 1}))
    if rules["min_coverage"] > 0:
        tries.append(("Drop the minimum coverage", {"min_coverage": 0.0}))
    for label, over in tries[:6]:
        rr = {**rules, **over}
        c = quick_counts(pairs, cls, rr, pop, filters_ok)
        out["relax"].append({"label": label, **c, "gain": c["candidates"] - current["candidates"]})
    out["relax"].sort(key=lambda x: -x["gain"])
    n = current["candidates"]
    if n == 0:
        best = max(out["relax"], key=lambda x: x["candidates"], default=None)
        out["advice"].append("No peptide passes these rules." + (f" The biggest single change is: {best['label'].lower()} ({best['candidates']} candidates)." if best and best["candidates"] else ""))
    elif current["regions"] > 60:
        out["advice"].append(f"{current['regions']} regions is a long list; the Strict preset is a better starting point for synthesis.")
    elif current["regions"] < 3 and any(p["regions"] >= 5 for p in out["presets"]):
        wide = next(p for p in out["presets"] if p["regions"] >= 5)
        out["advice"].append(f"Only {current['regions']} region(s). The {wide['label']} preset gives {wide['regions']}.")
    return out


def run_details(info, rules, tools, pop, warnings):
    return {**info, "rules": {k: v for k, v in rules.items()}, "tools_used": tools,
            "population": pop.info() if pop else None, "warnings": warnings}
