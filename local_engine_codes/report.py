"""Shared last step of both pipelines: consensus analysis, the master CSV and its companions, run details."""
import json
import os
import time

import pandas as pd

import consensus
from common import config, population

KEYS = ["Protein", "Start", "Peptide"]
PEP_ORDER_COLS = ["Priority_Rank", "Candidate", "Is_Representative", "In_Panel"]


def _round(df):
    return df.round(6)


def finish(base, cls, p, job, info, warns, filters_ok):
    """base: the pair table (all tool columns). Writes CSVs, returns the dict stored in job.result."""
    pop, msg = population.load(config.load().get("population"), config.search_root())
    if msg:
        warns.append(msg)
    rules = consensus.resolve_rules(cls, p.get("preset", "balanced"), p.get("rules"))
    res = consensus.analyse(base, cls, rules, pop, filters_ok)
    cur = {"candidates": int(res["peptides"]["Candidate"].sum()), "regions": int(res["peptides"]["Is_Representative"].sum()),
           "panel": len(res["panel"]), "panel_coverage": (res["panel"][-1]["cumulative_coverage"] if res["panel"] else 0.0)}
    sugg = consensus.suggestions(base, cls, rules, pop, filters_ok, cur)
    if pop:
        alleles = sorted(base["Allele"].unique())
        rep = pop.representativeness(alleles)
        loci = [r["locus"] for r in rep if r["locus"] in ("A", "B", "C", "DRB1")]
        miss = pop.missing_top(alleles, loci)
        sugg["alleles"] = {"representativeness": rep, "add": miss}
        for r in rep:
            if r["locus"] in ("A", "B", "C", "DRB1") and r["haplotype_share"] < 0.6:
                sugg["advice"].append(f"The chosen {r['locus']} alleles cover {r['haplotype_share']:.0%} of the population's {r['locus']} haplotypes; "
                                      "coverage numbers will understate what the whole population presents. See 'Alleles worth adding' below.")
        no_data = [a for a in alleles if pop.allele_info(a)[1] == "none"]
        if no_data:
            warns.append(f"No population data for {len(no_data)} allele(s) ({', '.join(no_data[:4])}{'...' if len(no_data) > 4 else ''}); "
                         "they count in agreement but not in coverage.")

    warns += res["notes"]
    pairs, pep = res["pairs"], res["peptides"]
    pcols = ["Priority_Rank", "Priority_Score", "Candidate", "Is_Representative", "Cluster_ID", "Region", "In_Panel", "Panel_Order",
             "Coverage_Binder", "Coverage_Candidate", "N_Alleles_Binder", "Agree_Fraction_PopW", "Bind_Consensus_PopW"]
    master = pairs.merge(pep[KEYS + pcols], on=KEYS, how="left")
    master = master.sort_values(["Priority_Rank", "Binder", "Bind_Consensus_%Rank", "Allele"],
                                ascending=[True, False, True, True], na_position="last", kind="stable").reset_index(drop=True)
    lead = ["Priority_Rank", "Protein", "Start", "Peptide", "Length", "Allele", "Allele_Freq", "Pair_Call", "Pair_Candidate",
            "Binder", "Strong", "Bind_Consensus_%Rank", "Agree_Votes", "Agree_Tools", "Agree_Fraction", "Agree_List"]
    master = master[[c for c in lead if c in master] + [c for c in master if c not in lead]]

    jobdir = os.path.join(config.tmp_dir(), "jobs")
    os.makedirs(jobdir, exist_ok=True)
    out = {}
    for name, df in (("master", master), ("peptides", pep), ("candidates", res["candidates"])):
        path = os.path.join(jobdir, f"{job.id}_{name}.csv")
        _round(df).to_csv(path, index=False)
        out[f"{name}_csv"] = path
    panel_path = os.path.join(jobdir, f"{job.id}_panel.csv")
    pd.DataFrame([{**x, "new_alleles": ", ".join(x["new_alleles"])} for x in res["panel"]]).to_csv(panel_path, index=False)
    out["panel_csv"] = panel_path

    details = {"job": job.id, "mhc_class": cls, "preset": p.get("preset", "balanced"), "finished": time.strftime("%Y-%m-%d %H:%M:%S"),
               "input": info, "rules": rules, "tools": res["tools"], "funnel": res["funnel"],
               "population": pop.info() if pop else None, "coverage_metric": res["coverage_metric"], "warnings": list(dict.fromkeys(warns))}
    dpath = os.path.join(jobdir, f"{job.id}_run_details.json")
    with open(dpath, "w") as f:
        json.dump(details, f, indent=2, default=str)
    out["details_json"] = dpath
    out["consensus"] = {"rules": rules, "funnel": res["funnel"], "panel": res["panel"], "curve": res["curve"], "suggestions": sugg,
                        "coverage_metric": res["coverage_metric"], "tools": res["tools"], "headline": cur,
                        "population": pop.info() if pop else None}
    return out
