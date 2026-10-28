"""One form, both T-cell arms: class I (CTL) and class II (HTL) on the same protein(s), then the regions where they meet.

Each class runs its own pipeline (alleles, fragment sizes, models) in the same job; the protein is read once and the
windows of every size are cut once per class and de-duplicated before any model sees them. Afterwards the candidate
peptides of the two classes are grouped into labelled regions (CTL+HTL, CTL only, HTL only), and a single HTML page,
per-model CSVs and a zip of everything are written.
"""
import json
import os
import time
import zipfile

import pandas as pd

import consensus
import report
import html_report
import pipeline
import pipeline_ii
import regions as regions_mod
from common import config, population
from common.jobview import JobView
import netmhcpan_app as nm

MODULE = {1: pipeline, 2: pipeline_ii}
SHARED = ("text", "mode", "preset", "rules", "title")


def split(params):
    """-> {1: params for class I or None, 2: ... }. Accepts {"c1": {...}, "c2": {...}, shared fields} (flat single-class
    params are still understood by the old endpoints)."""
    out = {}
    for c in (1, 2):
        sub = params.get(f"c{c}")
        if not sub or sub.get("enabled") is False:
            out[c] = None
            continue
        p = {k: v for k, v in sub.items() if k != "enabled"}
        for k in SHARED:
            if k in params:
                p.setdefault(k, params[k])
        p["mhc_class"] = c
        out[c] = p
    return out


def normalise(params):
    parts = split(params)
    if not any(parts.values()):
        raise ValueError("Switch on MHC class I, class II, or both.")
    if params.get("mode", "scan") == "peptides" and parts[2] is not None and parts[1] is None:
        pass
    return {c: (MODULE[c].normalise(p) if p else None) for c, p in parts.items()}, parts


def preview(params):
    norm, parts = normalise(params)
    out = {}
    for c, p in parts.items():
        if p:
            out[f"c{c}"] = MODULE[c].preview(p)
    tot = sum(v["n_pairs_netmhcpan"] + v["n_pairs_fast"] + v["n_pairs_bigmhc"] for v in out.values())
    warnings = list(dict.fromkeys(w for v in out.values() for w in v.get("warnings", [])))
    return {"classes": out, "warnings": warnings, "n_predictions": tot,
            "too_big": any(v.get("too_big") for v in out.values())}


def _seqs(text, mode):
    return dict(nm.read_sequences(text)) if mode == "scan" else {}


def run(params, job):
    t0 = time.time()
    norm, parts = normalise(params)
    classes = [c for c in (1, 2) if norm[c]]
    share = 1.0 / len(classes)
    views = {c: JobView(job, f"c{c}", "MHC I (CTL)" if c == 1 else "MHC II (HTL)", share) for c in classes}
    seqs = _seqs(parts[classes[0]].get("text", ""), parts[classes[0]].get("mode", "scan"))

    for c in classes:
        job.log(f"--- MHC class {'I' if c == 1 else 'II'} ---")
        MODULE[c].run(parts[c], views[c])
    _finalize(job, params, classes, {c: views[c].result for c in classes}, seqs, t0, parts[classes[0]].get("preset", "balanced"))


def _finalize(job, params, classes, class_results, seqs, t0, strategy):
    """Regions, report, per-model files and the zip, from the per-class results (fresh or re-ranked)."""
    job.stage("regions", "Regions and report", weight=0.3)
    job.update("regions", msg="grouping candidates into regions", state="running")

    results = {}
    for c in classes:
        r = class_results[c]
        results[c] = {"result": r, "peptides": pd.read_csv(r["peptides_csv"], low_memory=False), "consensus": r["consensus"],
                      "details": json.load(open(r["details_json"]))}
    rules = results[classes[0]]["consensus"]["rules"]
    reg = regions_mod.build(results.get(1, {}).get("peptides"), results.get(2, {}).get("peptides"), seqs,
                            float(rules.get("overlap_fraction", 0.7)))
    for c in classes:
        lab = regions_mod.label_peptides(results[c]["peptides"], reg, c)
        results[c]["peptides"] = lab
        lab.round(6).to_csv(results[c]["result"]["peptides_csv"], index=False)

    jobdir = config.run_dir(job.id)
    paths = {"regions_csv": os.path.join(jobdir, f"{job.id}_regions.csv"), "report_html": os.path.join(jobdir, f"{job.id}_report.html"),
             "bundle_zip": os.path.join(jobdir, f"{job.id}_all_outputs.zip")}
    reg.round(6).to_csv(paths["regions_csv"], index=False)

    warns = list(dict.fromkeys(w for c in classes for w in results[c]["result"]["summary"]["warnings"]))
    pop, _ = population.load(config.load().get("population"), config.search_root())
    info_text = "; ".join(
        f"class {'I' if c == 1 else 'II'}: {len(results[c]['details']['input']['alleles'])} alleles, fragment sizes "
        f"{', '.join(map(str, results[c]['details']['input'].get('lengths') or [])) or 'as given'}, models "
        f"{', '.join(results[c]['details']['input']['models'])}" for c in classes)
    html = html_report.build({"title": params.get("title") or "Epitope prediction report", "job": job.id, "strategy": strategy,
                              "regions": reg, "classes": results, "seqs": seqs, "warnings": warns, "rules": rules, "info_text": info_text,
                              "population_text": (f"{pop.label}: {len(pop.freqs)} alleles, regions weighted " +
                                                  ", ".join(f"{k} {v:.2f}" for k, v in pop.region_weights.items())) if pop else "none: coverage is the share of the chosen alleles"})
    with open(paths["report_html"], "w", encoding="utf-8") as f:
        f.write(html)

    permodel = _per_model(results, jobdir, job.id)
    models_by_key = {f"{c}:{t['key']}": path for c in classes for t in results[c]["consensus"]["tools"] for arc, path in permodel if path.endswith(f"_c{c}_{t['key']}.csv")}
    with zipfile.ZipFile(paths["bundle_zip"], "w", zipfile.ZIP_DEFLATED) as z:
        z.write(paths["report_html"], "report.html")
        z.write(paths["regions_csv"], "regions.csv")
        for c in classes:
            r, tag = results[c]["result"], f"class{'I' * c if c < 3 else ''}"
            for key, name in (("master_csv", "master"), ("peptides_csv", "peptides"), ("candidates_csv", "candidates"), ("panel_csv", "panel"), ("details_json", "run_details")):
                z.write(r[key], f"{tag}/{name}.{'json' if key.endswith('json') else 'csv'}")
        for arc, path in permodel:
            z.write(path, arc)

    n_dual = int((reg["Type"] == "CTL+HTL").sum()) if len(reg) else 0
    job.finish_stage("regions", f"{len(reg)} regions, {n_dual} CTL+HTL")
    s = {"peptides": sum(results[c]["result"]["summary"]["peptides"] for c in classes),
         "pairs": sum(results[c]["result"]["summary"]["pairs"] for c in classes), "regions": int(len(reg)), "regions_dual": n_dual,
         "seconds": round(time.time() - t0, 1), "warnings": warns, "chained": any(results[c]["result"]["summary"]["chained"] for c in classes),
         "pairs_passed": 0, "peptides_passed": 0}
    job.result = {"dual": True, "classes": {str(c): results[c]["result"] for c in classes}, **paths, "summary": s,
                  "consensus": {str(c): results[c]["consensus"] for c in classes}, "permodel": models_by_key}
    job.log(f"done in {s['seconds']}s: {len(reg)} regions ({n_dual} CTL+HTL)")


def _per_model(results, jobdir, jid):
    """One CSV per tool: peptide x allele rows with only that tool's columns (and the vote)."""
    out = []
    for c, d in results.items():
        master = pd.read_csv(d["result"]["master_csv"], low_memory=False)
        for t in d["consensus"]["tools"]:
            if t["col"] not in master:
                continue
            cols = [k for k in ["Protein", "Start", "Peptide", "Length", "Allele"] if k in master]
            extra = [k for k in master.columns if k.startswith(t["col"].split("_")[0] + "_") and k not in cols and not k.startswith(("Vote_", "TopPct_"))]
            keep = cols + list(dict.fromkeys([t["col"]] + extra)) + [k for k in (f"Vote_{t['key']}", f"TopPct_{t['key']}") if k in master]
            df = master[keep].dropna(subset=[t["col"]])
            path = os.path.join(jobdir, f"{jid}_c{c}_{t['key']}.csv")
            df.to_csv(path, index=False)
            out.append((f"per_model/class{'I' * c}_{t['label']}.csv", path))
    return out


LOOSER = {"any": 0, "majority": 1, "all": 2}


def rerank(src_job, params, job):
    """Same predictions, another strategy: rebuilds votes, coverage, ranking, regions and the report from the saved pair tables
    (seconds, no model runs). Pairs that a slow tool never scored because the original strategy ruled them out stay unscored."""
    t0 = time.time()
    src = src_job.result
    if not src or not src.get("dual"):
        raise ValueError("Only finished pipeline runs can be re-ranked.")
    classes = sorted(int(c) for c in src["classes"])
    pop_job_text = params.get("text", "")
    seqs = _seqs(pop_job_text, params.get("mode", "scan"))
    job.stage("rank", "Re-ranking the saved predictions", weight=1.0, total=len(classes))
    out, warns = {}, []
    for c in classes:
        r0 = src["classes"][str(c)]
        if not r0.get("base_pkl") or not os.path.exists(r0["base_pkl"]):
            raise ValueError("This run was made before re-ranking existed: its pair table was not kept. Run it again once.")
        base = pd.read_pickle(r0["base_pkl"])
        fok = base.pop("_filters_ok").astype(bool)
        details0 = json.load(open(r0["details_json"]))
        p = {"preset": params.get("preset", "balanced"), "rules": params.get("rules") or {}}
        new = consensus.resolve_rules(c, p["preset"], p["rules"])
        old = details0["rules"]
        looser = (new["bind_rank_max"] > old["bind_rank_max"] or LOOSER.get(str(new["bind_min_agree"]), 1) < LOOSER.get(str(old["bind_min_agree"]), 1))
        if looser and details0["input"].get("funnel", {}).get("mode", "auto") == "auto":
            warns.append(f"Class {'I' if c == 1 else 'II'}: this strategy is looser than the one the run was made with. The slow tools (BigMHC, CD4episcore) "
                         "scored only pairs that could pass the original rules, so some new binders have no slow-tool score. Run again for exact numbers.")
        view = JobView(job, f"c{c}", "MHC I (CTL)" if c == 1 else "MHC II (HTL)")
        res = report.finish(base, c, p, view, details0["input"], list(details0.get("warnings", [])), fok)
        res["summary"] = {**r0["summary"], "warnings": list(dict.fromkeys(details0.get("warnings", []) + warns)), "seconds": round(time.time() - t0, 1)}
        out[c] = res
        job.update("rank", done=len(out), msg=f"{len(out)}/{len(classes)} classes")
    job.finish_stage("rank", "done")
    _finalize(job, params, classes, out, seqs, t0, params.get("preset", "balanced"))

