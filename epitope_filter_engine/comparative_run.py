"""Run a comparison: score the panel and the query peptides, then read the queries against the panel.

The panel and the queries go through **the same pipeline, in the same run**, with the same
alleles, models and settings. That is the whole point: a query is only comparable with the
panel if nothing differs between them except the peptide. Scoring them in separate runs, or
against stored numbers from another day, would compare the peptides and the settings at once.

Peptides are used exactly as given (mode "peptides"), never cut: a validated epitope is a
specific peptide, and fragmenting it would compare something else.

One analysis per MHC class, so a peptide that could belong to either is never judged by a
mixed baseline. comparative.py holds the parsing and the statistics and knows nothing about
jobs; this module is only the plumbing.
"""
import json
import os
import re
import time

import pandas as pd

import comparative as cp
import pipeline
import pipeline_ii
from common import config
from common.jobview import JobView

MODULE = {1: pipeline, 2: pipeline_ii}
BOUNDS = {1: (8, 14), 2: (9, 30)}                  # loose-peptide lengths each class accepts
LABEL = {1: "MHC I (CTL)", 2: "MHC II (HTL)"}
MAX_PEPTIDES = 4000


def _queries(text):
    out, seen = [], set()
    for tok in re.split(r"[\s,;]+", (text or "").upper()):
        tok = tok.strip()
        if tok and tok not in seen:
            seen.add(tok)
            out.append(tok)
    return out


def normalise(params):
    """Validate the request and build one class-pipeline parameter set per class."""
    panel, warns = cp.parse_panel(params.get("panel") or "", default_class=params.get("mhc_class"))
    queries = _queries(params.get("queries"))
    if not len(panel):
        raise ValueError("The panel is empty. Give one peptide per line: MHC class, epitope, "
                         "alleles (optional), outcome (optional).")
    if not queries:
        raise ValueError("Give at least one query peptide to compare with the panel.")
    bad = [q for q in queries if set(q) - cp.AA]
    if bad:
        raise ValueError(f"Query peptides with non-standard residues: {', '.join(bad[:4])}.")

    want = params.get("classes") or sorted({int(c) for c in panel["Class"]})
    want = [c for c in (1, 2) if c in {int(x) for x in want}]
    if not want:
        raise ValueError("No MHC class to run: the panel has none and none was chosen.")

    out = {}
    for c in want:
        pan = [p for p in panel[panel["Class"] == c]["Peptide"]
               if BOUNDS[c][0] <= len(p) <= BOUNDS[c][1]]
        qs = [q for q in queries if BOUNDS[c][0] <= len(q) <= BOUNDS[c][1]]
        if not pan:
            warns.append(f"class {'I' if c == 1 else 'II'}: no panel peptide of an accepted length "
                         f"({BOUNDS[c][0]}-{BOUNDS[c][1]} aa); skipped")
            continue
        if not qs:
            warns.append(f"class {'I' if c == 1 else 'II'}: no query peptide of an accepted length "
                         f"({BOUNDS[c][0]}-{BOUNDS[c][1]} aa); skipped")
            continue
        peps = sorted(set(pan) | set(qs))
        if len(peps) > MAX_PEPTIDES:
            raise ValueError(f"{len(peps):,} peptides for class {'I' if c == 1 else 'II'}: "
                             f"at most {MAX_PEPTIDES:,} (panel plus queries).")
        src = {k: v for k, v in (params.get(f"c{c}") or {}).items()}
        p = {**src, "text": "\n".join(peps), "mode": "peptides",
             "preset": params.get("preset", src.get("preset", "balanced")),
             "rules": params.get("rules", src.get("rules")) or {},
             "filters": {},                                  # nothing is filtered out: every panel peptide must keep a score
             "funnel": {"mode": "all"},                      # and every model must score every peptide, or the baseline has holes
             "chain": False}
        try:
            norm_p = MODULE[c].normalise(p)
        except ValueError as e:                          # name the class: the message alone ("Select at
            tag = "I" if c == 1 else "II"                # least one model") does not say which arm failed
            raise ValueError(f"MHC class {tag}: {e} (models that are not installed are switched off "
                             f"automatically, so check Settings, or untick this class)") from e
        out[c] = {"p": norm_p, "queries": qs, "n_panel": len(pan), "n_peptides": len(peps)}
    if not out:
        raise ValueError("Nothing to run: " + " ".join(warns))
    return {"panel": panel, "queries": queries, "classes": out, "warnings": warns}


def preview(params):
    norm = normalise(params)
    per = {}
    for c, d in norm["classes"].items():
        pan = norm["panel"][norm["panel"]["Class"] == c]
        per[f"c{c}"] = {"n_panel": d["n_panel"], "n_queries": len(d["queries"]), "n_peptides": d["n_peptides"],
                        "n_positive": int((pan["Outcome"] == "positive").sum()),
                        "n_negative": int((pan["Outcome"] == "negative").sum()),
                        "n_alleles": len(d["p"]["alleles"]),
                        "n_predictions": d["n_peptides"] * len(d["p"]["alleles"])}
    no_neg = [k for k, v in per.items() if v["n_negative"] == 0]
    warns = list(norm["warnings"])
    if no_neg:
        warns.append("The panel has no peptide with outcome 'negative', so there is nothing to compare the "
                     "positives against and no confidence can be given. Add rows ending in 'negative' "
                     "(peptides that were tested and did not work).")
    return {"classes": per, "warnings": warns,
            "n_predictions": sum(v["n_predictions"] for v in per.values())}


def run(params, job):
    t0 = time.time()
    norm = normalise(params)
    classes = sorted(norm["classes"])
    views = {c: JobView(job, f"c{c}", LABEL[c]) for c in classes}

    for c in classes:
        p = norm["classes"][c]["p"]
        for st in MODULE[c].stage_plan(p, MODULE[c].n_unique(p)):
            views[c].stage(st["name"], st["label"], st["weight"], work=st.get("work"))
    job.stage("compare", "Comparison with the panel", weight=1.0)

    for c in classes:
        job.log(f"--- MHC class {'I' if c == 1 else 'II'}: {norm['classes'][c]['n_peptides']} peptides "
                f"({norm['classes'][c]['n_panel']} panel + {len(norm['classes'][c]['queries'])} query) ---")
        MODULE[c].run(norm["classes"][c]["p"], views[c])

    job.update("compare", msg="reading the queries against the panel", state="running")
    out = {"classes": {}, "warnings": list(norm["warnings"])}
    for c in classes:
        pep = pd.read_csv(views[c].result["peptides_csv"], low_memory=False)
        pep["Class"] = c
        pan_c = norm["panel"][norm["panel"]["Class"] == c]
        # A baseline that quietly shrank is the dangerous failure: the numbers still look
        # reasonable. Nothing is filtered or gated in this run, so a panel peptide with no
        # score means a model returned nothing for it, and that has to be said out loud.
        scored_set = set(pep["Peptide"].astype(str).str.upper())
        missing = [p for p in pan_c["Peptide"] if p not in scored_set]
        if missing:
            tag = "I" if c == 1 else "II"
            job.log(f"class {tag}: {len(missing)} panel peptide(s) got no score")
            out["warnings"].append(
                f"class {tag}: {len(missing)} of {len(pan_c)} panel peptide(s) got no score and are NOT in the "
                f"baseline ({', '.join(missing[:5])}{'...' if len(missing) > 5 else ''}). The baseline is smaller "
                f"than the panel you gave, so read the verdicts with that in mind.")
        res = cp.analyse(pep, pan_c, norm["classes"][c]["queries"])
        res["n_panel_scored"] = int(pep["Peptide"].isin(set(norm["panel"]["Peptide"])).sum())
        res["alleles"] = norm["classes"][c]["p"]["alleles"]
        res["peptides_csv"] = views[c].result["peptides_csv"]
        res["master_csv"] = views[c].result.get("master_csv")
        out["classes"][str(c)] = res

    jobdir = config.run_dir(job.id)
    os.makedirs(jobdir, exist_ok=True)
    path = os.path.join(jobdir, f"{job.id}_comparative.json")
    with open(path, "w") as f:
        json.dump({k: v for k, v in out.items()}, f, indent=1, default=str)

    rows = []
    for c, res in out["classes"].items():
        for q in res["queries"]:
            rows.append({"Class": c, "Peptide": q["peptide"], "Length": q.get("length"),
                         "Likelihood_Ratio": q.get("likelihood_ratio"), "Verdict": q.get("verdict", q.get("error")),
                         "Primary_Score": (q.get("evidence", {}).get(res["primary"], {}) or {}).get("score"),
                         "Pct_Of_Positives": (q.get("evidence", {}).get(res["primary"], {}) or {}).get("pct_of_positives"),
                         "Pct_Of_Negatives": (q.get("evidence", {}).get(res["primary"], {}) or {}).get("pct_of_negatives"),
                         "N_Positive": res["panel"]["n_positive"], "N_Negative": res["panel"]["n_negative"]})
    qcsv = os.path.join(jobdir, f"{job.id}_comparison.csv")
    pd.DataFrame(rows).to_csv(qcsv, index=False)

    verdicts = {}
    for r in rows:
        verdicts[r["Verdict"]] = verdicts.get(r["Verdict"], 0) + 1
    job.finish_stage("compare", f"{len(rows)} query result(s)")
    job.result = {"comparative": out, "comparative_json": path, "comparison_csv": qcsv,
                  "classes": {c: {"peptides_csv": res["peptides_csv"], "master_csv": res["master_csv"]}
                              for c, res in out["classes"].items()},
                  "summary": {"queries": len(rows), "verdicts": verdicts,
                              "panel": int(len(norm["panel"])),
                              "classes": [int(c) for c in classes],
                              "seconds": round(time.time() - t0, 1),
                              "warnings": out["warnings"]}}
    job.log(f"done in {time.time() - t0:.1f}s: " + ", ".join(f"{v} {k}" for k, v in verdicts.items()))
