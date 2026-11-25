#!/usr/bin/env python3
"""
Harvest measured vaccine immunogenicity from ClinicalTrials.gov.

CMI-PB gave us one formulation: a Tdap boost, recombinant protein on alum, in
166 adults. Everything the model says about OTHER platforms - mRNA, viral
vector, inactivated, adjuvant choice - was a structured prior with no measured
data behind it at all.

ClinicalTrials.gov publishes STRUCTURED RESULTS for completed trials, and one
of the posted outcome measures is the quantity this model actually predicts:

    GMFR - geometric mean fold rise of antibody titre over the pre-vaccination
    baseline, at a stated timepoint.

That is directly comparable to `boost_fold_over_baseline` in calibrate.py. It
is reported per arm, which means the platform can be read off the intervention
and the comparison can be stratified: mRNA trials against protein trials
against inactivated trials, all measured.

What this does NOT give: cell counts, time courses, or anything about CD8.
Only antibody fold-rise and its timepoint.

    python3 harvest_trials.py
    python3 harvest_trials.py --max-studies 1200
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

HERE = Path(__file__).parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

API = "https://clinicaltrials.gov/api/v2/studies"
UA = "immune-sim-research/1.0"

FIELDS = ",".join([
    "protocolSection.identificationModule.nctId",
    "protocolSection.identificationModule.briefTitle",
    "protocolSection.conditionsModule.conditions",
    "protocolSection.armsInterventionsModule.interventions",
    "resultsSection.outcomeMeasuresModule",
])

# Search terms. Several, because a trial only turns up if its posted outcome
# text matches, and different sponsors word the same measure differently.
TERMS = ["geometric mean fold rise", "geometric mean fold increase",
         "GMFR", "fold rise in antibody", "geometric mean titer ratio"]

FOLD_RE = re.compile(r"fold[\s-]?(rise|increase|change)|GMFR|fold[\s-]?ratio",
                     re.I)
# Reject outcome measures that are a PROPORTION or a COUNT, not a fold-rise
BAD_UNIT = re.compile(r"particip|percent|%|subject|count|propor|titer$|"
                      r"^titre|log2|log10", re.I)

PLATFORM_PATTERNS = [
    ("mrna_lnp",       r"\bmRNA\b|BNT162|mRNA-1\d{3}|modRNA|"
                       r"lipid nanoparticle|\bLNP\b"),
    ("sarna_lnp",      r"self[\s-]?amplif|\bsaRNA\b|replicon RNA|ARCT-"),
    ("adenovirus",     r"\bAd26|\bAd5\b|ChAdOx|adenovir|Gam-COVID|"
                       r"Sputnik|Ad-vector"),
    ("mva",            r"\bMVA\b|modified vaccinia|\bNYVAC\b"),
    ("dna",            r"\bDNA vaccine\b|plasmid DNA|INO-\d+"),
    ("live_attenuated", r"live[\s-]?attenuat|\bLAIV\b|FluMist|"
                        r"\bMMR\b|varicella|yellow fever|rotavirus"),
    ("inactivated",    r"inactivat|whole[\s-]?virion|\bCoronaVac\b|"
                       r"BBIBP|split[\s-]?virion|Vero cell"),
    ("protein",        r"recombinant|subunit|conjugat|polysacchar|"
                       r"\bNVX\b|adjuvanted protein|\bVLP\b|toxoid|"
                       r"virus[\s-]?like particle"),
]
ADJUVANT_PATTERNS = [
    ("as01", r"\bAS01\b"), ("as03", r"\bAS03\b"), ("mf59", r"\bMF59\b"),
    ("cpg1018", r"CpG[\s-]?1018|\bCpG\b"), ("alum", r"alum|aluminium|aluminum|"
                                                    r"hydroxide|phosphate"),
    ("saponin", r"Matrix-M|QS-?21|saponin"),
]
DAY_RE = re.compile(r"(?:day|d)\s*(\d{1,3})\b", re.I)
WEEK_RE = re.compile(r"(\d{1,2})\s*week", re.I)
MONTH_RE = re.compile(r"(\d{1,2})\s*month", re.I)


def get(url, tries=4, timeout=180):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read())
        except Exception as e:
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"{url[:120]} failed: {last}")


def classify(text, patterns):
    for key, pat in patterns:
        if re.search(pat, text, re.I):
            return key
    return None


def timepoint_days(*texts):
    """Best guess at the days-post-dose a measure refers to."""
    for t in texts:
        if not t:
            continue
        m = DAY_RE.search(t)
        if m:
            d = int(m.group(1))
            if 1 <= d <= 400:
                return d
        m = WEEK_RE.search(t)
        if m:
            return int(m.group(1)) * 7
        m = MONTH_RE.search(t)
        if m:
            return int(m.group(1)) * 30
    return None


def parse_study(st):
    ps = st.get("protocolSection", {})
    nct = ps.get("identificationModule", {}).get("nctId")
    title = ps.get("identificationModule", {}).get("briefTitle", "")
    conds = "; ".join(ps.get("conditionsModule", {}).get("conditions", []))
    ivs = ps.get("armsInterventionsModule", {}).get("interventions", []) or []
    iv_text = " | ".join(
        f"{i.get('name', '')} {i.get('description', '')}" for i in ivs)
    blob = f"{title} {iv_text}"
    platform = classify(blob, PLATFORM_PATTERNS)
    adjuvant = classify(blob, ADJUVANT_PATTERNS)

    oms = (st.get("resultsSection", {})
             .get("outcomeMeasuresModule", {})
             .get("outcomeMeasures", []) or [])
    out = []
    for om in oms:
        t = om.get("title", "") or ""
        unit = om.get("unitOfMeasure", "") or ""
        if not FOLD_RE.search(t):
            continue
        # a fold-rise reported as a count or a percentage is not a fold-rise
        if om.get("paramType") not in ("GEOMETRIC_MEAN", "MEAN", "MEDIAN"):
            continue
        if BAD_UNIT.search(unit) and not re.search(r"fold", unit, re.I):
            continue
        tf = om.get("timeFrame", "") or ""
        for cls in om.get("classes", []) or []:
            ctitle = cls.get("title", "") or ""
            for cat in cls.get("categories", []) or []:
                for meas in cat.get("measurements", []) or []:
                    v = meas.get("value")
                    try:
                        val = float(v)
                    except (TypeError, ValueError):
                        continue
                    if not (1.0 <= val <= 2000.0):
                        continue       # a GMFR below 1 or above 2000 is a
                                       # different quantity mislabelled
                    out.append({
                        "nct_id": nct, "condition": conds[:90],
                        "platform": platform, "adjuvant": adjuvant,
                        "outcome_title": t[:140],
                        "unit": unit[:40],
                        "param_type": om.get("paramType"),
                        "class_title": ctitle[:70],
                        "group_id": meas.get("groupId"),
                        "fold_rise": val,
                        "day": timepoint_days(ctitle, tf, t),
                        "intervention": iv_text[:150],
                    })
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-studies", type=int, default=1500)
    ap.add_argument("--page-size", type=int, default=50)
    args = ap.parse_args()

    print("Harvesting measured vaccine fold-rise from ClinicalTrials.gov")
    rows, seen = [], set()
    for term in TERMS:
        token, n = None, 0
        print(f"\n  term: {term!r}")
        while n < args.max_studies // len(TERMS) + args.page_size:
            q = {"query.intr": "vaccine", "query.term": f'"{term}"',
                 "aggFilters": "results:with", "fields": FIELDS,
                 "pageSize": args.page_size}
            if token:
                q["pageToken"] = token
            try:
                d = get(f"{API}?{urllib.parse.urlencode(q)}")
            except Exception as e:
                print(f"    stopped: {e}")
                break
            studies = d.get("studies", []) or []
            if not studies:
                break
            for st in studies:
                nct = (st.get("protocolSection", {})
                         .get("identificationModule", {}).get("nctId"))
                if nct in seen:
                    continue
                seen.add(nct)
                rows.extend(parse_study(st))
            n += len(studies)
            print(f"    {n} studies scanned, {len(rows)} measurements",
                  flush=True)
            token = d.get("nextPageToken")
            if not token:
                break

    if not rows:
        sys.exit("no fold-rise outcome measures found")

    df = pd.DataFrame(rows).drop_duplicates()
    out = DATA / "trial_foldrise.csv"
    df.to_csv(out, index=False)

    print("\n" + "=" * 74)
    print(f"  {len(df)} fold-rise measurements from "
          f"{df.nct_id.nunique()} trials")
    print(f"  overall median {df.fold_rise.median():.2f}x, "
          f"IQR {df.fold_rise.quantile(.25):.2f}-"
          f"{df.fold_rise.quantile(.75):.2f}")
    print("\n  by platform (as read from the intervention text):")
    g = (df.groupby(df.platform.fillna("unclassified"))
           .fold_rise.agg(["size", "median",
                           lambda s: s.quantile(.25),
                           lambda s: s.quantile(.75)]))
    g.columns = ["n", "median", "q25", "q75"]
    for k, r in g.sort_values("n", ascending=False).iterrows():
        print(f"    {k:16s} n={int(r['n']):5d}  median {r['median']:7.2f}x  "
              f"IQR {r['q25']:.2f}-{r['q75']:.2f}")
    dd = df[df.day.notna()]
    if not dd.empty:
        print(f"\n  timepoint recovered for {len(dd)}/{len(df)} "
              f"measurements, median day {dd.day.median():.0f}")
    print(f"\n  wrote {out}")
    print("\n  Platform labels come from pattern-matching the intervention")
    print("  text, so they are noisy: an 'unclassified' row is one the")
    print("  patterns could not read, not a trial without a platform.")


if __name__ == "__main__":
    main()
