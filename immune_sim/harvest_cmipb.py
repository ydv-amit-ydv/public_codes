#!/usr/bin/env python3
"""
Harvest real human vaccine-response time courses from CMI-PB.

CMI-PB (Computational Models of Immunity - Pertussis Boost, NIH-funded) openly
publishes longitudinal post-booster data: antibody titers, cell frequencies
and cytokines sampled before and at intervals after a Tdap boost in adults who
were primed in infancy with either acellular (aP) or whole-cell (wP) vaccine.

This is the wet-lab data the dynamic model needs. Everything else in this
project validates epitope PREDICTION; these are measured RESPONSE KINETICS -
how fast antibody rises after a boost, how large the fold-rise is, how cell
populations move - which is what dynamics.py actually simulates.

Important scope note: every CMI-PB subject was primed in infancy, so each
time course is a SECONDARY (memory recall) response, not a primary one. It
validates the boost arm of the model. It cannot validate a primary response,
and it says nothing about platform differences - Tdap is one formulation.

API: PostgREST at https://www.cmi-pb.org/api/v5/ , anonymous GET, no key.

    python3 harvest_cmipb.py
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

API = "https://www.cmi-pb.org/api/v5"
OUT = Path(__file__).parent / "data"
OUT.mkdir(exist_ok=True)
PAGE = 5000
PAUSE = 0.3


def fetch(endpoint, params=None):
    rows, offset = [], 0
    while True:
        p = dict(params or {}, limit=PAGE, offset=offset)
        url = f"{API}/{endpoint}?" + urllib.parse.urlencode(p, quote_via=urllib.parse.quote)
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        for attempt in range(4):
            try:
                with urllib.request.urlopen(req, timeout=120) as r:
                    batch = json.loads(r.read().decode())
                break
            except Exception as e:
                if attempt == 3:
                    print(f"    FAILED {endpoint} offset={offset}: {e}", file=sys.stderr)
                    return rows
                time.sleep(2 * (attempt + 1))
        if not batch:
            break
        rows.extend(batch)
        if len(batch) < PAGE:
            break
        offset += PAGE
        time.sleep(PAUSE)
    return rows


def main():
    print("Harvesting CMI-PB longitudinal vaccine-response data\n")

    tables = {}
    for name in ("subject", "specimen", "plasma_ab_titer", "pbmc_cell_frequency"):
        rows = fetch(name)
        tables[name] = pd.DataFrame(rows)
        print(f"  {name:24s} {len(rows):7,d} rows")

    subj, spec = tables["subject"], tables["specimen"]
    ab = tables["plasma_ab_titer"]
    if ab.empty or spec.empty:
        sys.exit("no data returned - the API may have moved")

    # Join titers to their sampling day and subject.
    df = (ab.merge(spec, on="specimen_id", how="left")
            .merge(subj, on="subject_id", how="left"))

    day = ("actual_day_relative_to_boost"
           if "actual_day_relative_to_boost" in df.columns
           else "planned_day_relative_to_boost")
    df["day"] = pd.to_numeric(df[day], errors="coerce")
    df = df.dropna(subset=["day", "subject_id", "antigen", "isotype"])

    value = "MFI_normalised" if "MFI_normalised" in df.columns else "MFI"
    df["titer"] = pd.to_numeric(df[value], errors="coerce")
    df = df.dropna(subset=["titer"])

    keep = [c for c in ["subject_id", "specimen_id", "day", "antigen", "isotype",
                        "is_antigen_specific", "titer", "unit", "infancy_vac",
                        "biological_sex", "dataset", "specimen_type", "visit"]
            if c in df.columns]
    out = df[keep].sort_values(["subject_id", "antigen", "isotype", "day"])
    f = OUT / "cmipb_ab_titer.csv"
    out.to_csv(f, index=False)
    print(f"\n  antibody time courses -> {f}")
    print(f"  {len(out):,} measurements, {out.subject_id.nunique()} subjects, "
          f"{out.antigen.nunique()} antigens")
    print(f"  sampling days: {sorted(out.day.unique())[:14]}")
    if "infancy_vac" in out.columns:
        print(f"  priming groups: {out.infancy_vac.value_counts().to_dict()}")

    cf = tables["pbmc_cell_frequency"]
    if not cf.empty:
        cf = (cf.merge(spec, on="specimen_id", how="left")
                .merge(subj, on="subject_id", how="left"))
        cf["day"] = pd.to_numeric(cf.get(day), errors="coerce")
        cf = cf.dropna(subset=["day", "cell_type_name"])
        kc = [c for c in ["subject_id", "specimen_id", "day", "cell_type_name",
                          "percent_live_cell", "infancy_vac", "dataset"]
              if c in cf.columns]
        f2 = OUT / "cmipb_cell_frequency.csv"
        cf[kc].sort_values(["subject_id", "cell_type_name", "day"]).to_csv(f2, index=False)
        print(f"\n  cell frequency time courses -> {f2}")
        print(f"  {len(cf):,} measurements, "
              f"{cf.cell_type_name.nunique()} cell types")


if __name__ == "__main__":
    main()
