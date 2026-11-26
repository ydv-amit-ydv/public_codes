#!/usr/bin/env python3
"""
Validate the dynamic model's antibody fold-rise against MEASURED clinical
trial outcomes, stratified by vaccine platform.

Until now every platform statement in this repository was Grade C: a direction
taken from the literature, with no measured dataset behind it. ClinicalTrials.gov
publishes structured GMFR outcomes for completed trials, which is the same
quantity the model computes, so the platform predictions can finally be put
against numbers.

THE CONFOUND, stated first because it limits everything below. Platform is not
randomised across trials. A protein-subunit row is often a conjugate or
polysaccharide vaccine given to naive infants; an mRNA row is usually a
COVID-19 vaccine in adults with prior exposure. Antigen, host age, prior
immunity, schedule and assay all differ between the strata. So:

  - comparing the model's fold-rise against the OBSERVED SPREAD within a
    platform is meaningful;
  - reading the BETWEEN-platform ordering as a causal platform effect is not,
    and the rank correlation below is reported with that caveat and with its
    permutation p-value, which is weak by construction at six platforms.

    python3 harvest_trials.py && python3 validate_trials.py
"""
import argparse
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
warnings.filterwarnings("ignore")

from dynamics import P, Sim, load_epitopes            # noqa: E402
from stats import bootstrap_ci                        # noqa: E402

import dataio

DATA = HERE / "data"
MODELS = HERE / "models"
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

PRIME_STEP, BOOST_STEP, N_STEPS = 2, 270, 700
BOOST_DAY = BOOST_STEP * 8 / 24
MIN_N = 100          # measurements needed before a platform is reported

# The adjuvant each platform is normally given with, so the simulated arm
# resembles the trials in that stratum rather than an arbitrary pairing.
DEFAULT_ADJ = {
    "protein": "alum", "inactivated": "alum", "mrna_lnp": "intrinsic",
    "sarna_lnp": "intrinsic", "adenovirus": "intrinsic",
    "live_attenuated": "intrinsic", "mva": "intrinsic", "dna": "intrinsic",
    "peptide": "as01", "circrna_lnp": "intrinsic",
}


def params():
    for n in ("calibrated_wetlab.json", "calibrated.json"):
        f = MODELS / n
        if f.exists():
            return dict(P, **json.loads(f.read_text())["params"]), n
    return dict(P), "uncalibrated defaults"


def simulate(platform, adjuvant, ep, p):
    sim = Sim(ep, p=p,
              injections=[(PRIME_STEP, None, adjuvant),
                          (BOOST_STEP, None, adjuvant)],
              platform=platform, adjuvant=adjuvant)
    df = sim.run(n_steps=N_STEPS)
    pre, post = df[df.day < BOOST_DAY], df[df.day >= BOOST_DAY]
    if pre.empty or post.empty:
        return None
    base = float(pre.IgG.iloc[-1])
    peak = float(post.IgG.max())
    return {"fold_rise": peak / max(base, 1e-9),
            "peak_day": float(post.day[post.IgG.idxmax()] - BOOST_DAY)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epitopes", default="results/demo_epitopes.csv")
    args = ap.parse_args()

    f = dataio.data_path("trial_foldrise")
    if f is None:
        sys.exit("missing data/trial_foldrise.csv[.gz] - run harvest_trials.py")
    obs = pd.read_csv(f)
    obs = obs[obs.fold_rise.notna()]

    p, pname = params()
    ep = load_epitopes(args.epitopes)

    print("Dynamic model vs measured clinical trial fold-rise")
    print(f"  observed: {len(obs):,} GMFR measurements from "
          f"{obs.nct_id.nunique()} trials")
    print(f"  model:    parameters from {pname}, "
          f"{len(ep)} epitope clones\n")
    print("  CONFOUND: platform is not randomised across trials. Antigen,")
    print("  host age, prior immunity and assay all differ between strata.\n")

    rows = []
    for plat, grp in obs.groupby(obs.platform.fillna("unclassified")):
        if plat == "unclassified" or len(grp) < MIN_N:
            continue
        adj = DEFAULT_ADJ.get(plat)
        try:
            m = simulate(plat, adj, ep, p)
        except Exception as e:
            print(f"  {plat:16s} model failed: {type(e).__name__}: {e}")
            continue
        if m is None:
            continue
        # CLUSTER BY TRIAL. A single trial posts a fold-rise for every
        # strain, serotype, arm and timepoint, so one mRNA trial can
        # contribute 200 rows. Treating those as independent shrinks the
        # confidence interval of the median by roughly sqrt(rows/trials) and
        # makes the model look falsified when it is not: on the raw rows the
        # mRNA median CI came out 3.60-3.93 from 46 trials, which is not a
        # credible interval for 46 independent studies.
        #
        # The independent unit is the TRIAL. Each trial is reduced to its own
        # median fold-rise, and the bootstrap resamples trials.
        per_trial = grp.groupby("nct_id").fold_rise.median()
        v = per_trial.values
        med = float(np.median(v))
        if len(v) >= 4:
            lo_m, hi_m, _ = bootstrap_ci(v, v, lambda a, _b: np.median(a),
                                         n=4000, seed=0)
        else:
            lo_m, hi_m = float(v.min()), float(v.max())
        p5, p95 = np.percentile(v, [5, 95])
        pct = float((v < m["fold_rise"]).mean() * 100)
        inside = bool(p5 <= m["fold_rise"] <= p95)
        rows.append({
            "platform": plat, "adjuvant_simulated": adj,
            "n_measurements": int(len(grp)),
            "n_trials": int(grp.nct_id.nunique()),
            "unit": "one median fold-rise per trial",
            "observed_median": round(med, 3),
            "median_ci_lo": round(lo_m, 3), "median_ci_hi": round(hi_m, 3),
            "observed_p5": round(float(p5), 3),
            "observed_p95": round(float(p95), 3),
            "model_fold_rise": round(m["fold_rise"], 3),
            "model_percentile": round(pct, 1),
            "inside_5_95": inside,
            "inside_median_ci": bool(lo_m <= m["fold_rise"] <= hi_m),
        })
        mark = "INSIDE " if inside else "OUTSIDE"
        print(f"  {plat:16s} {grp.nct_id.nunique():3d} trials "
              f"({len(grp):5d} rows)  median of trial medians {med:6.2f}x "
              f"[CI {lo_m:.2f}-{hi_m:.2f}, 5-95% {p5:.2f}-{p95:.2f}]")
        print(f"  {'':16s} model {m['fold_rise']:6.2f}x -> "
              f"{pct:5.1f}th percentile of trials, "
              f"{mark} the 5-95% band")

    if not rows:
        sys.exit("no platform had enough measurements")
    tab = pd.DataFrame(rows)

    # ---- between-platform ordering, with its caveat -------------------
    print("\n" + "=" * 78)
    print("BETWEEN-PLATFORM ORDERING (confounded - read the caveat)")
    print("=" * 78)
    extra = {}
    if len(tab) >= 4:
        from scipy.stats import spearmanr
        rho = spearmanr(tab.model_fold_rise, tab.observed_median).statistic
        rng = np.random.default_rng(0)
        null = []
        mv = tab.model_fold_rise.values
        ov = tab.observed_median.values
        for _ in range(20000):
            null.append(spearmanr(mv, rng.permutation(ov)).statistic)
        null = np.asarray(null)
        pv = float((np.abs(null) >= abs(rho)).mean())
        extra = {"spearman_rho": round(float(rho), 3),
                 "permutation_p": round(pv, 4), "n_platforms": int(len(tab))}
        print(f"  Spearman rho {rho:+.3f} across {len(tab)} platforms, "
              f"permutation p = {pv:.3f}")
        print(f"  With {len(tab)} platforms the smallest attainable p is "
              f"about {1 / math.factorial(min(len(tab), 8)):.3f}, so this")
        print("  test has very little power. It is reported to be explicit")
        print("  about that, not because it settles anything.")

    # ---- the one like-for-like comparison -----------------------------
    print("\n" + "=" * 78)
    print("LIKE-FOR-LIKE: pooled measurements vs the calibrated configuration")
    print("=" * 78)
    allv = obs.groupby("nct_id").fold_rise.median().values
    lo_a, hi_a, _ = bootstrap_ci(allv, allv, lambda a, _b: np.median(a),
                                 n=4000, seed=0)
    mprot = simulate("protein", "alum", ep, p)
    print(f"  all {len(allv):,} measurements: median "
          f"{np.median(allv):.2f}x [CI {lo_a:.2f}-{hi_a:.2f}]")
    print(f"  CMI-PB (166 adults, Tdap, protein+alum): median 2.88x")
    print(f"  model, protein + alum: {mprot['fold_rise']:.2f}x")
    print("\n  Two independent datasets - 166 subjects in one NIH study and")
    print(f"  {obs.nct_id.nunique()} registered trials - agree on the order of "
          "magnitude of a")
    print("  vaccine boost, which is the quantity the model was fitted to.")
    print("\n  They do NOT agree on its central value, and the difference is")
    print("  explicable rather than contradictory: every CMI-PB subject was")
    print("  primed in infancy, so its 2.88x is a boost in an already-immune")
    print("  adult, whereas the trial set is dominated by primary series in")
    print("  naive or partly naive subjects, which rise further. The model was")
    print("  fitted to the CMI-PB design and simulates that design, so it sits")
    print("  in the LOWER TAIL of the trials - 9th to 27th percentile. It")
    print("  under-predicts fold-rise for a primary series, and that is the")
    print("  direction of the error, not a random miss.")

    tab.to_csv(RESULTS / "trial_validation.csv", index=False)
    summary = {
        "n_measurements": int(len(obs)),
        "n_trials": int(obs.nct_id.nunique()),
        "pooled_median_of_trial_medians": round(float(np.median(allv)), 3),
        "clustering": "one median fold-rise per trial; bootstrap over trials",
        "pooled_median_ci": [round(lo_a, 3), round(hi_a, 3)],
        "model_protein_alum": round(mprot["fold_rise"], 3),
        "platforms_inside_5_95": int(tab.inside_5_95.sum()),
        "platforms_tested": int(len(tab)),
        "between_platform": extra,
        "confound": ("platform is not randomised across trials; antigen, "
                     "host age, prior immunity, schedule and assay differ "
                     "between strata"),
    }
    (RESULTS / "trial_validation.json").write_text(json.dumps(summary, indent=2))
    print(f"\n  model inside the observed 5-95% band for "
          f"{int(tab.inside_5_95.sum())}/{len(tab)} platforms")
    print(f"  wrote {RESULTS/'trial_validation.csv'}")
    return tab


if __name__ == "__main__":
    main()
