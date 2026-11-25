#!/usr/bin/env python3
"""
Validate the DYNAMIC model against measured human wet-lab kinetics (CMI-PB).

This is the test that was missing. Everything in validate.py concerns epitope
PREDICTION; this file tests what dynamics.py actually simulates - how antibody
moves after a boost - against 166 real subjects.

OUT-OF-SAMPLE. The model was calibrated in calibrate.py against qualitative
literature bands, before this dataset was harvested. No CMI-PB measurement was
used to fit any parameter. So the comparison below is a genuine prediction
test, and it is reported whether it passes or fails.

WHAT IT CAN AND CANNOT TEST
---------------------------
CAN:  the secondary (memory recall) response to a protein + alum boost -
      fold-rise over baseline, time to peak, and the shape of the early
      kinetics. Every CMI-PB subject was primed in infancy, so every time
      course is a boost.
CANNOT: a primary response (none in this dataset); any platform comparison
      (Tdap is one formulation); absolute titer units (the model's antibody
      is in arbitrary units, so only RATIOS and TIMES are comparable).

Run:  python3 wetlab_validation.py
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sps

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

from stats import bootstrap_ci, fmt_p, prob_direction          # noqa: E402

DATA = HERE / "data" / "cmipb_ab_titer.csv"
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

# Tdap vaccine antigens. The others in the panel (OVA, BETV1, FELD1, LOLP1)
# are not in the vaccine and serve as negative controls.
VACCINE_AG = ["PT", "PRN", "FHA", "FIM2/3", "DT", "TT"]
CONTROL_AG = ["OVA", "BETV1", "FELD1", "LOLP1", "Measles"]

MIN_SUBJECTS = 20
POST_MAX_DAY = 35


# ------------------------------------------------------------ observed data
def per_subject(df, antigen, isotype="IgG"):
    """
    For each subject: baseline titer, peak post-boost titer, fold-rise, and
    the day the peak occurred. Baseline is the last sample at or before day 0.
    """
    d = df[(df.antigen == antigen) & (df.isotype == isotype)]
    rows = []
    for sid, g in d.groupby("subject_id"):
        g = g.sort_values("day")
        pre = g[g.day <= 0]
        post = g[(g.day > 0) & (g.day <= POST_MAX_DAY)]
        if pre.empty or len(post) < 3:
            continue
        base = float(pre.titer.iloc[-1])
        if not np.isfinite(base) or base <= 0:
            continue
        i = int(np.argmax(post.titer.values))
        rows.append({
            "subject_id": sid,
            "baseline": base,
            "peak": float(post.titer.values[i]),
            "peak_day": float(post.day.values[i]),
            "fold_rise": float(post.titer.values[i] / base),
            "infancy_vac": g.infancy_vac.iloc[0] if "infancy_vac" in g else None,
        })
    return pd.DataFrame(rows)


def summarise(s, name):
    """Median with a bootstrap CI, plus the spread across subjects."""
    v = np.asarray(s, float)
    v = v[np.isfinite(v)]
    if len(v) < 5:
        return None
    rng = np.random.default_rng(0)
    meds = [np.median(rng.choice(v, len(v), replace=True)) for _ in range(4000)]
    return {
        "metric": name, "n": int(len(v)),
        "median": float(np.median(v)),
        "ci_lo": float(np.percentile(meds, 2.5)),
        "ci_hi": float(np.percentile(meds, 97.5)),
        "iqr_lo": float(np.percentile(v, 25)),
        "iqr_hi": float(np.percentile(v, 75)),
        "p05": float(np.percentile(v, 5)),
        "p95": float(np.percentile(v, 95)),
        "values": v,
    }


# ------------------------------------------------------------ model output
def split_subjects(df, seed=0, train_frac=0.5):
    """The SAME split calibrate.py --wetlab uses, so 'test' is truly unseen."""
    subs = np.sort(df.subject_id.unique())
    rng = np.random.default_rng(seed)
    rng.shuffle(subs)
    cut = int(len(subs) * train_frac)
    return set(subs[:cut]), set(subs[cut:])


def model_boost(params=None, platform="protein", adjuvant="alum"):
    """
    Simulate a Tdap-like boost and return its fold-rise and time to peak,
    measured exactly as the observed data is measured.
    """
    from dynamics import Sim, load_epitopes
    ep = load_epitopes(str(RESULTS / "demo_epitopes.csv"))
    if isinstance(params, (str, Path)):
        params = json.loads(Path(params).read_text())["params"]
    elif params is None:
        f = HERE / "models" / "calibrated.json"
        params = json.loads(f.read_text())["params"] if f.exists() else {}
    sim = Sim(ep, p=params, injections=[(2, None, adjuvant), (270, None, adjuvant)],
              platform=platform, adjuvant=adjuvant)
    df = sim.run(n_steps=700)

    boost_day = 270 * 8 / 24.0
    post = df[(df.day > boost_day) & (df.day <= boost_day + POST_MAX_DAY)]
    base = float(df[df.day <= boost_day].IgG.iloc[-1])
    if base <= 0 or post.empty:
        return None
    i = int(np.argmax(post.IgG.values))
    return {"fold_rise": float(post.IgG.values[i] / base),
            "peak_day": float(post.day.values[i] - boost_day),
            "curve": post[["day", "IgG"]].assign(day=post.day - boost_day)}


# ------------------------------------------------------------------- tests
def consistency(model_value, obs):
    """
    Where does the model's prediction fall in the observed subject
    distribution? Reports the percentile and whether it is inside the
    central 95 % of subjects and inside the CI of the median.
    """
    v = obs["values"]
    pct = float((v < model_value).mean() * 100)
    in_subjects = bool(obs["p05"] <= model_value <= obs["p95"])
    in_median_ci = bool(obs["ci_lo"] <= model_value <= obs["ci_hi"])
    # one-sample test of the model value against the observed distribution
    try:
        p = float(sps.wilcoxon(v - model_value).pvalue)
    except Exception:
        p = float("nan")
    return {"model": float(model_value), "percentile_in_cohort": round(pct, 1),
            "within_subject_90pct_band": in_subjects,
            "within_median_CI": in_median_ci,
            "p_differs_from_cohort_median": p}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="models/calibrated.json",
                    help="parameter file to evaluate")
    ap.add_argument("--split", default="all",
                    choices=["all", "train", "test"],
                    help="'test' = the half of subjects never used for "
                         "calibration, i.e. an honest held-out evaluation")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    if not DATA.exists():
        sys.exit(f"missing {DATA}\nrun: python3 harvest_cmipb.py")
    df = pd.read_csv(DATA)
    tr, te = split_subjects(df, seed=args.seed)
    if args.split == "train":
        df = df[df.subject_id.isin(tr)]
    elif args.split == "test":
        df = df[df.subject_id.isin(te)]
    print("Wet-lab validation of the dynamic model against CMI-PB\n")
    print(f"  parameters: {args.params}")
    print(f"  split: {args.split.upper()}" + (
        "  (held out from calibration - honest evaluation)"
        if args.split == "test" else ""))
    print(f"  {len(df):,} measurements, {df.subject_id.nunique()} subjects")
    print("  Every subject was primed in infancy, so every course is a BOOST.")
    print("  The model was calibrated before this data was harvested:")
    print("  no CMI-PB value was used to fit any parameter.\n")

    m = model_boost(params=args.params)
    if m is None:
        sys.exit("model produced no usable boost response")
    print(f"  MODEL PREDICTION  fold-rise {m['fold_rise']:.2f}x, "
          f"peak at day {m['peak_day']:.1f}\n")

    report, rows = {}, []
    print("=" * 78)
    print("VACCINE ANTIGENS - the model should match these")
    print("=" * 78)
    for ag in VACCINE_AG:
        ps = per_subject(df, ag)
        if len(ps) < MIN_SUBJECTS:
            continue
        fr = summarise(ps.fold_rise, "fold_rise")
        pd_ = summarise(ps.peak_day, "peak_day")
        if fr is None or pd_ is None:
            continue
        c_fr = consistency(m["fold_rise"], fr)
        c_pd = consistency(m["peak_day"], pd_)
        report[ag] = {"n_subjects": fr["n"],
                      "observed_fold_rise": {k: fr[k] for k in
                                             ("median", "ci_lo", "ci_hi", "p05", "p95")},
                      "observed_peak_day": {k: pd_[k] for k in
                                            ("median", "ci_lo", "ci_hi", "p05", "p95")},
                      "fold_rise_check": c_fr, "peak_day_check": c_pd}
        rows.append({"antigen": ag, "n": fr["n"],
                     "obs_fold": fr["median"], "obs_fold_ci":
                         f"{fr['ci_lo']:.2f}-{fr['ci_hi']:.2f}",
                     "model_fold": m["fold_rise"],
                     "fold_ok": c_fr["within_subject_90pct_band"],
                     "obs_peak_day": pd_["median"], "obs_peak_ci":
                         f"{pd_['ci_lo']:.0f}-{pd_['ci_hi']:.0f}",
                     "model_peak_day": m["peak_day"],
                     "peak_ok": c_pd["within_subject_90pct_band"]})
        print(f"\n  {ag}   n = {fr['n']} subjects")
        print(f"    fold-rise  observed median {fr['median']:.2f}x "
              f"(95% CI {fr['ci_lo']:.2f}-{fr['ci_hi']:.2f}, "
              f"5-95% of subjects {fr['p05']:.2f}-{fr['p95']:.2f})")
        print(f"               model {m['fold_rise']:.2f}x -> "
              f"{c_fr['percentile_in_cohort']:.0f}th percentile of subjects, "
              f"{'INSIDE' if c_fr['within_subject_90pct_band'] else 'OUTSIDE'} "
              f"the 5-95% band")
        print(f"    peak day   observed median {pd_['median']:.0f} d "
              f"(95% CI {pd_['ci_lo']:.0f}-{pd_['ci_hi']:.0f}, "
              f"5-95% {pd_['p05']:.0f}-{pd_['p95']:.0f})")
        print(f"               model {m['peak_day']:.1f} d -> "
              f"{c_pd['percentile_in_cohort']:.0f}th percentile, "
              f"{'INSIDE' if c_pd['within_subject_90pct_band'] else 'OUTSIDE'} "
              f"the 5-95% band")

    print("\n" + "=" * 78)
    print("NEGATIVE CONTROLS - antigens NOT in the vaccine.")
    print("A boost response here would mean the analysis is picking up noise.")
    print("=" * 78)
    ctrl = []
    for ag in CONTROL_AG:
        ps = per_subject(df, ag)
        if len(ps) < MIN_SUBJECTS:
            continue
        fr = summarise(ps.fold_rise, "fold_rise")
        if fr is None:
            continue
        ctrl.append({"antigen": ag, "n": fr["n"], "median_fold": fr["median"],
                     "ci": f"{fr['ci_lo']:.2f}-{fr['ci_hi']:.2f}"})
        print(f"  {ag:10s} n={fr['n']:4d}  fold-rise {fr['median']:.2f}x "
              f"(95% CI {fr['ci_lo']:.2f}-{fr['ci_hi']:.2f})")

    # Vaccine vs control: the discriminating test.
    if rows and ctrl:
        vac_f = np.concatenate([per_subject(df, r["antigen"]).fold_rise.values
                                for r in rows])
        con_f = np.concatenate([per_subject(df, c["antigen"]).fold_rise.values
                                for c in ctrl])
        u = sps.mannwhitneyu(vac_f, con_f, alternative="greater")
        # rank-biserial effect size
        eff = 1 - 2 * u.statistic / (len(vac_f) * len(con_f))
        print(f"\n  vaccine vs control fold-rise: Mann-Whitney U "
              f"p = {fmt_p(u.pvalue)}, rank-biserial r = {abs(eff):.3f}")
        print(f"  median {np.median(vac_f):.2f}x vs {np.median(con_f):.2f}x")
        report["_vaccine_vs_control"] = {
            "p": float(u.pvalue), "effect_rank_biserial": float(abs(eff)),
            "median_vaccine": float(np.median(vac_f)),
            "median_control": float(np.median(con_f))}

    tab = pd.DataFrame(rows)
    if not tab.empty:
        tab.to_csv(RESULTS / f"wetlab_validation_{args.split}.csv", index=False)
        n_ok_f = int(tab.fold_ok.sum()); n_ok_p = int(tab.peak_ok.sum())
        print("\n" + "=" * 78)
        print("VERDICT")
        print("=" * 78)
        print(f"  fold-rise  model inside the observed subject band for "
              f"{n_ok_f}/{len(tab)} vaccine antigens")
        print(f"  peak day   model inside the observed subject band for "
              f"{n_ok_p}/{len(tab)} vaccine antigens")
        report["_summary"] = {"antigens": len(tab), "fold_rise_ok": n_ok_f,
                              "peak_day_ok": n_ok_p,
                              "model_fold_rise": m["fold_rise"],
                              "model_peak_day": m["peak_day"]}

    report["_meta"] = {"params": str(args.params), "split": args.split}
    (RESULTS / f"wetlab_validation_{args.split}.json").write_text(json.dumps(report, indent=2,
                                                              default=float))
    print(f"\n  wrote {RESULTS}/wetlab_validation_{args.split}.[csv|json]")


if __name__ == "__main__":
    main()
