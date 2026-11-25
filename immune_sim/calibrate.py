#!/usr/bin/env python3
"""
Calibrate the dynamic model against quantitative immunological benchmarks.

Hand-tuning a 16-compartment ODE by eye does not converge: each fix displaces
the imbalance somewhere else. Instead we state what the response must do,
score any parameter set against those targets, and search.

TARGETS (textbook primary/secondary humoral response, and the qualitative
behaviour C-ImmSim is built to reproduce):

  1. primary IgG peaks  10-21 d after prime
  2. secondary IgG peaks 4-10 d after boost     (faster than primary)
  3. secondary/primary IgG peak ratio 5-50x     (stronger than primary)
  4. IgM:IgG ratio falls from primary to secondary (class switching)
  5. effector T cells contract >=80% from peak by 60 d after prime
  6. affinity matures 2-10x
  7. peaks stay in plausible ranges (no runaway, no collapse)

Score is the sum of squared log-distances outside each target band, so any
parameter set landing inside every band scores 0.

  python3 calibrate.py --iters 400
  python3 calibrate.py --iters 400 --out models/calibrated.json
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from dynamics import P, Sim, load_epitopes

warnings.filterwarnings("ignore")
HERE = Path(__file__).parent
MODELS = HERE / "models"
RESULTS = HERE / "results"

PRIME_STEP, BOOST_STEP, N_STEPS = 2, 270, 700
# CMI-PB is a Tdap boost: recombinant protein antigens on an aluminium
# adjuvant. Calibrate the model in that configuration so the comparison with
# the measured data is like for like.
PLATFORM, ADJUVANT = "protein", "alum"
PRIME_DAY = PRIME_STEP * 8 / 24
BOOST_DAY = BOOST_STEP * 8 / 24

def wetlab_targets(seed=0, train_frac=0.5):
    """
    Derive calibration bands from MEASURED human data (CMI-PB), using only a
    training half of the subjects. The held-out half is never seen here, so
    wetlab_validation.py --test-only remains an honest evaluation.

    Returns (targets, info) or (None, reason) if the data is not present.
    """
    import pandas as pd
    f = Path(__file__).parent / "data" / "cmipb_ab_titer.csv"
    if not f.exists():
        return None, f"missing {f} - run harvest_cmipb.py"
    sys.path.insert(0, str(Path(__file__).parent))
    from wetlab_validation import VACCINE_AG, per_subject

    df = pd.read_csv(f)
    subs = np.sort(df.subject_id.unique())
    rng = np.random.default_rng(seed)
    rng.shuffle(subs)
    train = set(subs[:int(len(subs) * train_frac)])
    d = df[df.subject_id.isin(train)]

    folds, peaks = [], []
    for ag in VACCINE_AG:
        ps = per_subject(d, ag)
        if len(ps) >= 15:
            folds.append(ps.fold_rise.values)
            peaks.append(ps.peak_day.values)
    if not folds:
        return None, "not enough training subjects"
    folds = np.concatenate(folds); peaks = np.concatenate(peaks)

    # Target the interquartile range of real subjects, not the extremes.
    f_lo, f_hi = np.percentile(folds, [25, 75])

    # Peak day needs care: CMI-PB samples at discrete days (0, 1, 3, 7, 14,
    # 30), so nearly every observed peak lands exactly on a sampling day and
    # the IQR collapses to a single point. A true peak anywhere between the
    # neighbouring sampling days would have been observed at the same day, so
    # the band is widened to that sampling resolution.
    sched = np.array([0, 1, 3, 7, 14, 30], float)
    obs_med = float(np.median(peaks))
    i = int(np.argmin(np.abs(sched - obs_med)))
    p_lo = float((sched[i - 1] + sched[i]) / 2) if i > 0 else float(sched[i])
    p_hi = (float((sched[i] + sched[i + 1]) / 2) if i + 1 < len(sched)
            else float(sched[i]))
    info = {"n_train_subjects": len(train), "n_courses": int(len(folds)),
            "fold_rise_IQR": [float(f_lo), float(f_hi)],
            "peak_day_band": [float(p_lo), float(p_hi)],
            "peak_day_band_note": "widened to the sampling resolution",
            "fold_rise_median": float(np.median(folds)),
            "peak_day_median": float(np.median(peaks))}
    return [("boost_fold_over_baseline", float(f_lo), float(f_hi)),
            ("secondary_peak_day", float(p_lo), float(p_hi))], info


# (name, low, high) - bands, not point targets
TARGETS = [
    ("primary_peak_day",   10.0,   21.0),
    ("secondary_peak_day",  4.0,   10.0),
    ("boost_ratio",         5.0,   50.0),
    ("switch_ratio",        2.0,  200.0),   # (IgG/IgM)_2nd / (IgG/IgM)_1st
    ("t_contraction",       0.80,   0.995),
    ("affinity_fold",       2.0,   10.0),
    ("peak_IgG",            1e3,    1e6),
    ("peak_TH",             1e4,    5e6),
    ("peak_CTL",            1e3,    2e6),
]

# parameters the search is allowed to move, with bounds (log-uniform)
TARGETS_ACTIVE = None      # set by main(); defaults to TARGETS

SEARCH = {
    "ag_deposit":   (0.005, 0.30),
    "agd_decay":    (0.03,  0.45),
    "dc_mature":    (0.005, 0.10),
    "dc_decay":     (0.05,  0.40),
    "th_prime":     (2e-4,  6e-3),
    "th_prolif":    (0.6,   2.2),
    "th_death":     (0.25,  0.95),
    "th_memory":    (0.01,  0.15),
    "th_mem_decay": (0.002, 0.05),
    "th_mem_recall":(3.0,   40.0),
    "tc_prolif":    (0.6,   2.2),
    "tc_death":     (0.25,  0.95),
    "b_activate":   (1e-3,  0.10),
    "b_prolif":     (0.6,   2.4),
    "b_death":      (0.10,  0.60),
    "b_memory":     (0.01,  0.20),
    "b_mem_decay":  (0.002, 0.05),
    "b_mem_recall": (3.0,   60.0),
    "plasma_diff":  (0.08,  0.60),
    "plasma_death": (0.15,  0.80),
    "igm_rate":     (1.0,   20.0),
    "igg_rate":     (3.0,   60.0),
    "switch_rate":  (0.05,  0.60),
    "aff_rate":     (0.02,  0.40),
    "k_ag_gc":      (5.0,   300.0),
    "k_switch":     (100.0, 5000.0),
    "k_aff_help":   (100.0, 5000.0),
    "k_carry":      (1e6,   5e7),
}


def measure(df):
    """Extract the benchmark quantities from one simulated time course."""
    m = {}
    pre = df[(df.day >= PRIME_DAY) & (df.day < BOOST_DAY)]
    post = df[df.day >= BOOST_DAY]
    if pre.empty or post.empty:
        return None

    i1 = int(np.argmax(pre.IgG.values))
    p1 = float(pre.IgG.values[i1])
    m["primary_peak_day"] = float(pre.day.values[i1] - PRIME_DAY)

    base = float(pre.IgG.values[-1])          # residual IgG at boost time
    i2 = int(np.argmax(post.IgG.values))
    p2 = float(post.IgG.values[i2])
    m["secondary_peak_day"] = float(post.day.values[i2] - BOOST_DAY)

    # rise above the residual, so a slowly-decaying primary cannot flatter it
    m["boost_ratio"] = (p2 - base) / p1 if p1 > 1e-9 else 0.0
    # The quantity wet-lab studies actually report: peak after the boost
    # divided by the pre-boost baseline in the same subject. This is NOT the
    # same as peak2/peak1 above, and conflating them overstates the response
    # by an order of magnitude.
    m["boost_fold_over_baseline"] = p2 / base if base > 1e-9 else 0.0
    m["peak_IgG"] = max(p1, p2)

    igm1 = float(pre.IgM.max()); igm2 = float(post.IgM.max())
    r1 = p1 / igm1 if igm1 > 1e-9 else 0.0
    r2 = p2 / igm2 if igm2 > 1e-9 else 0.0
    m["switch_ratio"] = r2 / r1 if r1 > 1e-9 else 0.0

    th_pk = float(pre.TH_active.max())
    w = pre[pre.day <= PRIME_DAY + 60]
    th_end = float(w.TH_active.values[-1]) if len(w) else th_pk
    m["t_contraction"] = 1.0 - th_end / th_pk if th_pk > 1e-9 else 0.0
    m["peak_TH"] = th_pk
    m["peak_CTL"] = float(pre.TC_active.max())
    m["affinity_fold"] = float(df.affinity.max())
    return m


def score(m, targets=None):
    if m is None:
        return 1e9, {}
    tot, detail = 0.0, {}
    for name, lo, hi in (targets if targets is not None else TARGETS):
        v = m.get(name, 0.0)
        if not np.isfinite(v) or v <= 0:
            pen = 25.0
        elif v < lo:
            pen = (np.log(lo / v)) ** 2
        elif v > hi:
            pen = (np.log(v / hi)) ** 2
        else:
            pen = 0.0
        detail[name] = (round(float(v), 4), round(float(pen), 3))
        tot += pen
    return tot, detail


# Cross-configuration constraints. Fitting a single configuration produced
# parameters that were badly wrong everywhere else - an unadjuvanted vaccine
# came out 800x STRONGER than an adjuvanted one, because nothing in the target
# set ever simulated it. These run alongside the primary fit.
CROSS_TARGETS = [
    ("adjuvant_effect",        5.0,   1000.0),   # alum IgG / no-adjuvant IgG
    ("boost_over_prime_peak",  1.2,     20.0),   # secondary > primary peak
    ("mrna_ctl_over_protein",  5.0,    1e4),     # endogenous routing -> CD8
]


def _run(params, ep, platform, adjuvant):
    sim = Sim(ep, p=params,
              injections=[(PRIME_STEP, None, adjuvant),
                          (BOOST_STEP, None, adjuvant)],
              platform=platform, adjuvant=adjuvant)
    return sim.run(n_steps=N_STEPS)


def evaluate(params, ep):
    try:
        df = _run(params, ep, PLATFORM, ADJUVANT)
        if not np.all(np.isfinite(df[["IgG", "IgM", "TH_active"]].values)):
            return 1e9, {}, None
        m = measure(df)
        if m is None:
            return 1e9, {}, None

        # other configurations the parameters must also get right
        none_df = _run(params, ep, "protein", "none")
        mrna_df = _run(params, ep, "mrna_lnp", "intrinsic")
        for d in (none_df, mrna_df):
            if not np.all(np.isfinite(d[["IgG", "TC_active"]].values)):
                return 1e9, {}, None

        m["adjuvant_effect"] = (float(df.IgG.max())
                                / max(float(none_df.IgG.max()), 1e-9))
        bd = BOOST_DAY
        pre, post = df[df.day < bd], df[df.day >= bd]
        m["boost_over_prime_peak"] = (float(post.IgG.max())
                                      / max(float(pre.IgG.max()), 1e-9))
        m["mrna_ctl_over_protein"] = (float(mrna_df.TC_active.max())
                                      / max(float(df.TC_active.max()), 1e-9))

        s, d = score(m, list(TARGETS_ACTIVE) + CROSS_TARGETS)
        return s, d, df
    except Exception:
        return 1e9, {}, None


def search(ep, iters, seed=0, out=None):
    rng = np.random.default_rng(seed)
    best_p = dict(P)
    best_s, best_d, _ = evaluate(best_p, ep)
    print(f"  baseline score {best_s:.3f}")

    cur_p, cur_s = dict(best_p), best_s
    for it in range(1, iters + 1):
        cand = dict(cur_p)
        # perturb a random subset, log-uniform within bounds
        for k in rng.choice(list(SEARCH), size=rng.integers(2, 7), replace=False):
            lo, hi = SEARCH[k]
            if rng.random() < 0.3:                       # occasional big jump
                cand[k] = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            else:                                        # local step
                v = cur_p[k] * float(np.exp(rng.normal(0, 0.35)))
                cand[k] = float(np.clip(v, lo, hi))

        s, d, _ = evaluate(cand, ep)
        if s < cur_s or rng.random() < 0.05:             # mild annealing
            cur_p, cur_s = cand, s
        if s < best_s:
            best_p, best_s, best_d = cand, s, d
            print(f"  iter {it:4d}  score {s:.4f}", flush=True)
            if out:                      # checkpoint every improvement, so a
                save(best_p, best_s, out)  # killed run still leaves its best
            if s == 0.0:
                print("  all targets met")
                break
        if it % 50 == 0:
            print(f"  iter {it:4d}  (best {best_s:.4f})", flush=True)
    return best_p, best_s, best_d


def save(params, sc, out):
    changed = {k: round(v, 6) for k, v in params.items()
               if abs(v - P.get(k, v)) > 1e-12}
    out = Path(out); out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps({"score": sc, "params": params,
                               "changed": changed}, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--epitopes", default="results/demo_epitopes.csv")
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="models/calibrated.json")
    ap.add_argument("--wetlab", action="store_true",
                    help="replace the literature bands for boost magnitude and "
                         "timing with bands measured from CMI-PB TRAINING "
                         "subjects (the held-out half stays unseen)")
    args = ap.parse_args()

    global TARGETS_ACTIVE
    TARGETS_ACTIVE = list(TARGETS)
    if args.wetlab:
        wl, info = wetlab_targets()
        if wl is None:
            sys.exit(f"--wetlab: {info}")
        names = {n for n, _, _ in wl}
        TARGETS_ACTIVE = [t for t in TARGETS_ACTIVE if t[0] not in names] + wl
        # peak2/peak1 is superseded by the measured fold-over-baseline
        TARGETS_ACTIVE = [t for t in TARGETS_ACTIVE if t[0] != "boost_ratio"]
        print("Calibrating against MEASURED human data (CMI-PB training half)")
        print(f"  {info['n_train_subjects']} training subjects, "
              f"{info['n_courses']} antigen time courses")
        print(f"  fold-rise  median {info['fold_rise_median']:.2f}x, "
              f"IQR {info['fold_rise_IQR'][0]:.2f}-{info['fold_rise_IQR'][1]:.2f}")
        print(f"  peak day   median {info['peak_day_median']:.0f} d, "
              f"band {info['peak_day_band'][0]:.1f}-{info['peak_day_band'][1]:.1f} d "
              f"({info['peak_day_band_note']})\n")

    ep = load_epitopes(args.epitopes)
    print(f"Calibrating on {len(ep)} epitope clones, {args.iters} iterations\n")
    bp, bs, bd = search(ep, args.iters, args.seed, out=args.out)

    print(f"\n  final score {bs:.4f}")
    print(f"\n  {'target':26s} {'value':>12s} {'band':>22s} {'penalty':>8s}")
    for name, lo, hi in list(TARGETS_ACTIVE) + CROSS_TARGETS:
        v, pen = bd.get(name, (float('nan'), float('nan')))
        ok = "ok" if pen == 0 else "  "
        print(f"  {name:26s} {v:12.4g} {f'[{lo:g}, {hi:g}]':>22s} {pen:8.3f} {ok}")

    save(bp, bs, args.out)
    changed = sum(1 for k, v in bp.items() if abs(v - P.get(k, v)) > 1e-12)
    print(f"\n  changed {changed} parameters -> {args.out}")


if __name__ == "__main__":
    main()
