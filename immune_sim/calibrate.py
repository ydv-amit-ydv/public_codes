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
PRIME_DAY = PRIME_STEP * 8 / 24
BOOST_DAY = BOOST_STEP * 8 / 24

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


def score(m):
    if m is None:
        return 1e9, {}
    tot, detail = 0.0, {}
    for name, lo, hi in TARGETS:
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


def evaluate(params, ep):
    try:
        sim = Sim(ep, p=params,
                  injections=[(PRIME_STEP, 1000.0, 100.0),
                              (BOOST_STEP, 1000.0, 100.0)])
        df = sim.run(n_steps=N_STEPS)
        if not np.all(np.isfinite(df[["IgG", "IgM", "TH_active"]].values)):
            return 1e9, {}, None
        m = measure(df)
        s, d = score(m)
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
    ap.add_argument("--epitopes", default="results/P0DTC2_epitopes.csv")
    ap.add_argument("--iters", type=int, default=300)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="models/calibrated.json")
    args = ap.parse_args()

    ep = load_epitopes(args.epitopes)
    print(f"Calibrating on {len(ep)} epitope clones, {args.iters} iterations\n")
    bp, bs, bd = search(ep, args.iters, args.seed, out=args.out)

    print(f"\n  final score {bs:.4f}")
    print(f"\n  {'target':22s} {'value':>12s} {'band':>22s} {'penalty':>8s}")
    for name, lo, hi in TARGETS:
        v, pen = bd.get(name, (float('nan'), float('nan')))
        ok = "ok" if pen == 0 else "  "
        print(f"  {name:22s} {v:12.4g} {f'[{lo:g}, {hi:g}]':>22s} {pen:8.3f} {ok}")

    save(bp, bs, args.out)
    changed = sum(1 for k, v in bp.items() if abs(v - P.get(k, v)) > 1e-12)
    print(f"\n  changed {changed} parameters -> {args.out}")


if __name__ == "__main__":
    main()
