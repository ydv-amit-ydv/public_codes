#!/usr/bin/env python3
"""
Run the local dynamic model on predicted epitopes and draw C-ImmSim-style plots.

  python3 simulate.py results/P0DTC2_epitopes.csv
  python3 simulate.py results/P0DTC2_epitopes.csv --inject 2,86 --steps 540
  python3 simulate.py results/P0DTC2_epitopes.csv --no-adjuvant   # LPS off
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from dynamics import Sim, load_epitopes
from platforms import ADJUVANTS, PLATFORMS, describe

RESULTS = Path(__file__).parent / "results"
RESULTS.mkdir(exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 130, "font.size": 9, "axes.grid": True,
    "grid.alpha": 0.25, "axes.spines.top": False, "axes.spines.right": False,
})

C = {"antigen": "#444444", "IgM": "#1f77b4", "IgG": "#d62728",
     "TH": "#2ca02c", "TC": "#9467bd", "B": "#ff7f0e",
     "mem": "#8c564b", "IL2": "#17becf", "danger": "#bcbd22"}


def plot_panel(df, ep, out, title, inj_days):
    fig, ax = plt.subplots(3, 2, figsize=(11, 10))

    def marks(a):
        for d in inj_days:
            a.axvline(d, color="k", ls=":", lw=0.9, alpha=0.6)

    # 1 antigen + antibody
    a = ax[0, 0]
    a.semilogy(df.day, np.clip(df.antigen, 1e-2, None), color=C["antigen"],
               label="antigen")
    a.semilogy(df.day, np.clip(df.IgM, 1e-2, None), color=C["IgM"], label="IgM")
    a.semilogy(df.day, np.clip(df.IgG, 1e-2, None), color=C["IgG"], label="IgG")
    a.set_title("Antigen and antibody"); a.set_ylabel("conc. (a.u.)")
    a.legend(fontsize=8); marks(a)

    # 2 T cells
    a = ax[0, 1]
    a.plot(df.day, df.TH_active, color=C["TH"], label="TH active")
    a.plot(df.day, df.TC_active, color=C["TC"], label="CTL active")
    a.plot(df.day, df.TH_memory, color=C["TH"], ls="--", label="TH memory")
    a.plot(df.day, df.TC_memory, color=C["TC"], ls="--", label="CTL memory")
    a.set_title("T cell response"); a.set_ylabel("cells")
    a.legend(fontsize=8); marks(a)

    # 3 B cells and plasma
    a = ax[1, 0]
    a.plot(df.day, df.B_active, color=C["B"], label="B active")
    a.plot(df.day, df.plasma, color=C["IgG"], label="plasma")
    a.plot(df.day, df.B_memory, color=C["mem"], ls="--", label="B memory")
    a.set_title("B cell compartment"); a.set_ylabel("cells")
    a.legend(fontsize=8); marks(a)

    # 4 APC and signals
    a = ax[1, 1]
    a.plot(df.day, df.DC_mature, color="#e377c2", label="mature DC")
    a.plot(df.day, df.MA_active, color="#7f7f7f", label="active MA")
    a.plot(df.day, df.danger, color=C["danger"], ls="--", label="danger signal")
    a2 = a.twinx(); a2.plot(df.day, df.IL2, color=C["IL2"], lw=1, label="IL-2")
    a2.set_ylabel("IL-2", color=C["IL2"]); a2.grid(False)
    a.set_title("Innate / presentation and cytokine"); a.set_ylabel("cells")
    a.legend(fontsize=8, loc="upper right"); marks(a)

    # 5 affinity maturation
    a = ax[2, 0]
    a.plot(df.day, df.affinity, color="#17becf")
    a.set_title("Antibody affinity maturation")
    a.set_ylabel("fold over naive"); a.set_xlabel("day"); marks(a)

    # 6 immunodominance: per-clone T cell expansion
    a = ax[2, 1]
    clone_cols = [c for c in df.columns if c.startswith(("TH_clone_", "TC_clone_"))]
    peak = {c: df[c].max() for c in clone_cols}
    for c in sorted(peak, key=peak.get, reverse=True)[:6]:
        lbl = c.replace("TH_clone_", "TH ").replace("TC_clone_", "CTL ")
        a.plot(df.day, df[c], lw=1.2, label=lbl)
    a.set_title("Immunodominance: top responding clones")
    a.set_ylabel("cells"); a.set_xlabel("day")
    a.legend(fontsize=6.5, ncol=1); marks(a)

    fig.suptitle(title, fontsize=11, y=0.995)
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)


def plot_prime_boost(df, inj_days, out):
    """Zoomed comparison of primary vs secondary response."""
    if len(inj_days) < 2:
        return
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
    for i, (d0, lab) in enumerate(zip(inj_days[:2], ["Primary", "Secondary"])):
        w = df[(df.day >= d0 - 2) & (df.day <= d0 + 40)]
        rel = w.day - d0
        ax[i].semilogy(rel, np.clip(w.IgM, 1e-2, None), color=C["IgM"], label="IgM")
        ax[i].semilogy(rel, np.clip(w.IgG, 1e-2, None), color=C["IgG"], label="IgG")
        pk = w.IgG.max()
        t_pk = rel.values[np.argmax(w.IgG.values)]
        ax[i].set_title(f"{lab}: peak IgG {pk:,.0f} at day {t_pk:.0f}")
        ax[i].set_xlabel("days after injection")
        ax[i].legend(fontsize=8)
    ax[0].set_ylabel("antibody (a.u.)")
    fig.tight_layout(); fig.savefig(out, bbox_inches="tight"); plt.close(fig)


def summarize(df, inj_days):
    out = {}
    for i, d0 in enumerate(inj_days, start=1):
        end = inj_days[i] if i < len(inj_days) else df.day.max()
        w = df[(df.day >= d0) & (df.day <= min(end, d0 + 60))]
        if w.empty:
            continue
        out[f"injection_{i}"] = {
            "day": round(float(d0), 1),
            "peak_IgG": round(float(w.IgG.max()), 1),
            "peak_IgG_day_rel": round(float(w.day.values[np.argmax(w.IgG.values)] - d0), 1),
            "peak_IgM": round(float(w.IgM.max()), 1),
            "peak_TH": round(float(w.TH_active.max()), 1),
            "peak_CTL": round(float(w.TC_active.max()), 1),
            "peak_plasma": round(float(w.plasma.max()), 1),
            "affinity_end": round(float(w.affinity.iloc[-1]), 2),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("epitopes", nargs="?", default=None,
                    help="results/<name>_epitopes.csv from predict.py")
    ap.add_argument("--engine", default=None,
                    help="epitope_filter_engine run dir or *_pairs.csv - "
                         "PREFERRED over our own model (11 fused predictors)")
    ap.add_argument("--platform", default="protein",
                    choices=sorted(PLATFORMS), help="vaccine platform")
    ap.add_argument("--adjuvant", default=None, choices=sorted(ADJUVANTS))
    ap.add_argument("--dose-ug", type=float, default=None)
    ap.add_argument("--heterologous", default=None, choices=sorted(PLATFORMS),
                    help="boost with a different platform (avoids anti-vector)")
    ap.add_argument("--list-platforms", action="store_true")
    ap.add_argument("--inject", default="2,86",
                    help="comma-separated injection steps (8h each)")
    ap.add_argument("--no-adjuvant", action="store_true",
                    help="LPS off - tests the danger-signal dependence")
    ap.add_argument("--steps", type=int, default=540)
    ap.add_argument("--replicating", action="store_true",
                    help="antigen replicates (live pathogen rather than vaccine)")
    ap.add_argument("--params", default="models/calibrated.json",
                    help="calibrated parameter file; 'none' for defaults")
    ap.add_argument("--tag", default=None)
    args = ap.parse_args()

    if args.list_platforms:
        for k in sorted(PLATFORMS):
            print(describe(k)); print()
        return
    if not args.epitopes and not args.engine:
        ap.error("give an epitope CSV or --engine <run dir>")
    steps = [int(s) for s in args.inject.split(",") if s.strip()]
    adj_key = "none" if args.no_adjuvant else args.adjuvant
    inj = [(st, args.dose_ug, adj_key) for st in steps]
    inj_days = [s * 8.0 / 24.0 for s in steps]

    if args.engine:
        import engine_adapter
        ep = engine_adapter.load(args.engine)
        name = Path(args.engine).stem.replace("_pairs", "").replace("_peptides", "")
        source = f"epitope_filter_engine ({args.engine})"
    else:
        ep = load_epitopes(args.epitopes)
        name = Path(args.epitopes).stem.replace("_epitopes", "")
        source = "immune_sim IEDB-trained model"
    tag = args.tag or (args.platform + ("_noadj" if args.no_adjuvant else ""))

    print(f"Epitope source: {source}")
    print(f"Clones       : {len(ep)}  "
          f"({(ep.mhc_class=='I').sum()} class I, {(ep.mhc_class=='II').sum()} class II)")
    print(f"Injections   : steps {steps} (days {[round(d,1) for d in inj_days]})")
    print()
    print(describe(args.platform, adj_key, args.dose_ug))
    print()
    cols = ["allele", "mhc_class", "peptide", "pred_IC50_nM", "immuno_p", "weight"]
    print(ep[cols].to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    p = {}
    if args.params and args.params.lower() != "none":
        pf = Path(args.params)
        if pf.exists():
            p.update(json.loads(pf.read_text())["params"])
            print(f"\nusing calibrated parameters from {pf}")
        else:
            print(f"\nWARNING: {pf} not found, using defaults")
    if args.replicating:
        p["ag_growth"] = 0.55
    sim = Sim(ep, p=p, injections=inj, platform=args.platform,
              adjuvant=adj_key, dose_ug=args.dose_ug,
              heterologous=args.heterologous)
    df = sim.run(n_steps=args.steps)

    csv = RESULTS / f"{name}_{tag}_timecourse.csv"
    df.to_csv(csv, index=False)

    title = (f"{name} - {PLATFORMS[args.platform].label} + "
             f"{sim.adj.label} - {len(ep)} epitope clones, steps {steps}")
    png = RESULTS / f"{name}_{tag}_response.png"
    plot_panel(df, ep, png, title, inj_days)
    pb = RESULTS / f"{name}_{tag}_primeboost.png"
    plot_prime_boost(df, inj_days, pb)

    s = summarize(df, inj_days)
    print("\n--- response summary ---")
    for k, v in s.items():
        print(f"{k}: peak IgG {v['peak_IgG']:>10,.0f} at +{v['peak_IgG_day_rel']:>4.0f} d | "
              f"IgM {v['peak_IgM']:>9,.0f} | TH {v['peak_TH']:>8,.0f} | "
              f"CTL {v['peak_CTL']:>8,.0f} | affinity x{v['affinity_end']:.2f}")
    if len(s) >= 2:
        a, b = list(s.values())[0], list(s.values())[1]
        if a["peak_IgG"] > 0:
            print(f"\nboost / prime IgG ratio : {b['peak_IgG']/a['peak_IgG']:.2f}x")
            print(f"time to peak  prime {a['peak_IgG_day_rel']:.0f} d "
                  f"-> boost {b['peak_IgG_day_rel']:.0f} d")

    print(f"\nwrote {csv}")
    print(f"wrote {png}")
    if pb.exists():
        print(f"wrote {pb}")


if __name__ == "__main__":
    main()
