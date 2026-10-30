"""Outputs: ranked CSV, FASTA, per-cycle ensemble record, comparison table and plots."""
from __future__ import annotations

import csv
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .config import Params
from .score import Scored


def _flat(d: dict, prefix="") -> dict:
    out = {}
    for k, v in d.items():
        if isinstance(v, dict):
            out.update(_flat(v, f"{prefix}{k}_"))
        elif isinstance(v, (list, tuple)):
            continue
        else:
            out[f"{prefix}{k}"] = v
    return out


def write_outputs(result: dict, p: Params, outdir: str, log=print) -> dict:
    os.makedirs(outdir, exist_ok=True)
    plotdir = os.path.join(outdir, "plots")
    os.makedirs(plotdir, exist_ok=True)
    paths: dict[str, str] = {}
    finalists: list[Scored] = result["finalists"]
    original: Scored = result["original"]
    extras: list[dict] = result["extras"]

    # ---- ranked table -------------------------------------------------------
    rows = []
    orow = original.row(p)
    orow.update({"rank": 0, "id": "original", "n_codon_changes": 0,
                 "pct_codons_changed": 0.0})
    rows.append(orow)
    base = original.genome
    for i, (s, ex) in enumerate(zip(finalists, extras), start=1):
        r = s.row(p)
        diff = sum(1 for a, b in zip(s.genome, base) if a != b)
        r.update({"rank": i, "id": f"variant_{i:02d}", "n_codon_changes": diff,
                  "pct_codons_changed": 100.0 * diff / max(len(base), 1)})
        r.update(_flat(ex))
        rows.append(r)
    cols = ["rank", "id", "total_score", "feasible", "violations"]
    rest = sorted({k for r in rows for k in r} - set(cols))
    paths["ranked_csv"] = os.path.join(outdir, "ranked_variants.csv")
    with open(paths["ranked_csv"], "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols + rest, extrasaction="ignore")
        wr.writeheader()
        for r in rows:
            wr.writerow(r)

    # ---- FASTA --------------------------------------------------------------
    paths["fasta"] = os.path.join(outdir, "variants.fasta")
    with open(paths["fasta"], "w") as fh:
        fh.write(f">{p.name}_original len={len(original.seq)} circular\n")
        for i in range(0, len(original.seq), 70):
            fh.write(original.seq[i:i + 70] + "\n")
        for i, s in enumerate(finalists, start=1):
            w = min(s.cai, key=s.cai.get)
            fh.write(f">{p.name}_variant_{i:02d} score={s.total:.4f} "
                     f"cai_gm={s.cai_gm:.3f} cai_weak({w})={s.cai[w]:.3f} "
                     f"mfe_per_nt={s.metrics.mfe_per_nt:.3f} dsrna_nt={s.metrics.dsrna_nt} "
                     f"feasible={s.feasible} circular\n")
            for j in range(0, len(s.seq), 70):
                fh.write(s.seq[j:j + 70] + "\n")

    # ---- best structure + per-cycle ensembles -------------------------------
    best: Scored = result["best"]
    paths["best_structure"] = os.path.join(outdir, "best_structure.txt")
    with open(paths["best_structure"], "w") as fh:
        fh.write(f"# {p.name} best variant (circular, {p.temperature} C, dangles={p.dangles})\n")
        fh.write(best.seq + "\n")
        fh.write(best.metrics.structure + f" ({best.metrics.mfe:.2f})\n")
        if best.metrics.centroid:
            fh.write(best.metrics.centroid + f" centroid (d={best.metrics.centroid_distance:.2f})\n")
    paths["ensembles"] = os.path.join(outdir, "ensemble_by_cycle.txt")
    with open(paths["ensembles"], "w") as fh:
        fh.write("# Boltzmann-sampled circular structures of each cycle's best variant.\n")
        fh.write("# This is the ensemble scan across optimisation cycles.\n")
        for rl in result["rounds"]:
            fh.write(f"\n## cycle {rl.round}  total={rl.best_total:.4f}  "
                     f"stats={json.dumps(rl.ensemble_stats)}\n")
            for s in rl.ensemble:
                fh.write(s + "\n")
        fh.write(f"\n## FINAL best  stats={json.dumps(result['best_ensemble_stats'])}\n")
        for s in result["best_ensemble"]:
            fh.write(s + "\n")

    # ---- comparison ---------------------------------------------------------
    paths["comparison_csv"] = os.path.join(outdir, "comparison_original_vs_best.csv")
    keys = ["total_score", "feasible", "cai_gm", f"cai_{p.primary_species}",
            f"cai_{p.secondary_species}", "gc", "gc_win_min", "gc_win_max", "mfe", "mfe_per_nt",
            "paired_fraction", "dsrna_nt", "dsrna_helices", "longest_stem", "long_range_pairs",
            "start_paired", "start_unpaired_prob", "early_paired_fraction", "efe",
            "ensemble_diversity", "mean_entropy", "red_pairs", "red_longest_run",
            "red_cds_crossing", "centroid_dist_to_mfe", "mfe_frequency", "motif_violations"]
    orow, brow = original.row(p), best.row(p)
    with open(paths["comparison_csv"], "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["metric", "original", "best", "delta"])
        for k in keys:
            a, b = orow.get(k), brow.get(k)
            try:
                dl = round(float(b) - float(a), 6)
            except (TypeError, ValueError):
                dl = ""
            wr.writerow([k, a, b, dl])

    # ---- run record ---------------------------------------------------------
    paths["run_json"] = os.path.join(outdir, "run.json")
    with open(paths["run_json"], "w") as fh:
        json.dump({
            "params": p.to_dict(), "design": result["design"],
            "elapsed_s": result["elapsed_s"], "weights_final": result["weights_final"],
            "rounds": [{"round": r.round, "best_total": r.best_total,
                        "elapsed_s": r.elapsed_s, "weights": r.weights,
                        "ensemble_stats": r.ensemble_stats, "notes": r.notes,
                        "best": _flat(r.best_row)} for r in result["rounds"]],
            "pareto": [{"total": s.total, "cai": s.cai, "mfe_per_nt": s.metrics.mfe_per_nt,
                        "dsrna_nt": s.metrics.dsrna_nt, "seq": s.seq} for s in result["pareto"]],
            "best": {"seq": best.seq, "structure": best.metrics.structure,
                     "row": _flat(brow)},
            "original": {"seq": original.seq, "row": _flat(orow)},
            "derna": result.get("derna", []),
        }, fh, indent=2, default=str)

    # ---- plots --------------------------------------------------------------
    paths["plots"] = plotdir
    try:
        cds = tuple(result["design"]["cds"])
        _plot_variant(best, p, os.path.join(plotdir, "best"), "best variant", cds)
        _plot_variant(original, p, os.path.join(plotdir, "original"), "original", cds)
        _plot_convergence(result, os.path.join(plotdir, "convergence.png"))
        _plot_pareto(result, os.path.join(plotdir, "pareto.png"))
    except Exception as e:                                  # plots must never kill a run
        log(f"warning: plotting failed: {type(e).__name__}: {e}")
    return paths


def _plot_variant(s: Scored, p: Params, stem: str, title: str, cds: tuple[int, int]) -> None:
    m = s.metrics
    if not m or not m.full:
        return
    n = len(s.seq)
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), height_ratios=[2, 1, 1])
    # base-pair probability map
    if m.bp_prob_pairs:
        xs = [i for i, _, _ in m.bp_prob_pairs]
        ys = [j for _, j, _ in m.bp_prob_pairs]
        cs = [pr for _, _, pr in m.bp_prob_pairs]
        sc = ax[0].scatter(xs, ys, c=cs, s=4, cmap="viridis", vmin=0, vmax=1)
        fig.colorbar(sc, ax=ax[0], label="pair probability")
    for a in ax:                       # CDS bounds on every panel (start codon in red)
        a.axvline(cds[0], color="#d62728", lw=0.9, ls="--")
        a.axvline(cds[1], color="#7f7f7f", lw=0.9, ls="--")
    ax[0].set_xlim(0, n)
    ax[0].set_ylim(0, n)
    ax[0].set_title(f"{title}: base-pair probabilities (circular fold)")
    ax[0].set_xlabel("i")
    ax[0].set_ylabel("j")
    # mountain plot
    ax[1].plot(m.mountain, lw=0.8)
    ax[1].set_ylabel("enclosing pairs")
    ax[1].set_title("mountain plot (MFE structure)")
    ax[1].set_xlim(0, n)
    # entropy + unpaired probability
    ax[2].plot(m.entropy, lw=0.7, label="positional entropy")
    ax[2].plot(m.unpaired_prob, lw=0.7, alpha=0.7, label="p(unpaired)")
    ax[2].set_xlim(0, n)
    ax[2].set_xlabel("position (nt)")
    ax[2].legend(loc="upper right", fontsize=8)
    ax[2].set_title("confidence tracks")
    fig.tight_layout()
    fig.savefig(stem + "_structure.png", dpi=130)
    plt.close(fig)


def _plot_convergence(result: dict, path: str) -> None:
    rounds = result["rounds"]
    if not rounds:
        return
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    x = [r.round for r in rounds]
    ax[0].plot(x, [r.best_total for r in rounds], "o-")
    ax[0].set_xlabel("cycle")
    ax[0].set_title("best composite score")
    # MFE/nt and dsRNA nt have incompatible scales -> twin axes
    ax[1].plot(x, [r.best_row.get("mfe_per_nt", 0) for r in rounds], "o-", color="#ff7f0e",
               label="MFE/nt (higher = less paired)")
    ax[1].set_ylabel("MFE per nt", color="#ff7f0e")
    ax[1].tick_params(axis="y", labelcolor="#ff7f0e")
    tw = ax[1].twinx()
    tw.plot(x, [r.best_row.get("dsrna_nt", 0) for r in rounds], "s--", color="#1f77b4",
            label="dsRNA nt")
    tw.set_ylabel("nt in long duplex", color="#1f77b4")
    tw.tick_params(axis="y", labelcolor="#1f77b4")
    tw.set_ylim(bottom=0)
    ax[1].set_xlabel("cycle")
    ax[1].set_title("structure of the cycle best")
    for k in [k for k in rounds[0].best_row if k.startswith("cai_")]:
        ax[2].plot(x, [r.best_row.get(k) for r in rounds], "o-", label=k)
    ax[2].legend(fontsize=8)
    ax[2].set_xlabel("cycle")
    ax[2].set_title("codon adaptation per species")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def _plot_pareto(result: dict, path: str) -> None:
    pts = result["pareto"]
    if not pts:
        return
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    x = [min(s.cai.values()) for s in pts]
    y = [s.metrics.mfe_per_nt for s in pts]
    c = [s.metrics.dsrna_nt for s in pts]
    sc = ax.scatter(x, y, c=c, cmap="magma", s=28)
    fig.colorbar(sc, ax=ax, label="dsRNA nt (>= threshold duplex)")
    b = result["best"]
    ax.scatter([min(b.cai.values())], [b.metrics.mfe_per_nt], marker="*", s=220,
               edgecolor="k", facecolor="none", label="selected")
    o = result["original"]
    ax.scatter([min(o.cai.values())], [o.metrics.mfe_per_nt], marker="s", s=60,
               edgecolor="r", facecolor="none", label="original")
    ax.set_xlabel("weaker-species CAI")
    ax.set_ylabel("MFE per nt (higher = less structured)")
    ax.legend(fontsize=8)
    ax.set_title("non-dominated set")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)
