"""Pick the parameters for the user, from the sequence and a time budget.

Two jobs:

1. *Scientific defaults.* Nearly every knob in this package is the same trade-off wearing a
   different hat -- codon adaptation against secondary structure. The defaults below encode
   the positions that have a mechanism behind them (PKR needs ~30 bp of contiguous duplex,
   initiation needs an accessible AUG) and leave the rest as soft pressure.

2. *Sizing the search to a budget.* Measured on this build, a fast candidate evaluation is
   O(n) at ~1.22 ms/nt and an exact circular partition function is O(n^2.7). Those two laws
   decide how much search fits in the time the user is willing to wait, so a 3000 nt design
   does not quietly turn into an overnight job.

Auto only fills values the user LEFT AT THE DEFAULT. Anything explicitly set is respected.
"""
from __future__ import annotations

import time
from dataclasses import fields

from .config import Params

# cost model, calibrated on the reference machine (ViennaRNA 2.7.2, one core)
COST_FAST_PER_NT = 1.22e-3        # seconds
COST_EXACT_COEF = 4.6e-8          # seconds, times n**COST_EXACT_EXP
COST_EXACT_EXP = 2.70

# (rounds, generations, population, elite_k, top_n, polish) from cheap to thorough
LADDER = [
    (2, 4, 12, 3, 2, 0),
    (3, 6, 16, 4, 3, 0),
    (3, 8, 24, 4, 3, 2),
    (4, 10, 32, 5, 4, 3),
    (4, 12, 40, 6, 6, 4),
    (5, 16, 48, 6, 6, 6),
    (6, 20, 64, 8, 8, 8),
    (8, 25, 96, 8, 8, 10),
    (10, 30, 128, 10, 10, 12),
]


def t_fast(n: int, factor: float = 1.0) -> float:
    return COST_FAST_PER_NT * n * factor


def t_exact(n: int, factor: float = 1.0) -> float:
    return COST_EXACT_COEF * (n ** COST_EXACT_EXP) * factor


def calibrate(seq: str, p: Params, cds) -> float:
    """Time one real fast evaluation and return this machine's speed factor (1.0 = reference)."""
    from .structure import fold
    t0 = time.time()
    fold(seq, p, cds, circ=True, full=False, local_start=True)
    measured = time.time() - t0
    predicted = t_fast(len(seq))
    if predicted <= 0:
        return 1.0
    return max(0.15, min(8.0, measured / predicted))


def estimate_seconds(cfg, n: int, workers: int, factor: float = 1.0) -> float:
    """Wall-clock estimate for one full run at configuration `cfg`."""
    rounds, gens, pop, elite, top_n, polish = cfg
    w = max(1, workers)
    ga = rounds * gens * pop * t_fast(n, factor) / w
    # exact circular partition functions: elites + one ensemble sample per cycle,
    # the final re-score of every unique elite, and the polish batches
    n_exact = rounds * (elite + 1) + rounds * elite + top_n * polish * max(4, w) * 0.5
    exact = n_exact * t_exact(n, factor) / min(w, max(1, elite + 1))
    return ga + exact


def choose_search(n: int, workers: int, budget_s: float, factor: float = 1.0):
    """Richest ladder rung that fits the budget (always at least the cheapest)."""
    best = LADDER[0]
    for cfg in LADDER:
        if estimate_seconds(cfg, n, workers, factor) <= budget_s:
            best = cfg
        else:
            break
    return best


def achievable_gc(design, cai_at: float) -> dict:
    """How low GC can go *while still translating well* -- the band must respect both.

    The unconstrained GC-poorest synonymous form is useless for setting a band: it is built
    from the rarest codons and its CAI collapses, so a band drawn around it can never be met
    at the CAI the run is actually aiming for. This instead starts from the codon-optimal
    sequence and gives GC away one codon at a time, stopping at `cai_at`. The result is
    reachable by construction, which is the whole point.
    """
    import math
    from .constraints import gc_stats
    model, p = design.model, design.p
    species = (p.primary_species, p.secondary_species)
    w = p.primary_weight

    def weak_cai(g):
        seq = design.build(g)
        return min(model.index(design.cds_codons(seq), sp) for sp in species)

    best = [max(range(len(o)), key=lambda i, o=o: w * math.log(model.w[species[0]][o[i]])
                + (1 - w) * math.log(model.w[species[1]][o[i]])) for o in design.options]
    hi_g = [max(range(len(o)), key=lambda i, o=o: sum(1 for b in o[i] if b in "GC"))
            for o in design.options]
    gc_hi = gc_stats(design.build(hi_g), p.gc_window, design.free_mask)["gc"]

    # give away GC in the order that costs the least codon fitness
    g = list(best)
    order = []
    for k, opts in enumerate(design.options):
        cur = opts[g[k]]
        for i, c in enumerate(opts):
            d_gc = sum(1 for b in cur if b in "GC") - sum(1 for b in c if b in "GC")
            if d_gc <= 0:
                continue
            d_cai = (w * (math.log(model.w[species[0]][cur]) - math.log(model.w[species[0]][c]))
                     + (1 - w) * (math.log(model.w[species[1]][cur]) - math.log(model.w[species[1]][c])))
            order.append((d_cai / d_gc, k, i))
    order.sort()
    for _, k, i in order:
        old = g[k]
        g[k] = i
        if weak_cai(g) < cai_at:
            g[k] = old
    st = gc_stats(design.build(g), p.gc_window, design.free_mask)
    base = gc_stats(design.build(best), p.gc_window, design.free_mask)
    return {"gc_lo_feasible": st["gc"], "win_lo_feasible": st["gc_win_max_free"],
            "gc_at_max_cai": base["gc"], "win_at_max_cai": base["gc_win_max_free"],
            "gc_hi": gc_hi, "cai_at": cai_at}


def _at_default(p: Params, name: str) -> bool:
    d = Params()
    return getattr(p, name) == getattr(d, name)


def autofill(p: Params, design, workers: int, log=print) -> Params:
    """Fill in whatever the user left at its default, for this sequence and budget."""
    n = design.n
    set_fields = []

    def put(name, value):
        if _at_default(p, name) and getattr(p, name) != value:
            setattr(p, name, value)
            set_fields.append(f"{name}={value}")

    # --- codon fitness: a floor that must hold, a target where it stops pushing ---
    put("cai_floor", 0.75)
    put("cai_target", 0.90)
    put("min_codon_freq", 0.08)       # drop codons too rare to translate well in either species

    # --- GC band: the lowest GC still reachable at the CAI this run aims for, plus headroom ---
    gcinfo = achievable_gc(design, p.cai_target)
    lo, win = gcinfo["gc_lo_feasible"], gcinfo["win_lo_feasible"]
    opt = gcinfo["gc_at_max_cai"]
    # `lo` is a single extreme point -- the GC-poorest sequence that still hits the target CAI.
    # A multi-objective search never sits exactly there, so a band drawn just above it is
    # unsatisfiable in practice. Put the ceiling halfway between that extreme and what pure
    # codon optimisation would give: a real compromise, and reachable.
    ceiling = lo + 0.5 * max(0.0, opt - lo)
    put("gc_min", round(max(0.25, lo - 0.12), 3))
    put("gc_max", round(min(0.65, max(ceiling, lo + 0.05)), 3))
    put("gc_window_max", round(min(0.85, win + 0.06), 3))
    put("gc_window_min", round(max(0.15, lo - 0.22), 3))

    # --- structure: thresholds with a mechanism behind them ---
    put("dsrna_min_bp", 30)          # PKR activation needs roughly this much contiguous duplex
    put("stem_target", 20)
    put("stem_cap", 30)
    put("min_paired_fraction", 0.10)  # a fully unstructured RNA hydrolyses faster
    put("start_unpaired_min_prob", 0.70)
    put("start_context_mean_min", 0.60)

    # --- folding cost settings ---
    put("inner_span", 0 if n <= 900 else 200)
    put("access_window", 200)
    put("access_span", 150)

    # --- size the search to the budget ---
    factor = calibrate(design.rna, p, design.cds)
    budget_s = max(30.0, p.time_budget_min * 60.0)
    cfg = choose_search(n, workers, budget_s, factor)
    rounds, gens, pop, elite, top_n, polish = cfg
    for name, value in (("rounds", rounds), ("generations_per_round", gens),
                        ("population", pop), ("elite_k", elite), ("top_n", top_n),
                        ("polish_codons", polish)):
        if _at_default(p, name):
            setattr(p, name, value)
    est = estimate_seconds(cfg, n, workers, factor)
    log(f"auto: machine is {1 / factor:.2f}x the reference speed; "
        f"{n} nt, {workers} worker(s), budget {p.time_budget_min:.0f} min")
    log(f"auto: search = {rounds} cycles x {gens} generations x {pop} candidates "
        f"(elites {elite}, report {top_n}, polish {polish}) -> estimated {est / 60:.1f} min")
    log(f"auto: at CAI {gcinfo['cai_at']:.2f} the lowest reachable GC is {lo:.3f} "
        f"(worst {p.gc_window} nt window {win:.3f}); codon-optimal would be "
        f"{gcinfo['gc_at_max_cai']:.3f} -> band {p.gc_min}-{p.gc_max}, "
        f"window max {p.gc_window_max}")
    if set_fields:
        log("auto: " + ", ".join(set_fields))
    p.validate()
    return p
