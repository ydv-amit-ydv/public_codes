"""The iterative engine.

Structure of a run (this is the "several cycles" part):

  round 1..R:
      * GA over synonymous codon choices, scored with a FAST fold
        (circular, max_bp_span = inner_span, MFE only)
      * the round's elites are re-folded EXACTLY (circular, unlimited span,
        MFE + partition function) and re-scored on the full objective
      * the circular Boltzmann ensemble of the round's best is sampled
        -> that is the per-cycle "ensemble of structures scanned"
      * weights are re-balanced from what the elites actually failed on
        (weak-species CAI, dsRNA, GC, motifs), then the next round starts
        from those elites
  final: top_n candidates re-folded exactly, optional noLP/noGU stress test
         and a linear-mode comparison, ranked, best returned.

The fast/exact split is what makes a 3000 nt molecule tractable: the GA does
thousands of cheap folds, and only a handful of exact partition functions run.
"""
from __future__ import annotations

import math
import os
import random
import time
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field

from .config import Params
from .constraints import Rules, gc_stats
from .design import Design
from .score import Scored, score
from .structure import fold

_W = {}   # per-process worker state


def _init_worker(p: Params):
    _W["p"] = p
    _W["design"] = Design(p)
    _W["rules"] = Rules(p)


def _eval_fast(args):
    genome, weights = args
    p, d, r = _W["p"], _W["design"], _W["rules"]
    seq = d.build(genome)
    m = fold(seq, p, d.cds, circ=True, full=False, local_start=True)   # span = p.inner_span
    s = score(d, genome, m, r, weights, seq=seq)
    return s.genome, s.total, s.feasible


def _eval_exact(args):
    genome, weights, n_sample = args
    p, d, r = _W["p"], _W["design"], _W["rules"]
    seq = d.build(genome)
    out = fold(seq, p, d.cds, circ=True, full=True, span=0, sample=n_sample)
    m, samples = out if n_sample else (out, [])
    s = score(d, genome, m, r, weights, seq=seq)
    return s, samples


def _gc_count(codon: str) -> int:
    return sum(1 for b in codon if b in "GC")


def gc_excess(seq: str, p: Params, free_mask) -> float:
    """Total GC-window violation mass over windows containing a designable position.

    The *max* window is a plateau -- several windows usually tie at the maximum, so no single
    codon swap lowers it and a greedy repair stalls immediately. Summing every window's excess
    over the configured band gives a strictly decreasing objective instead.
    """
    n = len(seq)
    w = min(p.gc_window, n)
    gc = [1 if b in "GC" else 0 for b in seq]
    fm = [1 if x else 0 for x in free_mask]
    ext_gc = gc + gc[: w - 1]
    ext_fm = fm + fm[: w - 1]
    cur = sum(ext_gc[:w])
    cfm = sum(ext_fm[:w])
    total = 0.0
    for i in range(n):
        if i:
            cur += ext_gc[i + w - 1] - ext_gc[i - 1]
            cfm += ext_fm[i + w - 1] - ext_fm[i - 1]
        if cfm:
            f = cur / w
            total += max(0.0, f - p.gc_window_max) + max(0.0, p.gc_window_min - f)
    return total


def repair(design, genome: list[int], rules: Rules, max_iter: int = 600) -> list[int]:
    """Greedily fix the LOCAL, deterministic constraints: GC windows and forbidden motifs.

    These are exactly the constraints a random GA is bad at -- they depend on a handful of
    specific positions, so blind mutation almost never lands the fix, while a directed swap
    always can. Folding is never called here, so this is cheap enough to run on every elite.
    The weaker species' CAI is held at or above `cai_floor` throughout, so repair cannot trade
    away the thing the floor exists to protect.
    """
    p = design.p
    g = list(genome)
    n = design.n
    model = design.model
    species = (p.primary_species, p.secondary_species)

    def weak_cai(seq):
        return min(model.index(design.cds_codons(seq), sp) for sp in species)

    def cost(seq):
        return (sum(rules.fixable(seq, design.free_mask).values()),
                gc_excess(seq, p, design.free_mask))

    slot_at = {}
    for k, pos in enumerate(design.slots):
        for off in range(3):
            slot_at[pos + off] = k

    seq = design.build(g)
    cur_cost = cost(seq)
    floor = min(weak_cai(seq), p.cai_floor)

    for _ in range(max_iter):
        if cur_cost == (0, 0.0):
            break
        # choose the positions to work on: a motif occurrence, else the worst GC window
        targets: list[int] = []
        prefer_low_gc = None
        motifs = rules.fixable(seq, design.free_mask)
        if motifs:
            for name, start, ln in rules.spans(seq):
                if name in motifs:
                    span = [(start + i) % n for i in range(ln)]
                    if any(design.free_mask[i] for i in span):
                        targets = [i for i in span if i in slot_at]
                        break
        if not targets:
            w = min(p.gc_window, n)
            gc = [1 if b in "GC" else 0 for b in seq]
            best_i, best_v, high = None, None, True
            for i in range(n):
                win = [(i + k) % n for k in range(w)]
                if not any(design.free_mask[x] for x in win):
                    continue
                f = sum(gc[x] for x in win) / w
                over = f - p.gc_window_max
                under = p.gc_window_min - f
                v = max(over, under)
                if v > 0 and (best_v is None or v > best_v):
                    best_i, best_v, high = i, v, over >= under
            if best_i is None:
                break
            prefer_low_gc = high
            targets = [x for x in ((best_i + k) % n for k in range(w)) if x in slot_at]

        improved = False
        for k in sorted({slot_at[i] for i in targets}):
            opts = design.options[k]
            cur = g[k]
            order = list(range(len(opts)))
            if prefer_low_gc is True:
                order.sort(key=lambda i: _gc_count(opts[i]))
            elif prefer_low_gc is False:
                order.sort(key=lambda i: -_gc_count(opts[i]))
            for alt in order:
                if alt == cur:
                    continue
                g[k] = alt
                seq2 = design.build(g)
                c2 = cost(seq2)
                if c2 < cur_cost and weak_cai(seq2) >= floor:
                    cur_cost, seq, improved = c2, seq2, True
                    break
                g[k] = cur
            if improved:
                break
        if not improved:
            break
    return g


@dataclass
class RoundLog:
    round: int = 0
    generations: int = 0
    best_total: float = 0.0
    best_row: dict = field(default_factory=dict)
    weights: dict = field(default_factory=dict)
    ensemble: list = field(default_factory=list)      # sampled circular structures of the best
    ensemble_stats: dict = field(default_factory=dict)
    elapsed_s: float = 0.0
    notes: list = field(default_factory=list)


class Optimizer:
    def __init__(self, p: Params, design: Design | None = None, log=print):
        self.p = p
        self.design = design or Design(p)
        self.rules = Rules(p)
        self.rng = random.Random(p.seed)
        self.weights = dict(p.weights)
        self.log = log
        self.rounds: list[RoundLog] = []
        self.pareto: list[Scored] = []
        self._pool: ProcessPoolExecutor | None = None

    # ---- population helpers -------------------------------------------------
    def _random_genome(self) -> list[int]:
        return [self.rng.randrange(len(o)) for o in self.design.options]

    def _best_codon_genome(self, species: str | None) -> list[int]:
        """Greedy per-codon: pick the highest-weight codon (one species, or the dual geometric mean)."""
        d, p = self.design, self.p
        g = []
        for k, opts in enumerate(d.options):
            if species:
                vals = [d.model.logw(species, c) for c in opts]
            else:
                vals = [p.primary_weight * d.model.logw(p.primary_species, c)
                        + (1 - p.primary_weight) * d.model.logw(p.secondary_species, c) for c in opts]
            g.append(max(range(len(opts)), key=lambda i: vals[i]))
        return g

    def _seed_population(self, extra: list[list[int]]) -> list[list[int]]:
        d = self.design
        pop = [d.genome_of()]                                  # the user's own sequence
        pop.append(self._best_codon_genome(None))              # dual-species greedy
        pop.append(self._best_codon_genome(self.p.primary_species))
        pop.append(self._best_codon_genome(self.p.secondary_species))
        for g in extra:
            pop.append(g)
        for s in self.p.seed_sequences:
            try:
                pop.append(d.genome_of(s))
            except Exception:
                self.log("  (skipped an unusable seed_sequence)")
        # mutated copies of the good seeds + random fill
        base = list(pop)
        while len(pop) < self.p.population:
            g = list(self.rng.choice(base))
            self._mutate(g, rate=max(self.p.mutation_rate, 0.05))
            pop.append(g)
        return pop[: self.p.population]

    def _mutate(self, g: list[int], rate: float | None = None) -> None:
        rate = self.p.mutation_rate if rate is None else rate
        n = len(g)
        if not n:
            return
        k = max(1, int(round(rate * n)))
        for _ in range(k):
            i = self.rng.randrange(n)
            opts = self.design.options[i]
            if len(opts) > 1:
                g[i] = self.rng.randrange(len(opts))

    def _crossover(self, a: list[int], b: list[int]) -> list[int]:
        n = len(a)
        if n < 2:
            return list(a)
        i, j = sorted(self.rng.sample(range(n), 2))
        return a[:i] + b[i:j] + a[j:]

    # ---- the GA ------------------------------------------------------------
    def _ga(self, pop: list[list[int]], generations: int, rnd: int) -> list[tuple[list[int], float, bool]]:
        p = self.p
        scored = self._map_fast(pop)
        for gen in range(generations):
            scored.sort(key=lambda t: (not t[2], t[1]))
            elites = [t[0] for t in scored[: p.elite_k]]
            kids: list[list[int]] = []
            while len(kids) < p.population - len(elites):
                a = self._tournament(scored)
                b = self._tournament(scored)
                c = self._crossover(a, b)
                self._mutate(c)
                kids.append(c)
            scored = self._map_fast(elites + kids)
            best = min(scored, key=lambda t: (not t[2], t[1]))
            if gen == generations - 1 or gen % 4 == 0:
                nfeas = sum(1 for t in scored if t[2])
                self.log(f"  round {rnd} gen {gen + 1}/{generations}: best={best[1]:.4f} feasible={nfeas}/{len(scored)}")
        scored.sort(key=lambda t: (not t[2], t[1]))
        return scored

    def _tournament(self, scored, k=3) -> list[int]:
        cand = [self.rng.choice(scored) for _ in range(k)]
        return list(min(cand, key=lambda t: (not t[2], t[1]))[0])

    def _map_fast(self, pop):
        args = [(g, self.weights) for g in pop]
        if self._pool:
            return list(self._pool.map(_eval_fast, args, chunksize=1))
        return [_eval_fast(a) for a in args]

    def _map_exact(self, genomes, n_sample=0):
        args = [(g, self.weights, n_sample) for g in genomes]
        if self._pool:
            return list(self._pool.map(_eval_exact, args, chunksize=1))
        return [_eval_exact(a) for a in args]

    # ---- weight re-balancing between rounds --------------------------------
    def _rebalance(self, elites: list[Scored], rl: RoundLog) -> None:
        """Push weight toward whatever the elites are still failing, so the rounds converge
        on a *balanced* solution instead of over-serving one objective."""
        p = self.p
        best = elites[0]
        notes = []
        weaker_sp = min(best.cai, key=best.cai.get)
        weaker = best.cai[weaker_sp]
        if weaker < p.cai_floor:
            self.weights["cai_floor"] *= 1.6
            notes.append(f"{weaker_sp} CAI {weaker:.3f} < floor {p.cai_floor} -> cai_floor w={self.weights['cai_floor']:.2f}")
        elif weaker > p.cai_floor + 0.08 and self.weights["cai_floor"] > p.weights["cai_floor"]:
            self.weights["cai_floor"] = max(p.weights["cai_floor"], self.weights["cai_floor"] / 1.3)
            notes.append("weak-species CAI comfortable -> relaxing cai_floor")
        m = best.metrics
        if m and m.dsrna_nt:
            self.weights["dsrna"] *= 1.5
            notes.append(f"{m.dsrna_nt} nt in >={p.dsrna_min_bp} bp duplex -> dsrna w={self.weights['dsrna']:.2f}")
        if m and m.longest_stem > p.stem_target:
            self.weights["stem"] *= 1.3
            notes.append(f"longest stem {m.longest_stem} > {p.stem_target} -> stem w={self.weights['stem']:.2f}")
        if best.motifs:
            self.weights["motif"] *= 2.0
            notes.append(f"motifs {best.motifs} -> motif w={self.weights['motif']:.2f}")
        if not (p.gc_min <= best.gc["gc"] <= p.gc_max):
            self.weights["gc"] *= 1.4
            notes.append(f"GC {best.gc['gc']:.3f} outside [{p.gc_min},{p.gc_max}] -> gc w={self.weights['gc']:.2f}")
        # Use the same measures the hard constraints use. (An earlier version read the strict
        # minimum over the whole start context, which is ~0 even in a good design, so these
        # weights inflated every round for no reason.)
        if m and m.start_codon_unpaired is not None and m.start_codon_unpaired < p.start_unpaired_min_prob:
            self.weights["start_access"] *= 1.5
            self.weights["start_prob"] *= 1.3
            notes.append(f"AUG p_unpaired {m.start_codon_unpaired:.2f} "
                         f"< {p.start_unpaired_min_prob} -> start weights up")
        if m and m.start_context_unpaired is not None and m.start_context_unpaired < p.start_context_mean_min:
            self.weights["start_context"] *= 1.4
            notes.append(f"start-context mean p_unpaired {m.start_context_unpaired:.2f} "
                         f"< {p.start_context_mean_min} -> start_context weight up")
        cap = 200.0
        for k in self.weights:
            self.weights[k] = min(self.weights[k], cap)
        rl.notes = notes
        for nt in notes:
            self.log(f"  rebalance: {nt}")

    def _update_pareto(self, cands: list[Scored]) -> None:
        """Keep the non-dominated set on (weak-species CAI up, dsRNA nt down, MFE/nt up)."""
        def key(s):
            return (min(s.cai.values()), -(s.metrics.dsrna_nt if s.metrics else 0),
                    s.metrics.mfe_per_nt if s.metrics else 0.0)
        pool = self.pareto + [c for c in cands if c.metrics]
        keep = []
        for a in pool:
            ka = key(a)
            if any(all(kb >= kaa for kb, kaa in zip(key(b), ka)) and key(b) != ka for b in pool):
                continue
            if any(c.seq == a.seq for c in keep):
                continue
            keep.append(a)
        keep.sort(key=lambda s: s.total)
        self.pareto = keep[:40]

    # ---- driver ------------------------------------------------------------
    def run(self, extra_seeds: list[list[int]] | None = None) -> dict:
        p, d = self.p, self.design
        self.log(f"design: {d.summary()}")
        if not d.slots:
            raise ValueError("no designable codons: the CDS is fully protected or has no synonyms")
        workers = p.workers or (os.cpu_count() or 1)
        pop = self._seed_population(extra_seeds or [])
        t_all = time.time()
        try:
            if workers > 1:
                self._pool = ProcessPoolExecutor(max_workers=workers,
                                                 initializer=_init_worker, initargs=(p,))
            else:
                _init_worker(p)
            all_elites: list[Scored] = []
            for rnd in range(1, p.rounds + 1):
                t0 = time.time()
                rl = RoundLog(round=rnd, generations=p.generations_per_round)
                scored = self._ga(pop, p.generations_per_round, rnd)
                elite_genomes = [t[0] for t in scored[: p.elite_k]]
                if p.repair:
                    fixed = [repair(self.design, g, self.rules) for g in elite_genomes]
                    n_changed = sum(1 for a, b in zip(elite_genomes, fixed) if a != b)
                    if n_changed:
                        self.log(f"  round {rnd}: repaired GC/motif violations in "
                                 f"{n_changed}/{len(fixed)} elites")
                    elite_genomes = fixed
                self.log(f"  round {rnd}: exact circular fold + partition function on {len(elite_genomes)} elites")
                res = self._map_exact(elite_genomes, n_sample=0)
                elites = [s for s, _ in res]
                elites.sort(key=lambda s: (not s.feasible, s.total))
                all_elites.extend(elites)
                self._update_pareto(elites)

                best = elites[0]
                rl.best_total = best.total
                rl.best_row = best.row(p)
                # per-cycle ensemble scan of the round's best
                _, samples = self._map_exact([best.genome], n_sample=10)[0]
                rl.ensemble = samples
                rl.ensemble_stats = _ensemble_stats(samples, best.metrics)
                rl.weights = dict(self.weights)
                rl.elapsed_s = time.time() - t0
                self.log(f"  round {rnd} best: total={best.total:.4f} "
                         f"CAI={ {k: round(v, 3) for k, v in best.cai.items()} } "
                         f"MFE/nt={best.metrics.mfe_per_nt:.3f} dsRNA={best.metrics.dsrna_nt}nt "
                         f"stem={best.metrics.longest_stem} feasible={best.feasible} "
                         f"({rl.elapsed_s:.1f}s)")
                self._rebalance(elites, rl)
                self.rounds.append(rl)
                # next round starts from the elites, diversified
                pop = [list(g) for g in elite_genomes]
                while len(pop) < p.population:
                    g = list(self.rng.choice(elite_genomes))
                    self._mutate(g, rate=max(p.mutation_rate, 0.03))
                    pop.append(g)

            # ---- final: rank, stress test, linear comparison ----
            # The adaptive weights differ from round to round, so totals carried over from
            # different rounds are not comparable. Everything that gets ranked or reported is
            # re-scored here on the USER's weights: the adaptation only steers the search.
            uniq: dict[str, Scored] = {}
            for s in all_elites:
                uniq.setdefault(s.seq, s)
            self.weights = dict(p.weights)
            rescored = [s for s, _ in self._map_exact([s.genome for s in uniq.values()])]
            rescored.sort(key=lambda s: (not s.feasible, s.total))
            self._update_pareto(rescored)
            finalists = rescored[: p.top_n]
            self.log(f"final: {len(finalists)} unique candidates; polishing + stress tests")
            finalists = [self._polish(s) for s in finalists]
            finalists.sort(key=lambda s: (not s.feasible, s.total))
            extras = [self._extra_checks(s) for s in finalists]
            original = self._evaluate_sequence(d.rna)
            best = finalists[0]
            best_ens = self._map_exact([best.genome], n_sample=20)[0][1]
            return {"design": d.summary(), "rounds": self.rounds, "finalists": finalists,
                    "extras": extras, "original": original, "best": best,
                    "best_ensemble": best_ens,
                    "best_ensemble_stats": _ensemble_stats(best_ens, best.metrics),
                    "pareto": self.pareto, "weights_final": dict(self.weights),
                    "elapsed_s": time.time() - t_all}
        finally:
            if self._pool:
                self._pool.shutdown()
                self._pool = None

    def _polish(self, s: Scored) -> Scored:
        """Greedy single-codon hill-climb on the exact objective (few, expensive evaluations)."""
        p = self.p
        if p.polish_codons <= 0:
            return s
        cur = s
        for _ in range(p.polish_codons):
            cands = []
            for k in self.rng.sample(range(len(cur.genome)), min(12, len(cur.genome))):
                opts = self.design.options[k]
                for alt in range(len(opts)):
                    if alt != cur.genome[k]:
                        g = list(cur.genome)
                        g[k] = alt
                        cands.append(g)
            if not cands:
                break
            self.rng.shuffle(cands)
            cands = cands[: max(4, (self.p.workers or os.cpu_count() or 4))]
            res = self._map_exact(cands)
            best = min((s2 for s2, _ in res), key=lambda s2: (not s2.feasible, s2.total))
            if (not best.feasible, best.total) < (not cur.feasible, cur.total):
                cur = best
            else:
                break
        return cur

    def _evaluate_sequence(self, seq: str) -> Scored:
        """Score an arbitrary full-length molecule (used for the original and for DERNA seeds)."""
        g = self.design.genome_of(seq)
        return self._map_exact([g])[0][0]

    def _extra_checks(self, s: Scored) -> dict:
        p, d = self.p, self.design
        out: dict = {}
        if p.linear_compare:
            m = fold(s.seq, p, d.cds, circ=False, full=True, span=0)
            out["linear"] = {"mfe": m.mfe, "mfe_per_nt": m.mfe_per_nt,
                             "paired_fraction": m.paired_fraction, "dsrna_nt": m.dsrna_nt,
                             "longest_stem": m.longest_stem, "efe": m.efe,
                             "ensemble_diversity": m.ensemble_diversity}
        if p.stress_test:
            m = fold(s.seq, p, d.cds, circ=True, full=False, span=0, noLP=True, noGUclosure=True)
            out["stress_noLP_noGUclosure"] = {
                "mfe": m.mfe, "mfe_per_nt": m.mfe_per_nt, "paired_fraction": m.paired_fraction,
                "dsrna_nt": m.dsrna_nt, "longest_stem": m.longest_stem}
        if not p.gquad:
            m = fold(s.seq, p, d.cds, circ=True, full=False, span=0, gquad=True)
            out["gquad_on"] = {"mfe": m.mfe, "mfe_per_nt": m.mfe_per_nt,
                               "paired_fraction": m.paired_fraction, "longest_stem": m.longest_stem}
        return out


def _ensemble_stats(samples: list[str], mfe_metrics) -> dict:
    """Spread of the sampled circular ensemble: how consistent the fold is."""
    if not samples:
        return {}
    import RNA
    n = len(samples[0])
    pf = [sum(1 for c in s if c != ".") / n for s in samples]
    d_mfe = [RNA.bp_distance(s, mfe_metrics.structure) for s in samples] if mfe_metrics else []
    pair = [RNA.bp_distance(samples[i], samples[j])
            for i in range(len(samples)) for j in range(i + 1, len(samples))]
    return {"n": len(samples),
            "paired_fraction_mean": sum(pf) / len(pf),
            "paired_fraction_min": min(pf), "paired_fraction_max": max(pf),
            "bp_distance_to_mfe_mean": (sum(d_mfe) / len(d_mfe)) if d_mfe else None,
            "pairwise_bp_distance_mean": (sum(pair) / len(pair)) if pair else 0.0}
