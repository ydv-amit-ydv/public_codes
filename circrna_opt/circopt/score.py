"""Composite objective: codon fitness + structure penalties + sequence-composition penalties.

Lower total == better. Every term is a non-negative penalty so weights read as "importance".
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .codon import NO_SCORE_AA, CODON_TO_AA
from .config import Params
from .constraints import Rules, gc_stats, inverted_repeat_count
from .design import Design
from .structure import Metrics


@dataclass
class Scored:
    genome: list = field(default_factory=list)
    seq: str = ""
    total: float = 0.0
    parts: dict = field(default_factory=dict)
    cai: dict = field(default_factory=dict)      # per species
    cai_gm: float = 0.0
    gc: dict = field(default_factory=dict)
    motifs: dict = field(default_factory=dict)
    metrics: Metrics | None = None
    feasible: bool = True
    violations: list = field(default_factory=list)

    def row(self, p: Params) -> dict:
        d = {"total_score": self.total, "feasible": self.feasible,
             "cai_gm": self.cai_gm, **{f"cai_{k}": v for k, v in self.cai.items()},
             **self.gc, "motif_violations": sum(self.motifs.values()),
             "motifs": ";".join(f"{k}x{v}" for k, v in sorted(self.motifs.items())),
             "violations": ";".join(self.violations)}
        if self.metrics:
            d.update(self.metrics.to_row())
        d.update({f"pen_{k}": v for k, v in sorted(self.parts.items())})
        return d


def codon_fitness(design: Design, seq: str) -> tuple[dict, float]:
    codons = design.cds_codons(seq)
    p = design.p
    cai = {sp: design.model.index(codons, sp) for sp in design.model.species}
    wp = p.primary_weight
    gm = math.exp(wp * math.log(max(cai[p.primary_species], 1e-9))
                  + (1 - wp) * math.log(max(cai[p.secondary_species], 1e-9)))
    return cai, gm


def score(design: Design, genome: list[int], metrics: Metrics, rules: Rules,
          weights: dict | None = None, seq: str | None = None) -> Scored:
    p = design.p
    w = weights or p.weights
    seq = seq if seq is not None else design.build(genome)
    n = len(seq)
    per100 = max(n / 100.0, 1.0)
    s = Scored(genome=list(genome), seq=seq, metrics=metrics)

    # --- codon fitness ---
    s.cai, s.cai_gm = codon_fitness(design, seq)
    weaker = min(s.cai.values())
    pen = {"codon": w["codon"] * max(0.0, p.cai_target - s.cai_gm),
           "cai_floor": w["cai_floor"] * max(0.0, p.cai_floor - weaker)}

    # --- structure (MFE level) ---
    # Long contiguous duplex is the immunogenicity/translation proxy that actually has a
    # mechanism behind it (PKR needs ~30+ bp, MDA5 senses long duplexes), so it carries the
    # weight; total pairing below is deliberately soft.
    pen["dsrna"] = w["dsrna"] * (metrics.dsrna_nt / per100 + metrics.dsrna_helices)
    goal = p.structure_goal
    if goal == "minimize_pairing":
        pen["mfe"] = w["mfe"] * max(0.0, -metrics.mfe_per_nt)
    elif goal == "minimize_mfe":
        pen["mfe"] = w["mfe"] * max(0.0, metrics.mfe_per_nt - p.mfe_per_nt_target)
    else:
        pen["mfe"] = w["mfe"] * abs(metrics.mfe_per_nt - p.mfe_per_nt_target)
    pen["paired"] = w["paired"] * metrics.paired_fraction
    span = max(p.stem_cap - p.stem_target, 1)
    pen["stem"] = w["stem"] * max(0, metrics.longest_stem - p.stem_target) / span
    pen["long_range"] = w["long_range"] * metrics.long_range_pairs / per100
    start_len = p.start_upstream + 3 + p.start_downstream
    pen["start"] = w["start"] * metrics.start_paired / max(start_len, 1)
    pen["early"] = w["early"] * metrics.early_paired_fraction
    pen["min_paired"] = w["min_paired"] * max(0.0, p.min_paired_fraction - metrics.paired_fraction)
    # initiation accessibility: use the exact PF values when we have them, else the fast ones
    codon_acc = (metrics.start_codon_unpaired if metrics.start_codon_unpaired is not None
                 else metrics.start_codon_unpaired_fast)
    ctx_acc = (metrics.start_context_unpaired if metrics.start_context_unpaired is not None
               else metrics.start_context_unpaired_fast)
    if codon_acc is not None:
        pen["start_access"] = w["start_access"] * max(0.0, p.start_unpaired_min_prob - codon_acc)
    if ctx_acc is not None:
        pen["start_context"] = w["start_context"] * max(0.0, p.start_context_mean_min - ctx_acc)
    if metrics.early_unpaired_fast is not None:
        pen["early_access"] = w["early_access"] * (1.0 - metrics.early_unpaired_fast)

    # --- ensemble (only when the partition function was computed) ---
    if metrics.full:
        pen["entropy"] = w["entropy"] * (metrics.mean_entropy or 0.0)
        pen["diversity"] = w["diversity"] * (metrics.ensemble_diversity or 0.0) / n
        pen["centroid"] = w["centroid"] * (metrics.centroid_dist_to_mfe or 0) / n
        gap = (metrics.mfe - (metrics.efe or metrics.mfe))   # >= 0; per nt
        pen["mfe_freq"] = w["mfe_freq"] * min(1.0, abs(gap) / max(0.02 * n, 1e-9))
        pen["red"] = w["red"] * (metrics.red_pairs / per100
                                 + metrics.red_longest_run / max(p.stem_target, 1)
                                 + metrics.red_cds_crossing)
        pen["start_prob"] = w["start_prob"] * max(
            0.0, p.start_unpaired_min_prob - (metrics.start_unpaired_prob or 1.0))

    # --- composition ---
    s.gc = gc_stats(seq, p.gc_window, design.free_mask)
    gc_dev = max(0.0, s.gc["gc"] - p.gc_max) + max(0.0, p.gc_min - s.gc["gc"])
    win_dev = (max(0.0, s.gc["gc_win_max_free"] - p.gc_window_max)
               + max(0.0, p.gc_window_min - s.gc["gc_win_min_free"]))
    pen["gc"] = w["gc"] * gc_dev + w["gc_low"] * win_dev
    s.motifs = rules.fixable(seq, design.free_mask)
    pen["motif"] = w["motif"] * sum(s.motifs.values())
    pen["inverted_repeat"] = w["inverted_repeat"] * inverted_repeat_count(
        seq, p.ir_kmer, design.free_mask) / per100
    pen["ggg"] = w["ggg"] * rules.ggg(seq) / per100

    # --- hard constraints ---
    v = []
    if metrics.longest_stem >= p.stem_cap:
        v.append(f"stem>={p.stem_cap}({metrics.longest_stem})")
    if s.motifs:
        v.append("forbidden_motif")
    if weaker < p.cai_floor:
        v.append(f"cai_{min(s.cai, key=s.cai.get)}<{p.cai_floor}({weaker:.3f})")
    if metrics.paired_fraction < p.min_paired_fraction:
        v.append("too_unstructured")
    if not (p.gc_min <= s.gc["gc"] <= p.gc_max):
        v.append(f"GC={s.gc['gc']:.3f} outside [{p.gc_min},{p.gc_max}]")
    # only windows the design can actually influence
    if s.gc["gc_win_max_free"] > p.gc_window_max:
        v.append(f"GC_window_max={s.gc['gc_win_max_free']:.3f}>{p.gc_window_max}")
    if s.gc["gc_win_min_free"] < p.gc_window_min:
        v.append(f"GC_window_min={s.gc['gc_win_min_free']:.3f}<{p.gc_window_min}")
    if codon_acc is not None and codon_acc < p.start_unpaired_min_prob:
        v.append(f"AUG_paired(p_unpaired={codon_acc:.2f}<{p.start_unpaired_min_prob})")
    if ctx_acc is not None and ctx_acc < p.start_context_mean_min:
        v.append(f"start_context_paired(mean={ctx_acc:.2f}<{p.start_context_mean_min})")
    s.violations = v
    s.feasible = not v
    s.parts = pen
    s.total = sum(pen.values())
    return s
