"""Ties the pieces together: DERNA seeds -> iterative optimiser -> reports."""
from __future__ import annotations

import os
import time

from . import derna as derna_mod
from . import seqio
from .auto import autofill
from .codon import CodonModel
from .config import Params
from .design import Design
from .optimize import Optimizer
from .report import write_outputs


def resolve_inputs(p: Params, log=lambda *a: None) -> Params:
    """Turn whatever was pasted into `rna` + CDS bounds. Safe to call on a throwaway copy
    to validate a form before a run is started."""
    if p.sequence and not p.rna:
        model = CodonModel(species=(p.primary_species, p.secondary_species),
                           metric=p.codon_metric, trna_file=p.trna_copy_file)
        resolved = seqio.resolve(p, model)
        for note in resolved.pop("notes"):
            log(f"input: {note}")
        for k, v in resolved.items():
            setattr(p, k, v)
        p.validate()
    return p


def run_design(p: Params, outdir: str, log=print) -> dict:
    p.validate()
    t0 = time.time()
    p.protein_was_given = bool(p.protein or p.rna)

    # 1. turn whatever was pasted into a molecule + CDS, saying what was assumed
    resolve_inputs(p, log)
    design = Design(p)

    # 2. fill in everything the user left alone, sized to the time budget
    if p.auto:
        workers = p.workers or (os.cpu_count() or 1)
        autofill(p, design, workers, log=log)
        design = Design(p)        # min_codon_freq can change the codon options
    log(f"molecule: {design.n} nt, CDS {design.cds}, "
        f"{len(design.slots)}/{(design.cds_end - design.cds_start) // 3} codons designable, "
        f"protein {len(design.protein.rstrip('*'))} aa")
    log(f"species: primary={p.primary_species} (w={p.primary_weight}) "
        f"secondary={p.secondary_species}; metric={p.codon_metric}; floor={p.cai_floor}")
    log(f"structure goal: {p.structure_goal} (dsRNA threshold {p.dsrna_min_bp} bp)")

    opt = Optimizer(p, design=design, log=log)
    extra_seeds, derna_info = [], []
    if p.use_derna:
        for d in derna_mod.run_derna(design, log=log):
            try:
                g = design.genome_of(d["full_seq"])
            except Exception as e:
                log(f"DERNA: could not graft a seed ({e})")
                continue
            extra_seeds.append(g)
            derna_info.append({k: v for k, v in d.items() if k != "full_seq"})

    result = opt.run(extra_seeds=extra_seeds)
    result["derna"] = derna_info
    if derna_info:
        # score DERNA's own designs on our objective so the comparison arm is explicit
        for d in derna_info:
            d["note"] = ("DERNA minimises MFE (more structure); scored here on this run's "
                         "objective for comparison only")
    paths = write_outputs(result, p, outdir, log=log)
    result["paths"] = paths
    log(f"done in {time.time() - t0:.1f}s -> {os.path.abspath(outdir)}")
    b = result["best"]
    o = result["original"]
    log("")
    # when the input was a protein there is no "original" construct -- only our seed
    label = "reverse-translated seed" if (p.sequence and not p.protein_was_given) else "original"
    log(f"== {label} -> best ==")
    log(f"  CAI {p.primary_species}:  {o.cai[p.primary_species]:.3f} -> {b.cai[p.primary_species]:.3f}")
    log(f"  CAI {p.secondary_species}: {o.cai[p.secondary_species]:.3f} -> {b.cai[p.secondary_species]:.3f}")
    log(f"  MFE/nt:        {o.metrics.mfe_per_nt:+.3f} -> {b.metrics.mfe_per_nt:+.3f}")
    log(f"  paired frac:   {o.metrics.paired_fraction:.3f} -> {b.metrics.paired_fraction:.3f}")
    log(f"  dsRNA nt:      {o.metrics.dsrna_nt} -> {b.metrics.dsrna_nt}")
    log(f"  longest stem:  {o.metrics.longest_stem} -> {b.metrics.longest_stem} bp")
    log(f"  GC:            {o.gc['gc']:.3f} -> {b.gc['gc']:.3f}")
    log(f"  AUG p(unpaired):     {o.metrics.start_codon_unpaired:.3f} -> {b.metrics.start_codon_unpaired:.3f}"
        f"  (need >= {p.start_unpaired_min_prob})")
    log(f"  start ctx mean:      {o.metrics.start_context_unpaired:.3f} -> {b.metrics.start_context_unpaired:.3f}"
        f"  (need >= {p.start_context_mean_min})")
    log(f"  longest high-prob stem: {o.metrics.red_longest_run} -> {b.metrics.red_longest_run} bp "
        f"(pairs p>={p.red_prob})")
    log(f"  feasible:      {o.feasible} -> {b.feasible} {b.violations}")
    return result
