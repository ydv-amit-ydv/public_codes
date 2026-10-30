"""ViennaRNA folding and structure metrics (circular by default, matching the RNAfold web server).

Two ViennaRNA 2.7 quirks this module works around:
  * `subopt()` SEGFAULTS in circular mode -- never call it with md.circ=1.
  * `pbacktrack()` silently returns nothing unless md.uniq_ML=1 (set via uniq_ML=True).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import RNA

from .config import Params

KT37 = 0.0019872041 * 310.15  # kcal/mol


def model_details(p: Params, circ=True, span=0, noLP=None, noGUclosure=None, gquad=None,
                  uniq_ML=False) -> "RNA.md":
    md = RNA.md()
    if uniq_ML:
        # Required for stochastic backtracking (pbacktrack); silently returns no
        # structures without it. Does not change MFE or ensemble energies.
        md.uniq_ML = 1
    md.temperature = p.temperature
    md.dangles = p.dangles
    md.circ = 1 if circ else 0
    md.noLP = 1 if (p.no_lonely_pairs if noLP is None else noLP) else 0
    md.noGU = 1 if p.no_gu else 0
    md.noGUclosure = 1 if (p.no_gu_closure if noGUclosure is None else noGUclosure) else 0
    md.gquad = 1 if (p.gquad if gquad is None else gquad) else 0
    if span:
        md.max_bp_span = span
    return md


def pair_table(db: str) -> list[int]:
    """0-based partner index per position, -1 if unpaired. Handles circular dot-bracket."""
    pt = [-1] * len(db)
    stack = []
    for i, ch in enumerate(db):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise ValueError("unbalanced dot-bracket string")
            j = stack.pop()
            pt[i], pt[j] = j, i
    if stack:
        raise ValueError("unbalanced dot-bracket string")
    return pt


def helices(pt: list[int]) -> list[tuple[int, int, int]]:
    """Maximal runs of stacked pairs as (start_i, start_j, length), circular-aware.

    A helix is (i,j), (i+1,j-1), ... -- i.e. a contiguous duplex. Length is in base pairs,
    so a helix of length L is L bp of double-stranded RNA.
    """
    n = len(pt)
    out, seen = [], set()
    for i in range(n):
        j = pt[i]
        if j < 0 or i in seen:
            continue
        prev = (i - 1) % n
        if pt[prev] >= 0 and pt[prev] == (j + 1) % n:
            continue                      # not the start of the helix
        length, a, b = 0, i, j
        while True:
            am, bm = a % n, b % n
            if pt[am] != bm or am in seen:
                break
            seen.add(am)
            seen.add(bm)
            length += 1
            a, b = a + 1, b - 1
            if length > n:
                break
        if length:
            out.append((i, j, length))
    return out


def longest_stem(pt: list[int]) -> int:
    """Longest contiguous duplex, in base pairs."""
    hs = helices(pt)
    return max((h[2] for h in hs), default=0)


def circ_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


@dataclass
class Metrics:
    """All structure numbers for one sequence. `full=False` means MFE-only (fast screening)."""
    mfe: float = 0.0
    mfe_per_nt: float = 0.0
    structure: str = ""
    paired_fraction: float = 0.0
    longest_stem: int = 0
    dsrna_nt: int = 0            # nt inside helices >= dsrna_min_bp (2 * bp)
    dsrna_helices: int = 0
    long_range_pairs: int = 0
    start_paired: int = 0
    early_paired_fraction: float = 0.0
    full: bool = False
    # partition-function metrics
    efe: float | None = None
    mfe_frequency: float | None = None
    ensemble_diversity: float | None = None
    centroid: str | None = None
    centroid_distance: float | None = None
    centroid_dist_to_mfe: int | None = None
    mean_entropy: float | None = None
    max_entropy: float | None = None
    red_pairs: int = 0
    red_longest_run: int = 0
    red_cds_crossing: int = 0
    # accessibility of the initiation region. *_fast comes from RNAplfold (every candidate),
    # the unsuffixed fields from the exact circular PF (elites only). Same definitions.
    start_codon_unpaired_fast: float | None = None
    start_context_unpaired_fast: float | None = None
    early_unpaired_fast: float | None = None
    start_codon_unpaired: float | None = None
    start_context_unpaired: float | None = None
    start_unpaired_prob: float | None = None    # min over the whole start context (strictest)
    early_unpaired_prob: float | None = None
    bp_prob_pairs: list = field(default_factory=list)   # (i, j, p) 0-based, p >= 1e-3
    unpaired_prob: list = field(default_factory=list)
    entropy: list = field(default_factory=list)
    mountain: list = field(default_factory=list)

    def to_row(self) -> dict:
        d = {k: v for k, v in self.__dict__.items()
             if k not in ("bp_prob_pairs", "unpaired_prob", "entropy", "mountain", "structure", "centroid")}
        return d


def circular_accessibility(seq: str, p: Params) -> list[float]:
    """Per-position p(unpaired) for the whole circle, RNAplfold-style (fast).

    A full circular partition function is far too slow to run on every GA candidate, so this
    uses `RNA.pfl_fold_up` (the RNAplfold algorithm: sliding window `access_window`, pair span
    `access_span`). The sequence is DOUBLED and the middle copy is read back, which gives every
    position -- including the ones at the junction -- a complete circular context instead of
    plfold's linear edge effects.

    Validated against the exact circular PF on the example construct: with window 600 / span 300
    the start-context profile matches to ~0.1 mean absolute error and reproduces near-zero
    accessibility correctly, in about a second for a 618 nt molecule.
    """
    n = len(seq)
    W = min(p.access_window, 2 * n)
    L = min(p.access_span, W)
    up = RNA.pfl_fold_up(seq + seq, 1, W, L)
    # pfl_fold_up is 1-based; take the middle copy so each position sees n/2 nt of context
    off = n // 2
    acc = [0.0] * n
    for k in range(n):
        v = up[off + k + 1][1]            # 1-based; position off+k of the doubled sequence
        acc[(off + k) % n] = min(max(v, 0.0), 1.0)
    return acc


def start_region_accessibility(seq: str, p: Params, cds: tuple[int, int]) -> tuple[float, float, float]:
    """(min over AUG, mean over start context, mean over early ORF) of p(unpaired).

    Only a slice around the initiation region is folded. RNAplfold never pairs further than
    `access_span`, so a slice reaching `access_window` nt beyond the context on both sides gives
    the same answer as folding the whole molecule -- but at a cost independent of molecule
    length, which is what makes a 3000 nt design affordable inside the GA.
    """
    n = len(seq)
    s0, e0 = cds
    pad = p.access_window
    lo = s0 - p.start_upstream - pad
    hi = s0 + 3 + p.start_downstream + max(pad, p.early_orf_nt)
    if hi - lo >= n:                       # short molecule: use the whole circle
        acc = circular_accessibility(seq, p)
        codon_r, ctx_r, early_r = start_regions(p, n, cds)
        return (min(acc[i] for i in codon_r),
                sum(acc[i] for i in ctx_r) / len(ctx_r),
                (sum(acc[i] for i in early_r) / len(early_r)) if early_r else 1.0)
    idx = [(i % n) for i in range(lo, hi)]
    win = "".join(seq[i] for i in idx)
    W = min(p.access_window, len(win))
    L = min(p.access_span, W)
    up = RNA.pfl_fold_up(win, 1, W, L)
    def at(k):                              # k = offset from `lo`
        return min(max(up[k + 1][1], 0.0), 1.0)
    base = p.start_upstream + pad           # offset of the start codon inside the slice
    codon = [at(base + k) for k in range(3)]
    ctx = [at(base - p.start_upstream + k) for k in range(p.start_upstream + 3 + p.start_downstream)]
    n_early = min(p.early_orf_nt, e0 - s0)
    early = [at(base + k) for k in range(n_early)]
    return (min(codon), sum(ctx) / len(ctx), (sum(early) / len(early)) if early else 1.0)


def start_regions(p: Params, n: int, cds: tuple[int, int]) -> tuple[list[int], list[int], list[int]]:
    """(start codon positions, start context positions, early ORF positions), all mod n."""
    s, e = cds
    codon = [(s + k) % n for k in range(3)]
    ctx = [(i % n) for i in range(s - p.start_upstream, s + 3 + p.start_downstream)]
    early_hi = min(e, s + p.early_orf_nt)
    early = [(i % n) for i in range(s, early_hi)]
    return codon, ctx, early


def _regions(p: Params, n: int, cds: tuple[int, int]) -> tuple[range, range]:
    s, e = cds
    start_lo = max(0, s - p.start_upstream)
    start_hi = min(n, s + 3 + p.start_downstream)
    early_hi = min(e, s + p.early_orf_nt)
    return range(start_lo, start_hi), range(s, early_hi)


def fold(seq: str, p: Params, cds: tuple[int, int], *, circ=True, full=False,
         span: int | None = None, noLP=None, noGUclosure=None, gquad=None,
         sample: int = 0, local_start=False) -> Metrics | tuple[Metrics, list[str]]:
    """Fold `seq` and collect metrics. full=True adds the partition function.

    span=None uses p.inner_span for screening folds; pass 0 for an unlimited (exact) fold.
    With sample>0 also returns `sample` Boltzmann-sampled structures (needs full=True).
    """
    n = len(seq)
    md = model_details(p, circ=circ, span=(p.inner_span if span is None else span),
                       noLP=noLP, noGUclosure=noGUclosure, gquad=gquad, uniq_ML=bool(sample))
    fc = RNA.fold_compound(seq, md)
    db, mfe = fc.mfe()
    pt = pair_table(db)
    hs = helices(pt)
    long_h = [h for h in hs if h[2] >= p.dsrna_min_bp]
    m = Metrics(mfe=mfe, mfe_per_nt=mfe / n, structure=db, full=full,
                paired_fraction=sum(1 for x in pt if x >= 0) / n,
                longest_stem=max((h[2] for h in hs), default=0),
                dsrna_nt=2 * sum(h[2] for h in long_h), dsrna_helices=len(long_h))
    m.long_range_pairs = sum(1 for i, j in enumerate(pt)
                             if j > i and circ_dist(i, j, n) > p.long_range_nt)
    start_r, early_r = _regions(p, n, cds)
    m.start_paired = sum(1 for i in start_r if pt[i] >= 0)
    m.early_paired_fraction = (sum(1 for i in early_r if pt[i] >= 0) / len(early_r)) if len(early_r) else 0.0

    if local_start:
        (m.start_codon_unpaired_fast, m.start_context_unpaired_fast,
         m.early_unpaired_fast) = start_region_accessibility(seq, p, cds)

    if not full:
        return (m, []) if sample else m

    fc.exp_params_rescale(mfe)
    _, efe = fc.pf()
    m.efe = efe
    m.mfe_frequency = math.exp((efe - mfe) / KT37)
    m.ensemble_diversity = fc.mean_bp_distance()
    cen, cdist = fc.centroid()
    m.centroid, m.centroid_distance = cen, cdist
    m.centroid_dist_to_mfe = RNA.bp_distance(db, cen)

    pl = [(x.i - 1, x.j - 1, x.p) for x in fc.plist_from_probs(1e-5) if x.i > 0 and x.j > 0]
    unp = [1.0] * n
    ent = [0.0] * n
    for i, j, pr in pl:
        if not (0 <= i < n and 0 <= j < n):
            continue
        unp[i] -= pr
        unp[j] -= pr
        c = -pr * math.log(pr) if pr > 0 else 0.0
        ent[i] += c
        ent[j] += c
    for i in range(n):
        u = min(max(unp[i], 0.0), 1.0)
        unp[i] = u
        if u > 0:
            ent[i] -= u * math.log(u)
    m.unpaired_prob, m.entropy = unp, ent
    m.mean_entropy = sum(ent) / n
    m.max_entropy = max(ent)
    m.bp_prob_pairs = [(i, j, pr) for i, j, pr in pl if pr >= 1e-3]

    red = [(i, j) for i, j, pr in pl if pr >= p.red_prob]
    m.red_pairs = len(red)
    redset = set(red)
    runs, seen = 0, set()
    for i, j in red:
        if (i - 1, j + 1) in redset:
            continue
        r, a, b = 0, i, j
        while (a, b) in redset and (a, b) not in seen:
            seen.add((a, b)); r += 1; a, b = a + 1, b - 1
        runs = max(runs, r)
    m.red_longest_run = runs
    s, e = cds
    m.red_cds_crossing = sum(1 for i, j in red if (s <= i < e) != (s <= j < e))
    codon_r, ctx_r, _ = start_regions(p, n, cds)
    m.start_codon_unpaired = min(unp[i] for i in codon_r)
    m.start_context_unpaired = sum(unp[i] for i in ctx_r) / len(ctx_r)
    m.start_unpaired_prob = min((unp[i] for i in start_r), default=1.0)
    m.early_unpaired_prob = (sum(unp[i] for i in early_r) / len(early_r)) if len(early_r) else 1.0

    # mountain plot: #pairs enclosing each position (from the MFE structure)
    h, mt = 0, []
    for ch in db:
        if ch == "(":
            h += 1
        mt.append(h)
        if ch == ")":
            h -= 1
    m.mountain = mt

    if sample:
        try:
            structs = list(fc.pbacktrack(sample))
        except Exception:
            structs = []
        return m, structs
    return m
