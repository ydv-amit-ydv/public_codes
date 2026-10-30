"""Forbidden-motif, composition and complementary-repeat checks (all circular-aware)."""
from __future__ import annotations

import re

from .config import Params

IUPAC = {"A": "A", "C": "C", "G": "G", "U": "U", "T": "U", "R": "[AG]", "Y": "[CU]", "S": "[GC]",
         "W": "[AU]", "K": "[GU]", "M": "[AC]", "B": "[CGU]", "D": "[AGU]", "H": "[ACU]",
         "V": "[ACG]", "N": "[ACGU]"}
_RC = {"A": "U", "U": "A", "G": "C", "C": "G"}


def revcomp(s: str) -> str:
    return "".join(_RC[b] for b in reversed(s))


def _rx(motif: str) -> str:
    return "".join(IUPAC[b] for b in motif.upper())


class Rules:
    """Compiled rule set. `counts(seq)` returns {rule: n_violations}, wrapping across the circle."""

    def __init__(self, p: Params):
        self.p = p
        self.rules: list[tuple[str, re.Pattern]] = []
        self.pad = 12
        if p.avoid_re_sites:
            for name, site in p.re_sites.items():
                s = site.upper().replace("T", "U")
                pats = {s, revcomp(s)} if set(s) <= set("ACGU") else {s}
                self.rules.append((f"RE:{name}", re.compile("(?=(" + "|".join(_rx(x) for x in pats) + "))")))
                self.pad = max(self.pad, len(s))
        for m in p.extra_forbidden:
            self.rules.append((f"motif:{m}", re.compile(f"(?=({_rx(m)}))")))
            self.pad = max(self.pad, len(m))
        self.rules.append(("polyU", re.compile(f"(?=(U{{{p.max_poly_u + 1}}}))")))
        self.rules.append(("homopolymer", re.compile(f"(?=(A{{{p.max_homopolymer + 1}}}|C{{{p.max_homopolymer + 1}}}|G{{{p.max_homopolymer + 1}}}))")))
        if p.avoid_slippery:
            self.rules.append(("slippery", re.compile(r"(?=(([ACGU])\2\2([AU])\3\3[ACU]))")))
        if p.avoid_g4:
            self.rules.append(("G4", re.compile(r"(?=(G{3,}[ACGU]{1,7}G{3,}[ACGU]{1,7}G{3,}[ACGU]{1,7}G{3,}))")))
            self.pad = max(self.pad, 40)
        if p.avoid_cryptic_splice:
            self.rules.append(("splice_donor", re.compile(r"(?=([AC]AGGU[AG]AGU))")))
            self.rules.append(("splice_acceptor", re.compile(r"(?=([CU]{10,}[ACGU]AG))")))
            self.pad = max(self.pad, 24)
        self._ggg = re.compile(r"(?=(GGG))")

    def spans(self, seq: str) -> list[tuple[str, int, int]]:
        """All violations as (rule, start, length); starts are inside [0,n), matches may wrap."""
        n = len(seq)
        s = seq + seq[: self.pad]
        out = []
        for name, rx in self.rules:
            for m in rx.finditer(s):
                if m.start() < n:   # start inside the molecule; wrap handled by padding
                    out.append((name, m.start(), len(m.group(1))))
        return out

    def counts(self, seq: str) -> dict[str, int]:
        out: dict[str, int] = {}
        for name, _, _ in self.spans(seq):
            out[name] = out.get(name, 0) + 1
        return out

    def fixable(self, seq: str, free_mask) -> dict[str, int]:
        """Violations that overlap at least one designable position (the rest cannot be fixed)."""
        n = len(seq)
        out: dict[str, int] = {}
        for name, st, ln in self.spans(seq):
            if any(free_mask[(st + k) % n] for k in range(ln)):
                out[name] = out.get(name, 0) + 1
        return out

    def ggg(self, seq: str) -> int:
        return sum(1 for m in self._ggg.finditer(seq + seq[:2]) if m.start() < len(seq))


def gc_stats(seq: str, window: int, free_mask=None) -> dict:
    """Total GC plus the extreme GC of any circular window of `window` nt.

    The `*_free` variants look only at windows that contain at least one designable
    position. A window lying entirely inside a protected element (a GC-rich IRES, say)
    cannot be changed by any synonymous edit, so holding the design responsible for it
    would report an unfixable violation.
    """
    n = len(seq)
    gc = [1 if b in "GC" else 0 for b in seq]
    tot = sum(gc) / n
    w = min(window, n)
    ext = gc + gc[: w - 1]                       # circular windows
    cur = sum(ext[:w])
    lo = hi = cur
    free = None
    if free_mask is not None:
        fm = [1 if x else 0 for x in free_mask]
        fext = fm + fm[: w - 1]
        fcur = sum(fext[:w])
        lo_f = hi_f = (cur if fcur else None)
    for i in range(1, n):
        cur += ext[i + w - 1] - ext[i - 1]
        lo, hi = min(lo, cur), max(hi, cur)
        if free_mask is not None:
            fcur += fext[i + w - 1] - fext[i - 1]
            if fcur:
                lo_f = cur if lo_f is None else min(lo_f, cur)
                hi_f = cur if hi_f is None else max(hi_f, cur)
    out = {"gc": tot, "gc_win_min": lo / w, "gc_win_max": hi / w}
    if free_mask is not None:
        out["gc_win_min_free"] = (lo_f / w) if lo_f is not None else tot
        out["gc_win_max_free"] = (hi_f / w) if hi_f is not None else tot
    return out


def inverted_repeat_count(seq: str, k: int, free_mask) -> int:
    """# k-mers whose reverse complement occurs elsewhere (non-overlapping), >=1 occurrence in designable region."""
    pos: dict[str, list[int]] = {}
    n = len(seq)
    for i in range(n - k + 1):
        pos.setdefault(seq[i:i + k], []).append(i)
    cnt = 0
    for km, lst in pos.items():
        rc = revcomp(km)
        if rc < km or rc not in pos:   # count each unordered pair once (palindromes: rc==km handled below)
            continue
        for i in lst:
            for j in pos[rc]:
                if abs(i - j) >= k and (free_mask[i] or free_mask[j]):
                    cnt += 1
    return cnt
