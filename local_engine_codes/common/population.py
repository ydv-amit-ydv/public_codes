"""Population allele frequencies and HLA coverage (Hardy-Weinberg), read from a weightage CSV.

The CSV is the output of the Epitope Workbench's hla_allele_weightage.py (one row per allele): `Allele` (A*11:01),
`Locus`, one `f_<region>` column per region, and `PopWeightedFreq` (the national frequency). Region weights are not
in the CSV, but PopWeightedFreq = sum_r w_r * f_r holds exactly, so they are recovered by least squares (checked
against the workbook's own Weights_Used sheet to 16 digits). Nothing else is needed.

How an allele of any HLA name is mapped to the data (the same rules as the Epitope Workbench):
  exact            the allele is typed (class I A/B/C, DRB1, DQB1)
  dqb1_component   DQA1*x/DQB1*y heterodimer: the DQB1 chain's frequency stands in for the pair (DQA1 is not typed)
  ld_proxy         DRB3/4/5: group-level estimate summed from the DRB1 alleles they are linked to
  none             no data (every HLA-DP pair, mouse alleles, anything not typed): weight 0, excluded from coverage

Coverage of a set of alleles = fraction of people carrying at least one of them. Within a locus the frequencies of
the alleles add; a person has two haplotypes, so P(carries none) = (1 - sum)^2 (Hardy-Weinberg); loci are treated as
independent (linkage disequilibrium between loci is ignored, which the Workbench documents as an approximation).
This is done per region and averaged with the region weights.
"""
import os
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

DRB345_LINKAGE = {
    "DRB3": ("DRB1*03", "DRB1*11", "DRB1*12", "DRB1*13", "DRB1*14"),
    "DRB4": ("DRB1*04", "DRB1*07", "DRB1*09"),
    "DRB5": ("DRB1*15", "DRB1*16"),
}


@dataclass
class Resolved:
    allele: str
    key: str                  # frequency-table entry it maps to ("A*02:01", "DQB1*02:01", "DRB3" ...)
    group: str                # locus group coverage is pooled within ("A", "DQ", "DRB345", "DP", ...)
    source: str               # exact | dqb1_component | ld_proxy | none
    freqs: dict = field(default_factory=dict)       # region -> frequency
    national: float = 0.0
    note: str = ""


def strip_hla(a):
    return a.replace("HLA-", "").strip()


def group_of(a):
    bare = strip_hla(a)
    if "/" in bare:
        return "DQ" if bare.startswith("DQA") else "DP"
    fam = bare.split("*")[0]
    return "DRB345" if fam in DRB345_LINKAGE else fam


class PopulationModel:
    def __init__(self, table, path="", label="India"):
        self.path, self.label = path, label
        self.regions = [c[2:] for c in table.columns if c.startswith("f_")]
        F = table[[f"f_{r}" for r in self.regions]].to_numpy(float)
        y = table["PopWeightedFreq"].to_numpy(float)
        w = np.linalg.lstsq(F, y, rcond=None)[0]
        self.fit_error = float(np.abs(F @ w - y).max())
        if self.fit_error > 1e-6 or w.min() < -1e-9:
            raise ValueError(f"{path}: PopWeightedFreq is not a weighted sum of the f_<region> columns "
                             f"(residual {self.fit_error:.2g}); is this a weightage CSV from hla_allele_weightage.py?")
        w = np.clip(w, 0, None)
        self.region_weights = dict(zip(self.regions, w / w.sum()))
        self.table = table
        self.freqs = {a: dict(zip(self.regions, row)) for a, row in zip(table["Allele"], F)}
        self.national = dict(zip(table["Allele"], y))
        self.locus = dict(zip(table["Allele"], table["Locus"]))
        self._cache = {}

    @classmethod
    def from_csv(cls, path, label="India"):
        t = pd.read_csv(path)
        need = {"Allele", "Locus", "PopWeightedFreq"}
        if need - set(t.columns) or not any(c.startswith("f_") for c in t.columns):
            raise ValueError(f"{path}: needs columns Allele, Locus, PopWeightedFreq and f_<region> (see common/population.py)")
        return cls(t, os.fspath(path), label)

    # ---- mapping an allele name to the data ------------------------------------------------------------
    def resolve(self, allele):
        bare = strip_hla(allele)
        grp = group_of(allele)
        if bare in self.freqs:
            return Resolved(allele, bare, grp, "exact", self.freqs[bare], self.national[bare])
        if "/" in bare:
            dqb = next((p for p in bare.split("/") if p.startswith("DQB1")), None)
            if dqb and dqb in self.freqs:
                return Resolved(allele, dqb, grp, "dqb1_component", self.freqs[dqb], self.national[dqb],
                                f"frequency of the {dqb} chain; DQA1 is not typed")
            return Resolved(allele, bare, grp, "none", note="no Indian frequency data for HLA-DP" if grp == "DP" else "not typed")
        fam = bare.split("*")[0]
        if fam in DRB345_LINKAGE:
            pref = DRB345_LINKAGE[fam]
            summed = {r: 0.0 for r in self.regions}
            for name, fr in self.freqs.items():
                if name.startswith(pref):
                    for r in self.regions:
                        summed[r] += fr[r]
            nat = sum(self.region_weights[r] * min(summed[r], 1.0) for r in self.regions)
            return Resolved(allele, fam, grp, "ld_proxy", summed, nat,
                            f"group-level estimate from linked {', '.join(pref)}; shared by all {fam} alleles")
        return Resolved(allele, bare, grp, "none", note="not typed in this dataset")

    def allele_info(self, allele):
        """(national frequency, source) for display; frequency 0 when there is no data."""
        r = self.resolve(allele)
        return r.national, r.source

    # ---- coverage of a set of alleles -------------------------------------------------------------------------
    def per_region(self, alleles):
        key = frozenset(alleles)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        groups = {}
        for a in alleles:
            r = self.resolve(a)
            if r.source != "none":
                groups.setdefault(r.group, {})[r.key] = r         # one entry per distinct table key (no double counting)
        out, national = {}, 0.0
        for region, w in self.region_weights.items():
            uncovered = 1.0
            for grp, members in groups.items():
                vals = [m.freqs.get(region, 0.0) for m in members.values()]
                if all(m.source == "ld_proxy" for m in members.values()):
                    pooled = max(vals)                              # DRB3/4/5 share a group-level estimate
                else:
                    pooled = sum(vals)
                uncovered *= (1.0 - min(pooled, 1.0)) ** 2
            out[region] = 1.0 - uncovered
            national += w * out[region]
        out["national"] = national
        if len(self._cache) < 200_000:
            self._cache[key] = out
        return out

    def coverage(self, alleles):
        return self.per_region(alleles)["national"]

    # ---- how representative is the user's allele set -----------------------------------------------------------------
    def representativeness(self, alleles):
        """Per locus group: share of Indian haplotypes the selected alleles account for, and the best alleles not selected."""
        sel = {}
        for a in alleles:
            r = self.resolve(a)
            if r.source in ("exact", "dqb1_component"):
                sel.setdefault(r.group, {})[r.key] = r.national
        rows = []
        for grp, members in sorted(sel.items()):
            rows.append({"locus": grp, "alleles": len(members), "haplotype_share": min(1.0, sum(members.values()))})
        return rows

    def missing_top(self, alleles, loci, per_locus=3):
        have = {self.resolve(a).key for a in alleles}
        out = []
        for loc in loci:
            cand = sorted(((a, f) for a, f in self.national.items() if self.locus.get(a) == loc and a not in have),
                          key=lambda x: -x[1])[:per_locus]
            out += [{"allele": ("HLA-" + a) if loc in ("A", "B", "C", "DRB1") else a, "locus": loc, "frequency": f} for a, f in cand]
        return out

    def alleles_for_coverage(self, locus, target=0.9):
        """Smallest set of the most frequent alleles of a locus whose frequencies add up to `target`."""
        items = sorted(((a, f) for a, f in self.national.items() if self.locus.get(a) == locus), key=lambda x: -x[1])
        out, acc = [], 0.0
        for a, f in items:
            out.append(a)
            acc += f
            if acc >= target:
                break
        return out

    def info(self):
        return {"label": self.label, "path": self.path, "regions": self.region_weights, "alleles": len(self.freqs),
                "loci": sorted(set(self.locus.values())), "fit_error": self.fit_error}


def find_csv(root):
    """A weightage CSV: ENGINE_POPULATION_CSV, else a *hla_allele_weightage*.csv beside/below the models folder."""
    env = os.environ.get("ENGINE_POPULATION_CSV")
    if env and os.path.isfile(os.path.expanduser(env)):
        return os.path.expanduser(env)
    import glob
    for pat in ("*hla_allele_weightage*.csv", "*/*hla_allele_weightage*.csv", "*/*/*hla_allele_weightage*.csv"):
        hits = sorted(glob.glob(os.path.join(root, pat)))
        if hits:
            return hits[0]
    return None


def load(cfg, root):
    """-> (PopulationModel | None, message). cfg is the config's `population` object."""
    cfg = cfg or {}
    path = cfg.get("weights_csv", "auto")
    path = find_csv(root) if path in ("", "auto", None) else os.path.expanduser(path)
    if not path or not os.path.isfile(path):
        return None, ("No population file found: coverage is counted as the share of the chosen alleles. "
                      "Set population.weights_csv in Settings to an HLA allele weightage CSV.")
    try:
        return PopulationModel.from_csv(path, cfg.get("label", "India")), ""
    except Exception as e:                                   # noqa: BLE001 - shown to the user, run continues
        return None, f"Population file not used: {e}"


def presets(pop):
    """Allele sets built from the population file: the most common alleles per locus up to 80 % of haplotypes.
    -> (class I preset or None, class II preset or None), each {"alleles": [...], "note": str}."""
    if pop is None:
        return None, None
    one, loci = [], []
    for loc in ("A", "B", "C"):
        got = pop.alleles_for_coverage(loc, 0.8)
        one += ["HLA-" + a for a in got]
        loci.append(f"{loc} {len(got)}")
    two = ["HLA-" + a for a in pop.alleles_for_coverage("DRB1", 0.8)]
    c1 = {"alleles": one, "note": f"{pop.label}: the most common alleles covering 80 % of each locus ({', '.join(loci)})"}
    c2 = {"alleles": two, "note": f"{pop.label}: DRB1 alleles covering 80 % of the locus (DQ and DP are not typed in the file)"}
    return c1, c2
