#!/usr/bin/env python3
"""Where the datasets live, and how to open them.

The repository ships every dataset gzipped under `data/`, so a clone can
train and validate with no network at all. A harvest script writes the plain
`.csv`; if one is present it wins, because it is by definition newer than the
bundled copy. Otherwise the bundled `.csv.gz` is used.

pandas reads `.gz` transparently, so this module is only about *which* path
to open and about saying clearly which copy you got.
"""
from pathlib import Path

import pandas as pd

DATA = Path(__file__).resolve().parent / "data"

# name -> (what it is, where it came from, what to cite)
SOURCES = {
    "binding_affinity": (
        "Measured peptide-MHC binding affinities",
        "IEDB (iedb.org)",
        "Vita R et al. The Immune Epitope Database (IEDB): 2018 update. "
        "Nucleic Acids Res. 2019;47(D1):D339-D343."),
    "tcell_immunogenicity": (
        "T cell assay outcomes for peptide-MHC pairs",
        "IEDB (iedb.org)",
        "Vita R et al. Nucleic Acids Res. 2019;47(D1):D339-D343."),
    "eluted_ligands": (
        "Mass-spectrometry eluted ligands plus proteome-sampled decoys",
        "IEDB elution assays; decoy parents from UniProt",
        "Vita R et al. Nucleic Acids Res. 2019;47(D1):D339-D343. "
        "UniProt Consortium. Nucleic Acids Res. 2023;51(D1):D523-D531."),
    "cmipb_ab_titer": (
        "Antibody titers from 166 adults given a Tdap boost",
        "CMI-PB (cmi-pb.org)",
        "CMI-PB Consortium. The Computational Models of Immunity - "
        "Pertussis Boost resource."),
    "cmipb_cell_frequency": (
        "PBMC cell frequencies for the same subjects",
        "CMI-PB (cmi-pb.org)",
        "CMI-PB Consortium."),
    "trial_foldrise": (
        "Posted antibody fold-rise results from registered vaccine trials",
        "ClinicalTrials.gov API v2",
        "ClinicalTrials.gov, National Library of Medicine."),
}


def data_path(name):
    """Path to a dataset, preferring a freshly harvested CSV over the bundle.

    `name` is the stem, e.g. "binding_affinity". Returns None if neither a
    harvested nor a bundled copy exists.
    """
    stem = name[:-4] if name.endswith(".csv") else name
    plain = DATA / f"{stem}.csv"
    packed = DATA / f"{stem}.csv.gz"
    if plain.exists():
        return plain
    if packed.exists():
        return packed
    return None


def have(name):
    return data_path(name) is not None


def is_bundled(name):
    p = data_path(name)
    return p is not None and p.suffix == ".gz"


def read(name, **kw):
    """Read a dataset by stem, from whichever copy is present."""
    p = data_path(name)
    if p is None:
        raise FileNotFoundError(
            f"{name}: no data/{name}.csv and no data/{name}.csv.gz. "
            f"Run the harvest for it, or restore the bundled copy "
            f"with `git checkout data/`.")
    return pd.read_csv(p, **kw)


def inventory():
    """One row per dataset: present, bundled or harvested, and its source."""
    out = []
    for name, (what, where, cite) in SOURCES.items():
        p = data_path(name)
        out.append({
            "name": name, "what": what, "source": where, "citation": cite,
            "present": p is not None,
            "copy": None if p is None else ("bundled" if p.suffix == ".gz"
                                            else "harvested"),
            "path": None if p is None else str(p.relative_to(DATA.parent)),
            "bytes": None if p is None else p.stat().st_size,
        })
    return out


if __name__ == "__main__":
    for row in inventory():
        state = row["copy"] or "MISSING"
        size = f"{row['bytes'] / 1e6:.1f} MB" if row["bytes"] else "-"
        print(f"  {row['name']:24s} {state:10s} {size:>9s}  {row['source']}")
