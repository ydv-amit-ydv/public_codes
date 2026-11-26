# Where the bundled data came from

This repository ships a gzipped copy of every dataset under `data/`, about
8.6 MB in total, so a clone can train and validate everything **with no
network access at all**. `dataio.py` opens whichever copy is present and
prefers a freshly harvested `.csv` over the bundled `.csv.gz`.

Each dataset is redistributed here in the form the harvest script produced:
selected columns, de-duplicated, and in the case of the eluted ligands with
decoys added. None of it is modified beyond that. Every one of these sources
is openly available and asks to be cited rather than licensed, but **they are
other people's databases, and if you redistribute this repository onward you
are redistributing them too** — check each source's current terms for your
use, particularly if that use is commercial.

| File | Rows | What it is | Source | Harvest script |
|---|---|---|---|---|
| `binding_affinity.csv.gz` | 67,365 | Measured peptide–MHC binding affinities (IC50, nM) | [IEDB](https://www.iedb.org) | `harvest_iedb.py` |
| `tcell_immunogenicity.csv.gz` | 60,820 | T cell assay outcomes for peptide–MHC pairs | [IEDB](https://www.iedb.org) | `harvest_iedb.py` |
| `eluted_ligands.csv.gz` | 852,820 | Mass-spectrometry eluted ligands + proteome-sampled decoys | [IEDB](https://www.iedb.org) elution assays; decoy parent proteins from [UniProt](https://www.uniprot.org) | `harvest_eluted.py` |
| `cmipb_ab_titer.csv.gz` | 52,576 | Antibody titers, 166 adults given a Tdap boost | [CMI-PB](https://www.cmi-pb.org) | `harvest_cmipb.py` |
| `cmipb_cell_frequency.csv.gz` | 25,047 | PBMC cell frequencies for the same subjects | [CMI-PB](https://www.cmi-pb.org) | `harvest_cmipb.py` |
| `trial_foldrise.csv.gz` | 16,940 | Posted antibody fold-rise results from registered vaccine trials | [ClinicalTrials.gov](https://clinicaltrials.gov) API v2 | `harvest_trials.py` |

## Please cite

**IEDB** — Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database
(IEDB): 2018 update. *Nucleic Acids Research* 2019;47(D1):D339–D343.

**UniProt** — The UniProt Consortium. UniProt: the Universal Protein
Knowledgebase in 2023. *Nucleic Acids Research* 2023;51(D1):D523–D531.

**CMI-PB** — The Computational Models of Immunity – Pertussis Boost
consortium, `cmi-pb.org`.

**ClinicalTrials.gov** — U.S. National Library of Medicine. Individual trials
carry their own NCT identifiers, which are in the `nct_id` column.

## Refreshing it

The Setup page has a **Re-download every dataset** button, or:

```bash
python3 harvest_iedb.py && python3 harvest_cmipb.py && python3 harvest_eluted.py && python3 harvest_trials.py
```

That writes plain `.csv` files, which `dataio.py` then prefers over the
bundle. Note that **every number in VALIDATION.md was produced from the
bundled copies**: a fresh harvest returns whatever those databases hold today,
so results may move. If they move a lot, that is worth knowing about and is
not necessarily a bug in this code.

To go back to the bundled copies, delete the plain CSVs:

```bash
rm -f data/*.csv
```
