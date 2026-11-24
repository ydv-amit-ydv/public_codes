# public_codes

## epitope_filter_engine

A local, offline web GUI that runs **eleven MHC epitope predictors** side by side and fuses
them into one ranked peptide list — both T-cell arms in one job:

- **MHC class I (CTL/CD8):** NetMHCpan 4.2, MixMHCpred 3.0, MHCflurry 2, PRIME 2.1,
  BigMHC EL + IM, NetCTLpan 1.1, NetChop 3.1, IEDB/Calis immunogenicity (built in)
- **MHC class II (HTL/CD4):** NetMHCIIpan 4.3/4.2, MixMHC2pred 2.1, CD4episcore,
  ImmuScope-IM, TLimmuno2

Each tool votes per peptide-allele pair; the engine reports agreement, population-weighted
coverage (HLA frequencies), a priority score, CTL+HTL regions, a host-proteome self-check,
a minimum allele panel, and a self-contained HTML report. It finds the models on disk
automatically and runs each one in its own interpreter.

**The predictors are not included and not redistributed here** — each has its own licence.

Quick start on Linux / WSL — clone, one command, done:

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git
```

```bash
cd public_codes/epitope_filter_engine && bash start.sh
```

The first run installs everything it can and opens the engine on <http://localhost:5000>;
later runs just open it. The one exception is **NetMHCpan 4.2**, which needs DTU's licence
form — the engine runs without it and tells you where to put it.

Full details, and every step by hand: [`epitope_filter_engine/INSTALL.md`](epitope_filter_engine/INSTALL.md).
What each model gives and how the consensus works: [`epitope_filter_engine/README.md`](epitope_filter_engine/README.md).

## circrna_opt

`circopt` — a dual-species, structure-aware CDS optimiser for **circRNA vaccines**. It
redesigns the coding sequence of a circular RNA using synonymous codon changes only, so the
protein never changes, while balancing **mouse and human** codon adaptation against secondary
structure folded in **circular mode** (ViennaRNA `--circ`). IRES, UTRs and spacers are declared
as protected regions and are never edited.

Plain Flask + HTML (no Streamlit), a CLI, and 28 unit tests.

Quick start on Linux / WSL:

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git
cd public_codes/circrna_opt

conda create -y -n circopt python=3.11 && conda activate circopt
conda install -y -c conda-forge -c bioconda derna      # optional seed/comparison arm
pip install -r requirements.txt

python -m circopt.cli selftest
python -m circopt.app                                   # http://localhost:5000
```

Full instructions: [`circrna_opt/README.md`](circrna_opt/README.md) and
[`circrna_opt/INSTALL.md`](circrna_opt/INSTALL.md).

## home_designer

**Home Designer**: our own floor-plan and interior design tool. You draw walls (exterior and interior, typed
exact lengths) and rooms, furnish them from a 337-item catalog, annotate places with notes, and see it all in 3D with real
photo textures, ambient occlusion and real CC0 3D models from Poly Haven. You can also trace over a sketch or photo.
Designs are saved as JSON in `home_designer/designs/` so they're versioned in git. Only Python 3 is needed.

```bash
cd public_codes/home_designer && bash start.sh      # http://127.0.0.1:8765
```

Details: [`home_designer/README.md`](home_designer/README.md).

## OpenPrism (moved)

**OpenPrism** — the open-source scientific graphing and statistics workbench — now has its own repository:
**https://github.com/ydv-amit-ydv/openprism**

```bash
git clone https://github.com/ydv-amit-ydv/openprism.git
cd openprism
python3 start.py --new --new --new --newserveruserdata ~/openprism-data
```
