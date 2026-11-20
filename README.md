# public_codes

## epitope_filter_engine

*(renamed from `local_engine_codes`; same code, plus an install guide)*

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

Quick start on Linux / WSL:

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git ~/public_codes
```

```bash
mkdir -p ~/local_engines && cd ~/public_codes/epitope_filter_engine && ENGINE_MODELS_DIR=~/local_engines bash start.sh
```

That already works (<http://127.0.0.1:5000>) with the built-in IEDB/Calis model. To add the
predictors, run `bash install_tools.sh` for the open-licence ones and follow
[`epitope_filter_engine/INSTALL.md`](epitope_filter_engine/INSTALL.md) for NetMHCpan and
the IEDB bundles.

- **[`INSTALL.md`](epitope_filter_engine/INSTALL.md)** — full install, every tool, verification, troubleshooting
- **[`README.md`](epitope_filter_engine/README.md)** — what each model gives and how the consensus works

Update with `git -C ~/public_codes pull`; your runs and settings live outside the code folder.

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

## OpenPrism (moved)

**OpenPrism** — the open-source scientific graphing and statistics workbench — now has its own repository:
**https://github.com/ydv-amit-ydv/openprism**

```bash
git clone https://github.com/ydv-amit-ydv/openprism.git
cd openprism
python3 start.py --new --new --new --newserveruserdata ~/openprism-data
```
