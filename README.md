# public_codes

## local_engine_codes

A local web GUI that runs MHC-I epitope predictors on your own machine:
**NetMHCpan 4.2** (binding affinity / eluted ligand) and **BigMHC** (presentation +
immunogenicity), with a Settings tab. It finds the models automatically.

The predictors themselves are not included: NetMHCpan is DTU academic-licence
software (download it from DTU), BigMHC is a separate open-source repository.

Quick start on Linux / WSL:

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git
mkdir -p ~/local_engines/runner_codes
cp -r public_codes/local_engine_codes/. ~/local_engines/runner_codes/
# put the NetMHCpan extraction in ~/local_engines/netmhc/ and the BigMHC clone in ~/local_engines/bigmhc/
cd ~/local_engines/runner_codes && bash start.sh      # http://localhost:5000
```

Full instructions: [`local_engine_codes/README.md`](local_engine_codes/README.md).
Update later with `git -C ~/public_codes pull` (or re-clone) and copy again.

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
