# Install (WSL / Ubuntu / Linux)

Everything below was run and verified on Linux x86-64 with Python 3.13 and ViennaRNA 2.7.2.
DERNA is a compiled C++ binary and only ships through bioconda, so it needs conda; the rest is pip.

## 1. Conda (skip if you already have conda/mamba)

Miniforge is the conda-forge installer and already has the right channel defaults:

```bash
cd ~
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash          # so `conda activate` works in new shells
```

## 2. The environment

```bash
conda create -y -n circopt python=3.11
conda activate circopt

# DERNA (C++ binary, bioconda). conda-forge must come first or the solve picks wrong deps.
conda install -y -c conda-forge -c bioconda derna

# Python side
pip install -r requirements.txt
```

`requirements.txt` pulls: `viennarna` (the `RNA` Python module), `flask`, `numpy`, `pandas`,
`matplotlib`, `python_codon_tables` (Kazusa human + mouse tables), `pytest`.

If you prefer mamba/micromamba (faster solve), this is the exact command used here:

```bash
micromamba create -y -n circopt -c conda-forge -c bioconda derna
```

## 3. Verify

```bash
conda activate circopt
python -c "import RNA; print('ViennaRNA', RNA.__version__)"
python -c "import python_codon_tables as p; print(len(p.get_codons_table('m_musculus_10090')), 'aa in mouse table')"
derna | head -3                       # prints its usage block
python -m circopt.cli selftest        # folds a short circle + checks the scorer
```

Expected: `ViennaRNA 2.7.2`, `21 aa in mouse table`, a DERNA usage block, and `selftest OK`.

## 4. Run

```bash
# one-shot CLI run
python -m circopt.cli run examples/example_params.json --out runs/demo

# web runner (Flask)
python -m circopt.app            # then open http://127.0.0.1:5000
```

## Notes / gotchas found while building this

* **DERNA's `-c` custom codon table is unusable.** Six plausible CSV layouts were tried and the
  binary aborts with `std::invalid_argument` ("no char found" / "stof") on all of them; the format
  is undocumented on GitLink and the source is not mirrored on GitHub. So DERNA runs with its
  built-in (human) table and is used only as a *seed/reference* generator. Human and mouse codon
  usage are nearly identical (e.g. Leu CUG 0.40 vs 0.39), and the dual-species balancing is done
  by this package's own optimiser, so nothing is lost.
* **DERNA optimises the opposite structure direction** to the default here: it minimises MFE
  (maximally structured, the LinearDesign direction). It is therefore run as a comparison arm and
  as a high-CAI seed, not as the optimiser. Set `structure_goal: "minimize_mfe"` to test that
  assumption instead.
* **DERNA is O(n^3) and single-threaded**: ~27 s for a 70 aa protein at fixed lambda, and the
  lambda sweep (`-s 3`) is several minutes for the same input. It is off by default for long CDSs;
  see `derna_max_aa` in the params.
* **ViennaRNA 2.7.2 `subopt()` segfaults with `circ=1`** -- this package never calls it in
  circular mode. Ensemble structures come from `pbacktrack()`, which needs `uniq_ML=1` or it
  silently returns nothing.
