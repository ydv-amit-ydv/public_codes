# MHC-I immunogenicity

Two models, run separately for now:

| Model | Source | Status |
|---|---|---|
| **BigMHC EL + IM** (Albert et al. 2023) | KarchinLab, GitHub | wrapper written: `bigmhc_runner.py` |
| **IEDB Class I pMHC Immunogenicity** (Calis et al. 2013; has the "Positions to mask" option) | IEDB standalone | install steps below; wrapper waits for its `--help` |

Neither is stored in this repo; install them elsewhere and point the scripts at them.

## Download links

**BigMHC**
- Code: <https://github.com/KarchinLab/bigmhc> (`git clone https://github.com/karchinlab/bigmhc.git`, about 5 GB because the trained models are in the repo)
- Paper: <https://www.nature.com/articles/s42256-023-00694-6>
- Training/benchmark data: <https://doi.org/10.17632/dvmz6pkzvb>

**IEDB Class I immunogenicity**
- Download page (follow its licence prompt): <https://tools.iedb.org/immunogenicity/download/>
- Package README (v3.0, Python 3.6+): <https://downloads.iedb.org/tools/immunogenicity/LATEST/README>
- Newer bundle that also covers binding, elution and processing: <https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/README>
- Web version: <https://nextgen-tools.iedb.org>

Both IEDB packages are free for academic/non-profit use; the legacy one asks you to
accept the NetMHC 3.0 licence and the Non-Profit OSL 3.0. Industry use needs a
licence (license@iedb.org).

## Install BigMHC (WSL / Ubuntu, conda)

BigMHC needs torch, so give it its own env and let `engine_config.json` point the
server at it (the server never imports torch).

```bash
conda create -n bigmhc python=3.9 -y        # the authors' version; see "Tested" below
conda activate bigmhc
pip install numpy pandas psutil
pip install torch --index-url https://download.pytorch.org/whl/cpu    # no GPU needed

mkdir -p ~/local_engines/bigmhc && cd ~/local_engines/bigmhc
git clone --depth 1 https://github.com/karchinlab/bigmhc.git
cd bigmhc/src
python predict.py -i=../data/example1.csv -m=el -t=2 -d="cpu"        # README's own example
which python                                  # <- goes in models.bigmhc.python
```

The engine finds this clone and the env by itself when it sits beside `runner_codes/`
(the env needs torch, numpy, pandas, psutil). Open Settings and press **Test run**; set
`models.bigmhc.python` / `home` only to override.

**Tested:** the authors' example reproduces their reference file (max difference
5e-7) on Linux x86_64 with CPU-only Python 3.13 and torch 2.14, so the pinned
Python 3.9 / torch 1.13 stack is not required. 3.9 remains the authors' tested
configuration; a 3.12 env is likely fine but I have not run it.

## Run

From the engine GUI: **BigMHC EL/IM** tab (paste FASTA or peptides, pick alleles, Run).

From the command line (the interpreter and install come from `engine_config.json`,
or pass `--python` / `--bigmhc-home`):

```bash
python3 bigmhc_runner.py --fasta ../../../HEV/Raw/AF444002.1.fa --lengths 9 10 \
    --alleles HLA-A02:01 HLA-B07:02 -o bigmhc_results.csv
```

Output is one row per peptide x allele with `BigMHC_EL` and `BigMHC_IM`; alleles are
written `HLA-A*02:01`. Mouse alleles are not offered: BigMHC's allele panel is
human HLA-I.

`--workbench-csv FILE` also writes the IM score in the column layout of
`IEDB/HEV/Raw/MHC_1_Immunogenicity_*.csv`. BigMHC IM and IEDB's Calis score are
different scales, so keep them in separate folders rather than mixing them in one
Workbench `Raw/` directory.

## Status and caveats

- Verified against real BigMHC: the `.prd` output is `mhc,pep,[tgt,]len,BigMHC_EL|IM`
  and the runner keeps input order (it checks the echoed peptide column and refuses
  to align by position if it ever differs). GILGFVFTL on A*02:01 scores EL 0.48, IM 0.76.
- BigMHC IM is trained on immunogenicity assays; treat the score as a ranking
  signal, not a probability.
- The IEDB tool's "Positions to mask" setting hides selected peptide positions
  (the MHC anchor residues by default, as I understand its description; check the
  package README) from the scoring. Its wrapper will expose that as an option once
  the package is installed and I can see its `--help` and output.
