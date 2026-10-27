# Epitope Filter Engine (local)

Local, offline prediction tools that produce IEDB-style tables for epitope filtering
(an Epitope Workbench ingests them). Each model is its own tool;
`engine_server.py` hosts them side by side and will become the combined pipeline.

| What | Status |
|---|---|
| **Pipeline** tab: one form, five models, chained, filters, per-peptide summary | working |
| NetMHCpan 4.2: binding affinity + eluted ligand (`MHC1/BA_EL_scripts/`) | working, runs in parallel |
| MixMHCpred 3.0: presentation (`MHC1/BA_EL_scripts/`) | working (reproduces the authors' test output) |
| PRIME 2.1: immunogenicity (`MHC1/Immunogenicity_scripts/`) | working (reproduces the authors' test output) |
| BigMHC EL + IM (`MHC1/Immunogenicity_scripts/`) | working |
| IEDB Class I immunogenicity (Calis 2013, "Positions to mask") as a baseline (`MHC1/Immunogenicity_scripts/calis_immunogenicity.py`) | working, built in, nothing to install |
| MHC class II BA/EL and CD4 immunogenicity | planned, see below |

Third-party predictors are licensed separately and are not stored here.

## The Pipeline tab

One form: paste a protein/FASTA (cut into fragments of the sizes you tick) or a list of
peptides, choose alleles (presets, search, or any allele NetMHCpan knows), tick the models,
set the filters, press **Run**.

- **Progress:** an overall bar with elapsed time and time left, a bar per model with what
  it is doing ("loading network 3/7", "scoring batch 1/2"), and a live log. **Cancel** stops
  the running programs. Reloading the page picks the job back up. One job runs at a time;
  others wait ("Waiting for the previous job...").
- **Chain:** NetMHCpan, MixMHCpred and PRIME are fast and score every pair. BigMHC is slow
  (it reloads its seven networks on every call), so with Chain on it scores only the pairs that
  pass the NetMHCpan / MixMHCpred filters. That is what keeps runs short.
- **A pair passes** when every model you ticked scored it and met its filter. A model that
  does not know an allele (MixMHCpred and PRIME have their own allele lists; BigMHC is human
  only) cannot pass it, and the page names the skipped alleles.
- **Output:** *By peptide* (for each peptide: how many of your alleles it passes and which,
  best NetMHCpan rank, best BigMHC scores) and *Every peptide-allele pair*. Rows that pass
  every filter are green. Both download as CSV (the page shows the best 5,000 rows; the CSV
  has all of them).
- **Before you run:** the form shows "4,901 peptides x 3 alleles -> 14,703 predictions", and
  a model that is not ready is refused immediately with the reason.
- BigMHC covers human HLA-I only; mouse alleles are skipped for it with a note. Runs over
  3,000,000 predictions are refused (use fewer alleles or sizes).
- Stray characters in a FASTA (for example a `"` left by GenBank formatting) make the
  windows containing them skip, and the page says so.

## Copy, paste, run

Put the folders for the models and this code side by side:

```
~/local_engines/
    runner_codes/     <- copy the CONTENTS of this directory here (from GitHub)
    netmhc/           <- the NetMHCpan tarball extracted anywhere below (even nested)
    bigmhc/           <- the BigMHC clone anywhere below
```

```bash
cd ~/local_engines/runner_codes
bash start.sh                     # http://127.0.0.1:5000
```

`start.sh` makes its own `.venv` (flask + pandas) on first run, so no conda env is
needed for the server. With no config file the server **discovers everything**:

| What | How it is found |
|---|---|
| NetMHCpan | a file named `netMHCpan` that has a `Linux_x86_64/` folder beside it, anywhere under `runner_codes/..` (up to 6 levels; shallowest wins) |
| BigMHC | a folder with `src/predict.py` and `models/` |
| BigMHC's Python | `BIGMHC_PYTHON` if set, else the conda env (named like `bigmhc` first) that can import torch, numpy, pandas and psutil |
| Scratch dir | `/tmp/epitope_engine` |

- **NetMHCpan's script is never edited.** DTU's script ships with a placeholder
  `NMHOME`; if it is wrong the server runs a patched copy from the scratch dir
  instead. (Editing it yourself is fine too.)
- To search somewhere else: `ENGINE_MODELS_DIR=/path bash start.sh`.
- The Settings tab shows "In use right now" (what was discovered) and the readiness
  checklist for each model.

## Installing MixMHCpred 3.0 and PRIME 2.1

Both are open GitHub repositories (free for academic use; for-profit use needs a licence from
the Ludwig Institute, see their READMEs). Put them beside the other tools and the engine finds them:

```bash
sudo apt install -y g++ perl                 # PRIME's C++ part and script
mkdir -p ~/local_engines/gfeller && cd ~/local_engines/gfeller
git clone --depth 1 https://github.com/GfellerLab/MixMHCpred.git
git clone --depth 1 https://github.com/GfellerLab/PRIME.git
```

MixMHCpred needs only numpy and pandas, which the engine's own environment already has. PRIME
ships a binary for other systems, so the engine rebuilds `lib/PRIME.x` for Linux the first time
PRIME runs. Settings -> **Test run** checks each one with GILGFVFTL on A*02:01. PRIME calls
MixMHCpred, so selecting PRIME needs MixMHCpred installed.

## Settings (optional override, edited live in the GUI)

`engine_config.json` is only needed to override discovery or change options such as
`device`. Copy `engine_config.example.json`; every path there defaults to `"auto"`.

- Settings shows the file in an editor. **Save** validates it (bad JSON or wrong
  types are rejected and the file is left untouched) and applies it on the next
  request, no restart.
- **Test run** predicts GILGFVFTL on A*02:01 end to end. The first BigMHC call loads
  the models from disk and can take minutes on a cold cache; later calls are quicker.
- `ENGINE_CONFIG=/path/file.json` points the server at a different file.

```json
"bigmhc": {"python": "~/miniconda3/envs/bigmhc/bin/python", "home": "~/local_engines/bigmhc/bigmhc"}
```

Each model runs as a subprocess with **its own interpreter**, so the server never
imports torch and BigMHC can live in a separate conda env.

## Layout

```
start.sh                    first-run venv + start
engine_server.py            server + Settings + health checks
common/config.py            engine_config.json load / validate / save; resolves 'auto'
common/discover.py          finds the models beside this folder; patched NetMHCpan copy
common/alleles.py           allele presets and spelling (HLA-A*02:01 vs HLA-A02:01)
MHC1/BA_EL_scripts/         NetMHCpan app (standalone or mounted)
MHC1/Immunogenicity_scripts/ BigMHC runner (CLI) + app
test_engine.py              config, server, both wrappers (stubs; optional real BigMHC)
```

Run the tests: `python3 test_engine.py`. With real BigMHC:
`BIGMHC_PYTHON=<env python> BIGMHC_HOME=<clone> python3 test_engine.py`.

## Known limits

- Jobs live in server memory: restarting the server forgets them (finished CSVs stay in
  `/tmp/epitope_engine/jobs`). Only one job runs at a time.

## IEDB Class I immunogenicity (baseline)

IEDB's "Class I pMHC Immunogenicity" model (Calis et al. 2013, *PLoS Comput Biol* 9:e1003266) is
tiny: a 20-value amino-acid scale times nine position weights, with some positions masked. It is
implemented natively here (`calis_immunogenicity.py`), so there is **nothing to download**.

- **Positions to mask**, as on the IEDB site: *1,2,C Terminal (Default)*, *By allele* (IEDB's
  per-allele table; an allele without an entry falls back to the default and the page says so) or
  *Custom* (1-based positions such as `2,5,9`). Default and Custom give the same score for every
  allele; By allele differs per allele.
- **Verified** two ways: it equals IEDB's own `immunogenicity_predictor.py` (next-generation T Cell
  Class I package v0.2-beta) exactly on 28,000 random peptides over all lengths and mask choices,
  and it reproduces all 11,781 peptide scores in the HEV export produced by IEDB's website
  (maximum difference 0.0). `test_engine.py` keeps reference values from both.
- Higher means more immunogenic. IEDB defines no cut-off, so the filter is off by default.
- The constants come from IEDB's package (Non-Profit Open Software License 3.0, La Jolla Institute for
  Immunology); the code here is an independent implementation of the published model.

IEDB's full standalone bundle (<https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/>, 743 MB,
`uv` + Python 3.11) also contains NetMHCpan 4.1/4.2, MHCflurry, NetChop and NetCTL for antigen
processing. It is not needed for anything above; it is the place to look if you later want those.

## MHC class II: official sources found

Not built yet; these are the current official tools, from their own pages:

| Purpose | Tool | Where | Notes |
|---|---|---|---|
| BA + EL | **NetMHCIIpan 4.3** (package 4.3k) | <https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/> | DTU, same kind of download as NetMHCpan; free for academic use. HLA-DR/DQ/DP plus mouse H-2 and bovine BoLA-II. |
| EL (second opinion) | **MixMHC2pred 2.0** | <https://github.com/GfellerLab/MixMHC2pred> (releases page, `MixMHC2pred-2.0.zip`) | Open GitHub release; web version at <https://mixmhc2pred.gfellerlab.org>. |
| CD4 immunogenicity | **CD4episcore** | <http://tools.iedb.org/CD4episcore> | IEDB; no HLA typing needed. IEDB release 3.0 (March 2025) lists a class II T-cell standalone; whether it contains CD4episcore is unconfirmed. |

IEDB's own class II bundle also wraps NetMHCIIpan 4.2/4.3.
