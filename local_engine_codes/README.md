# Epitope Filter Engine (local)

Local, offline prediction tools that produce IEDB-style tables for epitope filtering
(an Epitope Workbench ingests them). Each model is its own tool;
`engine_server.py` hosts them side by side and will become the combined pipeline.

| Directory | Tool | Status |
|---|---|---|
| `MHC1/BA_EL_scripts/` | NetMHCpan 4.2 binding affinity / eluted ligand | working (tested against real NetMHCpan output) |
| `MHC1/Immunogenicity_scripts/` | BigMHC EL + IM | working (tested against real BigMHC) |
| `MHC1/Immunogenicity_scripts/` | IEDB Class I immunogenicity (with "Positions to mask") | not installed yet |
| `MHC2/...` | MHC-II BA/EL and CD4 immunogenicity | planned |

Third-party predictors are licensed separately and are not stored here.

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

- Requests are synchronous. A large BigMHC scan on CPU can run for minutes; async
  jobs with progress polling are the next step before full-protein scans.
- The IEDB immunogenicity tool and MHC-II models are not wired in yet; the config
  already has an `iedb_immunogenicity` slot.
