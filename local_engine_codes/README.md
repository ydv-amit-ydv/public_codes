# Epitope Filter Engine (local)

Local, offline prediction tools that produce the IEDB-style tables the Epitope
Workbench (`IEDB/HEV/Scripts/epitope_app.py`) filters. Each model is its own tool;
`engine_server.py` hosts them side by side and will become the combined pipeline.

| What | Status |
|---|---|
| **Pipeline** tab: one form, nine models, chained, filters, per-peptide summary | working |
| NetMHCpan 4.2: binding affinity + eluted ligand (`MHC1/BA_EL_scripts/`) | working, runs in parallel |
| MixMHCpred 3.0: presentation (`MHC1/BA_EL_scripts/`) | working (reproduces the authors' test output) |
| MHCflurry 2: affinity + processing + presentation (`MHC1/BA_EL_scripts/`) | working (real model checked) |
| PRIME 2.1: immunogenicity (`MHC1/Immunogenicity_scripts/`) | working (reproduces the authors' test output) |
| IEDB Class I immunogenicity (Calis 2013, "Positions to mask") baseline | working, built in |
| BigMHC EL + IM (`MHC1/Immunogenicity_scripts/`) | working |
| NetCTLpan 1.1 + NetChop 3.1 antigen processing, via IEDB's bundle (`MHC1/Processing_scripts/`) | working (real tools checked) |
| MHC class II: NetMHCIIpan EL + BA (4.3 / 4.2 / IEDB default), MixMHC2pred 2.1, CD4episcore (`MHC2/`) | working (real tools checked) |

Third-party predictors are licensed separately and are not stored here.

## The Pipeline tab: both T-cell arms in one run

One form. Paste a protein / FASTA (or a list of peptides), and the page runs **MHC class I (CTL, CD8)** and **MHC class II (HTL, CD4)**
on it in one job. Each arm has its own fragment sizes (all the sizes the models accept are ticked: 8-14 for class I, 11-25 for class II),
alleles (the India population preset comes first when the population file is found), models and chain switch; an arm can be switched off.
A small **?** beside every parameter explains it and gives its default.

1. **Windows are cut once and de-duplicated.** Identical windows (the same sequence at several positions or in several proteins, or the same
   peptide from two sizes of a loose list) are predicted once and mapped back; the preview shows how many were merged.
2. **Each arm** runs its models, then the consensus step (see below): votes, agreement, population-weighted coverage, priority score.
3. **Regions.** Candidates of both classes are grouped by position into labelled regions: **CTL+HTL** (a class I and a class II candidate
   overlap), **CTL only**, **HTL only**; CTL+HTL regions rank first. Each peptide table carries its `Region_ID`.
4. **Outputs.** The page shows the regions, the candidates and every value per class, the funnel and suggestions, and a downloads tab:
   `report.html` (one self-contained summary page with a map of the protein), `regions.csv`, per class the master CSV (every tool value, vote,
   agreement count and consensus for every peptide x allele), peptides, candidates, minimum panel, run details, one CSV per model, and one zip
   of everything.
5. **Automation.** *Save settings* writes the form as json; `python3 run_engine.py epitope_settings.json -o results/` runs the same thing
   without the page (and `--fasta other.fa` swaps the protein); the page's API takes the same json (`POST /api/pipeline/run`).

How it stays fast (measured on a 4-core machine, HEV genome 2,486 aa, class I sizes 8-14 x 27 alleles = 462,000 predictions):
MixMHCpred 9 s, PRIME 7 s. BigMHC costs about 26 s to start plus 9 ms per pair, for each of EL and IM, so it is chained: EL scores only pairs
that can still become binders (the gate is exact for the agreement rule, so no result changes), and IM scores only pairs already called
binders. Class II chains CD4episcore the same way.

- **Progress:** an overall bar with elapsed time and time left, a bar per model, and a live log. **Cancel** stops the running programs.
  Reloading the page picks the job back up. One job runs at a time.
- **A model that does not know an allele** (MixMHCpred and PRIME have their own allele lists; BigMHC is human only) casts no vote for it;
  the page names the skipped alleles.
- Runs over 3,000,000 predictions per class are refused (use fewer alleles or sizes). Stray characters in a FASTA make the windows containing
  them skip, and the page says so.

## The funnel: which peptides go on to the next model

Each arm has a **Funnel** setting:

- **Automatic** (default): the quick tools (MixMHCpred, NetMHCpan, MHCflurry, PRIME, IEDB, NetCTLpan, NetChop; MixMHC2pred, NetMHCIIpan) score every
  pair. BigMHC EL then scores only the pairs that can still become binders under the agreement rule, BigMHC IM only the pairs already called binders,
  CD4episcore only peptides that bind some allele. The gate is exact for the rule, so nothing a full run would have found is lost.
- **No funnel**: every model scores every pair; the filtering happens afterwards (slowest).
- **Manual**: you list stages. Stage *k* runs its models on the pairs still alive; its **pass rule** decides what goes on: everything, the pairs the
  stage's tools call a hit (any / majority / all of the tools, with your own cut-off or the strategy's), the best N peptides (or pairs), or the best P %.
  Models of one stage run side by side. A model appears in one stage only. Start from a template (cheap tools first, slow ones on the survivors) and edit.

Every pair in the tables carries the flags: `S1_In, S1_Pass, S2_In, ...` (scored in stage k / went on), `Funnel_Reached` (stages that scored it),
`Stopped_At` (the stage that stopped it, 0 = came through) and `Funnel_Alive`. A pair that was dropped is never a candidate, and the run details list
the real counts of each stage (pairs and peptides in and out, seconds). NetCTLpan and NetChop work on whole proteins, so putting them late saves nothing.

**Re-rank**: a finished run keeps its pair table. *Re-rank* (results page) applies another strategy or other rules to the same predictions in seconds
and saves the result as a new run. If the new rules are looser than the ones the slow tools were gated with, the page says some pairs have no slow-tool score.

**Clusters**: overlapping candidate windows (default 70 % of the shorter peptide, anchored on the best one) collapse onto one representative. Each
peptide table has `Cluster_ID, Cluster_Representative, Cluster_Overlap`; a representative also has `Cluster_Size` and `Cluster_Members` ("SEQ@start; ...").
Click a candidate to see the peptides clustered into it.

## Performance and no hanging

Settings > Performance (also `performance` in the json): `threads` per tool (0 = the tool decides), `parallel_models` (how many models of a stage run at
once, sharing the threads), `nice` (process priority, 5 by default so the browser stays smooth), `stall_seconds` (a tool with no CPU, disk or output
activity for this long is stopped and the run says which one; a merely slow tool keeps computing and is never stopped) and `timeout_minutes` (optional hard
limit per call). Cancel stops every running program. Standard input of every tool is closed, so none can wait for a keypress.

## Single model and Projects

- **Single model** runs one tool on its own (any of the twelve) and keeps everything it reports, with each peptide's position, sorted by the tool's
  main score. Nothing is combined or filtered, and the code path is the one the pipeline uses, so it is the way to double-check a number. The table
  can be cut by a threshold on any column and downloaded as CSV.
- **Projects** lists every run (pipeline or single) with its date, input, result, state, time and size. Each run is stored in the workspace with the
  exact settings and input, so after a restart you can **Open** it (all tabs, report and downloads work again), **Use settings** as a template for a new
  run, rename it, add a note, star it, download it or delete it. A run that was going when the server stopped is shown as *interrupted*.
- **Where the files are:** `workspace_dir` in Settings (default `<models folder>/engine_workspace`, beside the tools, so replacing the code folder
  does not touch it; `ENGINE_WORKSPACE` overrides it). Each run is a folder `runs/<id>/` with `run.json` (input, settings, dates, state, tools used,
  engine version) and its outputs.

## Consensus, agreement and population-weighted ranking

Every run ends with a consensus step (`consensus.py`, tested in `test_consensus.py`) that turns the tools' outputs into one ranked answer.

- **Votes.** Each tool casts one vote per peptide-allele pair. Tools with a percentile against random human peptides (NetMHCpan, MixMHCpred,
  MHCflurry, NetMHCIIpan, MixMHC2pred, PRIME, NetCTLpan) vote at an absolute cut-off (binding %Rank <= 2 for class I, 5 for class II by default).
  BigMHC-EL votes at score >= 0.2. Tools with no absolute scale for immunogenicity (BigMHC-IM, IEDB/Calis, CD4episcore) vote relatively: top q %
  among the plausible binders of the run (needs >= 8 of them, otherwise the tool gives no vote and says so).
- **Agreement** = votes / tools that scored the pair, in total and per group (binding, immunogenicity, processing). A pair is a **binder** when
  enough binding tools agree (strict majority by default); the binding consensus is the **median** %Rank of the percentile tools.
- **Population weighting.** Point `population.weights_csv` (Settings) at an HLA allele weightage CSV (`india_hla_allele_weightage.csv` from the
  Epitope Workbench, or drop it beside the engine folder; env `ENGINE_POPULATION_CSV` also works). A peptide's **coverage** is the share of the
  population carrying at least one allele it is a binder for (Hardy-Weinberg per region, loci independent, regions weighted as in the file).
  The weighted agreement and weighted binding rank use allele frequencies as weights. No file: coverage falls back to the share of your
  chosen alleles, and the page says so. Alleles with no data (HLA-DP, mouse) count in agreement but not in coverage.
- **Priority score** = weighted mean of the peptide's percentile standing in the run on coverage, agreement, binding and immunogenicity
  (0.35 / 0.25 / 0.20 / 0.20). Every component is a column, so any other ranking can be rebuilt in a spreadsheet.
- **Strategies** (Balanced, Strict, Max coverage, IEDB cut-offs, Wide net) state every rule; the "Adjust the rules" box overrides single rules.
  Overlapping windows collapse to one region (70 % overlap with the best one); a greedy **minimum panel** reaches the coverage goal.
- **Funnel & suggestions** show where peptides drop out, what each strategy would give on the same data, what single relaxations would add,
  and which common alleles are missing from your selection.
- **Downloads.** `Master CSV`: one row per peptide x allele, every tool's raw value, every vote, agreement counts, consensus, peptide-level
  coverage / rank / cluster / panel marks, sorted best first. Also peptides, candidates, panel and `run details` (json: the exact rules,
  tools, population and warnings behind the run).

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

## Installing MHCflurry

Open source (Apache 2.0), runs on torch, needs Python 3.10+, so it gets its own environment:

```bash
conda create -n mhcflurry python=3.11 -y && conda activate mhcflurry
pip install mhcflurry --extra-index-url https://download.pytorch.org/whl/cpu     # CPU-only torch is enough
mhcflurry-downloads fetch models_class1_presentation                              # 160 MB of trained models
```

The engine finds the environment (one named like `mhcflurry` is tried first), checks that the models are
downloaded, and shows the reason in Settings if not. For peptides cut from a protein it passes 15 residues of
flanking sequence on each side, which is what its processing predictor is designed for.

## Installing NetCTLpan and NetChop (IEDB's bundle)

Both are DTU programs that IEDB redistributes inside its next-generation T Cell Class I standalone
(free for academic and non-profit use; the bundle carries DTU's licence text). The engine drives the
bundle's own `tcell_mhci.py` in its own environment:

```bash
sudo apt install -y tcsh gawk
pip install uv                                    # or: curl -LsSf https://astral.sh/uv/install.sh | sh
mkdir -p ~/local_engines/iedb && cd ~/local_engines/iedb
wget https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/IEDB_NG_TC1-v0.2-beta.tar.gz    # 743 MB; the README in that folder names the current file
tar -xzf IEDB_NG_TC1-v0.2-beta.tar.gz
cd ng_tc1-v0.2-beta && uv sync && source .venv/bin/activate && ./configure
```

Then Settings -> **Test run** on `iedb_bundle`. Notes from testing the real tools:

- **They need the protein.** Cleavage depends on the residues around a peptide, so NetCTLpan and NetChop apply when a
  protein is cut into fragments, not to loose peptides (the form switches them off in that mode).
- **NetCTLpan scores 8-11-mers only**; other lengths get no NetCTLpan score and the page says so.
- **One bundle job per allele.** With several alleles in one job the bundle labels every row with the whole allele
  list, so the engine runs alleles separately (up to 4 in parallel).
- **Stray characters** (a `"` or `*` in a FASTA) make NetCTLpan return nothing for the whole protein, so they are
  replaced by `X` first, with a warning. The bundle's programs also hang if their standard input is left open;
  the engine always closes it.
- Speed (warm, 4 cores): the HEV FASTA (4,901 peptides, 3 alleles) takes about 2.5 minutes through all nine models;
  MHCflurry about 45 s, NetCTLpan about 45 s for 3 alleles (it grows with the number of alleles), BigMHC about 25 s
  per model on the chained pairs. The first run is slower while disk caches fill.

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

## MHC class II

Switch **MHC class II** at the top of the Pipeline form. Everything is the same idea as class I (one form, progress,
chaining, filters, per-peptide summary), with its own alleles, models and filters:

| Model | Source | What it gives |
|---|---|---|
| **NetMHCIIpan EL + BA** | DTU's program, run through IEDB's T Cell Class II standalone | EL score and %Rank; BA IC50 (nM) and %Rank. Versions 4.3 (default, newest), 4.2, or *IEDB default* (what IEDB's website returns: it reproduces the scores in IEDB web exports) |
| **MixMHC2pred 2.1** | Gfeller lab, GitHub | %Rank per allele (lower = better). Uses the 12-residue context around peptides cut from a protein |
| **CD4episcore** | IEDB (Dhanda et al. 2018), same standalone | CD4 T-cell immunogenicity percentile of the peptide alone, no HLA needed. **Lower = more immunogenic** |

- **Alleles:** DR (DRB1/3/4/5), DQ and DP as alpha/beta pairs (`HLA-DQA1*05:01/DQB1*02:01`), and mouse I-A / I-E. Presets: the
  IEDB class II reference panel (27 alleles: 15 DR, 6 DQ, 6 DP), DR only, the 7-allele DR set CD4episcore is built on, and mouse. The
  search box looks through MixMHC2pred's full allele table on the server.
- **Lengths:** 11-25 (default 15). Loose peptides may be 9-30 residues.
- **Chain:** CD4episcore only runs for peptides that pass the binding filters for at least one allele.
- A pair passes when every ticked model scored it and met its filter; an allele a tool does not know cannot pass, and the page
  names it. Class I and class II selections are remembered separately.

### Installing the class II tools

```bash
# MixMHC2pred (open source; free for academic use, for-profit use needs a licence, see its README)
cd ~/local_engines/gfeller
git clone --depth 1 https://github.com/GfellerLab/MixMHC2pred.git && cd MixMHC2pred
mkdir -p PWMdef && for sp in Human Mouse; do
  curl -LO https://mixmhc2pred.unil.ch/data/PWMdef/PWMdef_$sp.zip && unzip -oj PWMdef_$sp.zip -d PWMdef
done                                           # the git repository has the program but not the allele definitions

# IEDB T Cell Class II standalone: NetMHCIIpan and CD4episcore (1.6 GB; free for academic and non-profit use)
sudo apt install -y tcsh gawk perl r-base-core    # R and perl are needed by CD4episcore's NN-align
pip install uv
mkdir -p ~/local_engines/iedb && cd ~/local_engines/iedb
wget https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/IEDB_NG_TC2-v0.2.3-beta.tar.gz   # the folder's README names the current file
tar -xzf IEDB_NG_TC2-v0.2.3-beta.tar.gz && cd ng_tc2-v0.2.3-beta && uv sync
```

Then Settings -> **Test run** on `mixmhc2pred` and `iedb_bundle_ii`. Notes from testing the real tools:

- **`peptide_list` input fails** in the class II standalone (v0.2.3-beta), even in its own README example, and CD4episcore only
  accepts a protein. The engine sends every loose peptide as its own one-window sequence, grouped by length, which gives
  identical numbers (checked: 120 peptides reproduce IEDB's website CD4episcore scores exactly).
- **An unknown allele returns an empty table and no error** from NetMHCIIpan, and MixMHC2pred stops with a message and no output
  file; the engine checks alleles first and reports the ones it skipped.
- IEDB's "combined" immunogenicity method fails in this standalone and is not used.
- Speed (HA protein, 205 peptides, 8 alleles, three models): about 14 seconds.

