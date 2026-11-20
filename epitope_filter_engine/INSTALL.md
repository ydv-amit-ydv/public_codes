# Installing the Epitope Filter Engine on a new system

Complete, ordered instructions: the engine, then each of the eleven external predictors,
then how to verify every one of them.

The engine is a small Flask app. **It contains no predictors.** Each one is third-party
software with its own licence, downloaded separately, and the engine finds it on disk and
drives it as a subprocess with its own interpreter. So the engine itself installs in about
a minute, and the rest of this page is the predictors.

Tested on Ubuntu 22.04 / 24.04 and WSL2. Any Linux with `python3` ≥ 3.9 works
(developed and tested up to 3.13). macOS runs the engine but not DTU's Linux binaries
(NetMHCpan, NetMHCIIpan, NetCTLpan, NetChop).

---

## 0. The short version

Nothing but the engine, working in one minute, with one model (IEDB/Calis immunogenicity,
which is built in):

```bash
sudo apt update && sudo apt install -y python3-venv git
```

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git ~/public_codes
```

```bash
mkdir -p ~/local_engines && cd ~/public_codes/epitope_filter_engine && ENGINE_MODELS_DIR=~/local_engines bash start.sh
```

Open <http://127.0.0.1:5000>. Settings shows every model as not installed except
*IEDB immunogenicity*. Add predictors with the steps below, in any order; each one appears
in Settings as soon as its files are on disk — no reinstall, no restart.

Then install the predictors you want: **§3 is the one-command, open-licence set**
(MixMHCpred, PRIME, MixMHC2pred, BigMHC, MHCflurry, ImmuScope, TLimmuno2) and **§4 is the
set that needs a licence form or a large download** (NetMHCpan, the two IEDB bundles).

---

## 1. Where things go

One folder holds the tools; the engine is discovered relative to it. Nothing is installed
system-wide, and nothing is written inside the repository clone.

```
~/local_engines/                 <- ENGINE_MODELS_DIR: searched for the tools
    netmhc/                      NetMHCpan 4.2            (§4.1, DTU licence form)
    bigmhc/                      BigMHC clone             (§3.3)
    gfeller/
        MixMHCpred/              MixMHCpred 3.0           (§3.1)
        PRIME/                   PRIME 2.1                (§3.2)
        MixMHC2pred/             MixMHC2pred 2.1          (§3.4)
    iedb/
        ng_tc1-v0.2-beta/        NetCTLpan + NetChop      (§4.2, 743 MB)
        ng_tc2-v0.2.3-beta/      NetMHCIIpan + CD4episcore (§4.3, 1.6 GB)
    ImmuScope/                   ImmuScope-IM             (§3.6)
    TLimmuno2/                   TLimmuno2                (§3.7)

~/public_codes/epitope_filter_engine/    <- the engine (this folder); `git pull` to update
    .venv/                       made by start.sh on first run (flask + pandas)
    engine_workspace/            runs, caches, learned timings  (workspace_dir)
```

**The folder names above do not matter.** Discovery searches for each tool by the files
inside it (`netMHCpan` beside `Linux_x86_64/`, `src/predict.py` beside `models/`, and so on),
up to 6 levels deep, shallowest match wins. Put a tool anywhere under `~/local_engines` and
it is found. An explicit path in `engine_config.json` always overrides discovery.

Two ways to point the engine at that folder — pick one:

| | How | When |
|---|---|---|
| **A (recommended)** | Run in place from the repo clone with `ENGINE_MODELS_DIR=~/local_engines` | `git -C ~/public_codes pull` updates the engine. Put the variable in `~/.bashrc` and plain `bash start.sh` works. |
| **B** | Copy the engine *into* `~/local_engines/runner_codes/` and run it there with no variable | Discovery then defaults to the parent folder. Matches the older `local_engine_codes` instructions; update means copying again. |

For **B**: `mkdir -p ~/local_engines/runner_codes && cp -r ~/public_codes/epitope_filter_engine/. ~/local_engines/runner_codes/`

---

## 2. System packages

Install these once. The right-hand column says which predictors need each one, so you can
skip the rows for tools you are not installing.

```bash
sudo apt update && sudo apt install -y \
    python3-venv git curl wget unzip \
    g++ perl tcsh gawk r-base-core
```

| Package | Needed by |
|---|---|
| `python3-venv`, `git`, `curl`/`wget`, `unzip` | the engine, and fetching everything else |
| `g++` | PRIME (the engine rebuilds its C++ binary for Linux on first use) |
| `perl` | PRIME's script, and CD4episcore's NN-align |
| `tcsh`, `gawk` | the IEDB bundles (NetCTLpan, NetChop, NetMHCIIpan, CD4episcore) |
| `r-base-core` | CD4episcore's NN-align |
| `qemu-user` | only if CD4episcore segfaults on a very new Intel CPU — see §7 |

`uv` is needed by both IEDB bundles (§4.2, §4.3):

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Conda is **optional**. It is one convenient way to give BigMHC and MHCflurry their own
Python; `python3 -m venv` works just as well, and §3 shows both. The engine never imports
torch or TensorFlow itself, so these environments never collide with each other or with the
server.

> **If you use conda, read §7's note on shell tools.** With a conda environment active,
> tcsh/perl/R can load conda's libraries instead of the system's and fail oddly. The engine
> already starts those programs with a cleaned environment, which handles it.

---

## 3. The open-licence predictors

Each block is independent and safe to re-run. All of them are free for academic use; the
Gfeller lab tools require a licence from the Ludwig Institute for for-profit use (see each
repository's README). Start here:

```bash
mkdir -p ~/local_engines/gfeller && cd ~/local_engines/gfeller
```

### 3.1 MixMHCpred 3.0 — class I presentation

```bash
git clone --depth 1 https://github.com/GfellerLab/MixMHCpred.git
```

Needs only numpy and pandas, which the engine's own environment already has. Nothing to
build.

### 3.2 PRIME 2.1 — class I immunogenicity

```bash
git clone --depth 1 https://github.com/GfellerLab/PRIME.git
```

PRIME **calls MixMHCpred**, so §3.1 is a prerequisite. The repository ships a binary for
other systems; the engine rebuilds `lib/PRIME.x` for Linux the first time PRIME runs
(hence `g++`).

### 3.3 BigMHC — class I eluted-ligand and immunogenicity

```bash
cd ~/local_engines
git clone --depth 1 https://github.com/KarchinLab/bigmhc.git
```

BigMHC needs torch. Give it its own environment:

```bash
cd ~/local_engines/bigmhc
python3 -m venv .venv
.venv/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu   # CPU build is enough
.venv/bin/pip install numpy pandas psutil
```

The engine looks for an interpreter that can import `torch, numpy, pandas, psutil`: a
`.venv` inside or beside the clone, then a conda environment (one named like `bigmhc`
first). Override with `BIGMHC_PYTHON=/path/to/bin/python`.

Equivalent with conda:
`conda create -n bigmhc python=3.11 -y && conda run -n bigmhc pip install torch --index-url https://download.pytorch.org/whl/cpu numpy pandas psutil`

### 3.4 MixMHC2pred 2.1 — class II presentation

The git repository ships the program but **not** the allele definitions, which are a
separate download:

```bash
cd ~/local_engines/gfeller
git clone --depth 1 https://github.com/GfellerLab/MixMHC2pred.git && cd MixMHC2pred
mkdir -p PWMdef
for sp in Human Mouse; do
  curl -LO https://mixmhc2pred.unil.ch/data/PWMdef/PWMdef_$sp.zip
  unzip -oj PWMdef_$sp.zip -d PWMdef
done
```

Without `PWMdef/` the engine reports MixMHC2pred as not ready and names this URL, and a run
that selects it is refused with that reason rather than failing silently.

### 3.5 MHCflurry 2 — class I affinity, processing, presentation

Open source (Apache 2.0), needs Python ≥ 3.10 and torch, so give it its own environment:

```bash
conda create -n mhcflurry python=3.11 -y && conda activate mhcflurry
pip install mhcflurry --extra-index-url https://download.pytorch.org/whl/cpu
mhcflurry-downloads fetch models_class1_presentation          # 160 MB of trained models
conda deactivate
```

Without conda:

```bash
python3 -m venv ~/local_engines/mhcflurry-env
~/local_engines/mhcflurry-env/bin/pip install mhcflurry --extra-index-url https://download.pytorch.org/whl/cpu
~/local_engines/mhcflurry-env/bin/mhcflurry-downloads fetch models_class1_presentation
export MHCFLURRY_PYTHON=~/local_engines/mhcflurry-env/bin/python
```

The engine tries an environment named like `mhcflurry` first and checks that the models are
downloaded, saying so in Settings if they are not. `MHCFLURRY_PYTHON` overrides.

### 3.6 ImmuScope-IM — class II immunogenicity per peptide-allele pair

Code (GPL-3.0) plus weights (CC-BY-4.0) from the authors' Zenodo record 14810445:

```bash
cd ~/local_engines
git clone --depth 1 https://github.com/shenlongchen/ImmuScope.git && cd ImmuScope

# weights (149.8 MB) -> weights/IM, weights/EL, weights/CD4
curl -L -o w.tgz https://zenodo.org/api/records/14810445/files/ImmuScope-weights.tar.gz/content
tar -xzf w.tgz && rm w.tgz

# only one file is needed out of the 3 GB data archive: stream it, do not keep the archive
mkdir -p data/raw
curl -sL https://zenodo.org/api/records/14810445/files/ImmuScope-data.tar.gz/content \
  | tar -xz --wildcards '*pseudosequence.2023.dat'
find . -name pseudosequence.2023.dat      # must end up at data/raw/ — move it if it landed in data/data/raw/

python3 -m venv .venv
.venv/bin/pip install torch --index-url https://download.pytorch.org/whl/cpu
.venv/bin/pip install numpy pandas h5py scikit-learn ruamel.yaml tqdm click
```

Override the interpreter with `IMMUSCOPE_PYTHON`.

### 3.7 TLimmuno2 — class II immunogenicity, second network

The models are in the repository. It needs **TensorFlow 2.15 (Keras 2)** — Keras 3 cannot
load these models:

```bash
cd ~/local_engines
git clone --depth 1 https://github.com/XSLiuLab/TLimmuno2.git && cd TLimmuno2
python3 -m venv .venv
.venv/bin/pip install "tensorflow-cpu==2.15.1" "numpy<2" pandas
```

Override with `TLIMMUNO2_PYTHON`.

> TensorFlow 2.15 wants Python 3.9–3.11. On a newer system, make this one environment with
> an older interpreter (`python3.11 -m venv .venv`, or a conda env named `tlimmuno2`) — the
> engine only ever runs TLimmuno2 through it, so the rest of the install is unaffected.

### 3.8 IEDB class I immunogenicity (Calis 2013) — already installed

Nothing to download. The model is a 20-value amino-acid scale times nine position weights,
implemented natively in `MHC1/Immunogenicity_scripts/calis_immunogenicity.py`, and verified
against IEDB's own `immunogenicity_predictor.py` on 28,000 random peptides and against
11,781 scores from an IEDB website export (maximum difference 0.0).

---

## 4. The predictors that need a licence form or a large download

### 4.1 NetMHCpan 4.2 — class I binding affinity and eluted ligand

DTU academic-licence software: you must accept the licence on their site, so this step
cannot be scripted.

1. Go to <https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/> and request the
   **Linux** version. You receive the software tarball and a link to the data tarball.
2. Extract both into `~/local_engines/netmhc/` and follow the `netMHCpan-4.2.readme`
   inside: it tells you to unpack the data file into the same folder.

```bash
mkdir -p ~/local_engines/netmhc && cd ~/local_engines/netmhc
tar -xzf ~/Downloads/netMHCpan-4.2b.Linux.tar.gz
cd netMHCpan-4.2 && tar -xzf ~/Downloads/data.tar.gz          # file names depend on your download
```

**You do not need to edit the script.** DTU's `netMHCpan` ships with a placeholder `NMHOME`
path; rather than modify your copy, the engine runs a patched copy from its scratch
directory whenever the placeholder is wrong. Editing it yourself is fine too.

The engine finds any folder holding a file named `netMHCpan` beside a `Linux_x86_64/`
directory, nested at any depth, so the extra `netMHCpan-4.2/` level above is harmless.

### 4.2 NetCTLpan 1.1 and NetChop 3.1 — class I antigen processing

DTU programs that IEDB redistributes inside its next-generation T Cell Class I standalone
(free for academic and non-profit use; the bundle carries DTU's licence text). The engine
drives the bundle's own `tcell_mhci.py` in the bundle's own environment.

```bash
mkdir -p ~/local_engines/iedb && cd ~/local_engines/iedb
wget https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/IEDB_NG_TC1-v0.2-beta.tar.gz   # 743 MB
tar -xzf IEDB_NG_TC1-v0.2-beta.tar.gz
cd ng_tc1-v0.2-beta && uv sync && ./configure
```

Check the `README` in that downloads folder for the current file name if the one above is
gone. Needs `tcsh` and `gawk` (§2).

Two things worth knowing before you use them:

- **They need the whole protein.** Cleavage depends on the residues around a peptide, so
  these two apply when a protein is cut into fragments, not to a loose peptide list; the
  form switches them off in that mode.
- **NetCTLpan scores 8–11-mers only.** Other lengths get no score and the page says so.

### 4.3 NetMHCIIpan and CD4episcore — class II binding and immunogenicity

IEDB's T Cell Class II standalone, same licence terms, 1.6 GB.

```bash
cd ~/local_engines/iedb
wget https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/IEDB_NG_TC2-v0.2.3-beta.tar.gz
tar -xzf IEDB_NG_TC2-v0.2.3-beta.tar.gz
cd ng_tc2-v0.2.3-beta && uv sync
```

Needs `tcsh`, `gawk`, `perl` and `r-base-core` (§2). This bundle provides NetMHCIIpan
EL + BA in versions 4.3 (default), 4.2 and *IEDB default*, plus CD4episcore.

**Apply the fix in §7.1 now** — in this release CD4episcore does not work without it.

---

## 5. Verify every model

```bash
cd ~/public_codes/epitope_filter_engine && ENGINE_MODELS_DIR=~/local_engines bash start.sh
```

In the browser, **Settings** shows *In use right now* (every path that was discovered) and a
readiness checklist per model, with the reason when something is missing. Then press
**Test run** for each of these keys — it predicts a known peptide end to end through the real
tool:

| Key | Covers | Installed in |
|---|---|---|
| `netmhcpan` | NetMHCpan 4.2 | §4.1 |
| `bigmhc` | BigMHC EL + IM | §3.3 |
| `mixmhcpred` | MixMHCpred 3.0 | §3.1 |
| `prime` | PRIME 2.1 | §3.2 |
| `mhcflurry` | MHCflurry 2 | §3.5 |
| `iedb_bundle` | NetCTLpan, NetChop | §4.2 |
| `mixmhc2pred` | MixMHC2pred 2.1 | §3.4 |
| `iedb_bundle_ii` | NetMHCIIpan, CD4episcore | §4.3 |
| `immuscope` | ImmuScope-IM | §3.6 |
| `tlimmuno2` | TLimmuno2 | §3.7 |

The first BigMHC call loads its models from disk and can take minutes on a cold cache;
later calls are quick.

Two numbers to check by hand if you suspect a tool is wrong:

- **GILGFVFTL on HLA-A\*02:01** — what *Test run* uses for class I.
- **PKYVKQNTLKLATGM on DRB1\*01:01 must give CD4episcore 96.32.** If it does not, §7.1.

And the test suite, which needs no predictors at all (it stubs them):

```bash
.venv/bin/python test_engine.py        # config, server, API, funnel, workspace, both wrappers
.venv/bin/python test_consensus.py     # votes, agreement, coverage, priority score
.venv/bin/python test_host.py          # host-proteome search against brute force
```

All three must print `... passed`. `test_engine.py` skips its NetCTLpan and class II
sections when tcsh/gawk/perl/R are absent and says so. To exercise real BigMHC instead of
the stub: `BIGMHC_PYTHON=<env python> BIGMHC_HOME=<clone> .venv/bin/python test_engine.py`.

---

## 6. Two optional data files

Both are optional; the engine says so in the page when either is missing.

**Population HLA frequencies** — needed for population-weighted coverage and for the
priority score's coverage term. Point `population.weights_csv` (Settings) at an HLA allele
weightage CSV, for example `india_hla_allele_weightage.csv` from the Epitope Workbench, or
drop it beside the engine folder. `ENGINE_POPULATION_CSV` also works. With no file,
coverage falls back to the share of the alleles *you* selected, and the page says so.

**Host proteome** — for the self-similarity check that flags candidates identical, or one
residue away, from a protein the host already makes. Settings → *Host proteome* →
**Download** (about 8 MB, kept in the workspace), or set `host_proteome.fasta` to your own
FASTA (plain or `.gz`). Set `host_proteome.organism` to `mouse` for mouse studies. A run
whose rules must *remove* self-peptides will not start until this is downloaded; *flag only*
(the default) just notes that it was not checked.

---

## 7. Troubleshooting the real faults

These are the failures actually hit while testing the tools, with the fix.

### 7.1 CD4episcore returns nothing: "loop of ufunc does not support argument 0 of type float"

Two independent causes in IEDB's class II bundle v0.2.3-beta, both seen on a
Core Ultra 9 / WSL machine. Run these from inside `ng_tc2-v0.2.3-beta/`.

**(a)** `src/imscore.py` builds the R path from bytes, giving `-R b/usr/bin/R`:

```bash
sed -i 's|R_path = subprocess.check_output("which R", shell=True).rstrip()$|R_path = subprocess.check_output("which R", shell=True).rstrip().decode()|' src/imscore.py
grep -n "R_path =" src/imscore.py        # must show exactly one .decode()
```

Run that `sed` **once only**.

**(b)** NN-align's statically linked binaries (`src/nnalign/bin/Linux_x86_64/*`, built for
kernel 2.6) segfault on very new Intel CPUs — "Segmentation fault" before
"Uploaded 0 peptides". Run them under qemu:

```bash
sudo apt install -y qemu-user
cd src/nnalign/bin/Linux_x86_64
for f in *; do
  [ -f "$f" ] && file "$f" | grep -q ELF || continue
  mv "$f" "$f.real"
  printf '#!/bin/sh\nexec qemu-x86_64 "$(dirname "$0")/%s.real" "$@"\n' "$f" > "$f"
  chmod +x "$f"
done
```

Confirm with the hand check: `PKYVKQNTLKLATGM` + `DRB1*01:01` → **96.32**.

### 7.2 PRIME hangs

PRIME 2.1's C++ program never returns for **HLA-C\*07:06** (it spins at 100 % CPU). The
engine skips that allele for PRIME with a note. PRIME is also run one allele per call under
a time limit, so any other allele that hangs is given up on instead of stalling the run.

### 7.3 tcsh prints "libtinfo.so.6: no version information available", or perl/R fail silently

A conda environment is active and these programs are loading conda's libraries instead of
the system's. The engine already starts tcsh, gawk, perl and R with a cleaned environment
(no `LD_LIBRARY_PATH`, no conda variables, no perl/R paths), which is what fixed
CD4episcore failing with the ufunc error above. If you still see it, `perl -v` and
`R --version` must work in a plain terminal outside conda.

### 7.4 An allele gets no score

A model that does not know an allele casts no vote for it and the page names the skipped
alleles: MixMHCpred and PRIME have their own allele lists, BigMHC is human only, and
HLA-DP and mouse alleles have no population frequency data (they count in agreement but
not in coverage). An unknown allele makes NetMHCIIpan return an empty table with no error
and MixMHC2pred stop with a message and no output file, so the engine checks alleles first
and reports which it skipped.

### 7.5 A FASTA with stray characters

A `"` or `*` in a sequence makes NetCTLpan return nothing for the whole protein, so such
characters are replaced by `X` first, with a warning. Windows containing other stray
characters are skipped and the page says so.

### 7.6 A model fails mid-run

The run logs it (`✗ CD4episcore failed and is skipped: ...`), closes that stage with the
reason, warns on the results page and continues with the other models. A class that
finished with a failed model is not saved for reuse, so pressing Run again retries it.
Only if *every* selected model fails does the run fail.

### 7.7 `python3 -m venv failed`

`sudo apt install python3-venv`. `start.sh` needs nothing else.

### 7.8 Port 5000 is taken

`PORT=5001 bash start.sh`. The server binds 127.0.0.1 only.

---

## 8. Environment variables

All optional; every one of them can be set in `engine_config.json` instead.

| Variable | Effect |
|---|---|
| `ENGINE_MODELS_DIR` | Folder searched for the tools (default: the engine's parent folder) |
| `ENGINE_CONFIG` | Path of the settings JSON (default: `engine_config.json` beside the engine) |
| `ENGINE_WORKSPACE` | Where runs, caches and timings are kept (default: `<models folder>/engine_workspace`) |
| `ENGINE_POPULATION_CSV` | HLA frequency CSV (§6) |
| `PORT` | Server port (default 5000) |
| `BIGMHC_PYTHON`, `BIGMHC_HOME` | BigMHC interpreter and clone |
| `MHCFLURRY_PYTHON` | MHCflurry interpreter |
| `IMMUSCOPE_PYTHON` | ImmuScope interpreter |
| `TLIMMUNO2_PYTHON` | TLimmuno2 interpreter |
| `ENGINE_DEBUG` | Verbose server logging |

`engine_config.json` is only needed to override discovery or change options such as
`device`, thread counts or the performance profile. Copy `engine_config.example.json` —
every path in it is `"auto"`. Settings edits it live: **Save** validates it (bad JSON or a
wrong type is rejected and the file left untouched) and applies it on the next request,
with no restart.

---

## 9. Updating and removing

```bash
git -C ~/public_codes pull                      # layout A: the engine is updated, nothing else to do
```

With layout B, copy again: `cp -r ~/public_codes/epitope_filter_engine/. ~/local_engines/runner_codes/`

Your runs and settings live in the workspace (`engine_workspace/`), outside the code
folder, so replacing the code never touches them. To remove the engine, delete the clone
and `~/local_engines`. Nothing was installed system-wide except the apt packages in §2.

---

## 10. Headless use

No browser needed:

```bash
# save a run's form as JSON from the page (Save settings), then
python3 run_engine.py epitope_settings.json -o results/
python3 run_engine.py epitope_settings.json --fasta other.fa -o results/   # same settings, another protein
```

The page's API takes the same JSON: `POST /api/pipeline/run`.

---

## Licences

The engine's own code is in this repository. **Every predictor is licensed separately and
none is redistributed here**: NetMHCpan, NetMHCIIpan, NetCTLpan and NetChop are DTU
academic-licence software (the last two via IEDB's bundles, which carry DTU's licence text);
MixMHCpred, PRIME and MixMHC2pred are free for academic use and need a Ludwig Institute
licence for for-profit use; BigMHC, MHCflurry (Apache 2.0), ImmuScope (GPL-3.0, weights
CC-BY-4.0) and TLimmuno2 are open source. The IEDB bundles are free for academic and
non-profit use. The Calis constants come from IEDB's package (Non-Profit Open Software
License 3.0, La Jolla Institute for Immunology); the code here is an independent
implementation of the published model. Check each tool's own terms before any commercial
use.
