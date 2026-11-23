# Installing the Epitope Filter Engine

## Clone it, run one command, done

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git
```

```bash
cd public_codes/epitope_filter_engine && bash start.sh
```

That is the whole installation. The first run:

- installs the system packages the shell-based predictors need (it asks once, with sudo),
- makes the engine's own environment (conda if you have it, otherwise a virtualenv),
- installs `uv`, then downloads and sets up **every predictor that can be installed without
  a licence form** — MixMHCpred, PRIME, MixMHC2pred, BigMHC, MHCflurry, ImmuScope, TLimmuno2
  and both IEDB bundles (NetMHCIIpan, CD4episcore, NetCTLpan, NetChop),
- applies the CD4episcore fix that release needs,
- and opens the engine on <http://localhost:5000>.

Every run after that just opens the engine, in a second. Re-running is always safe: anything
already installed is skipped.

**One thing you must do by hand: NetMHCpan 4.2.** DTU require you to accept their licence on
their own site, so no script may download it. Get the Linux version from
<https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/>, extract it anywhere under
`~/local_engines`, and the engine finds it. Everything else works without it, and the first
run tells you this again at the end. See §6.1.

```
Options you will probably never need:
  bash start.sh --light        skip the two big IEDB bundles (~2.4 GB)
  bash start.sh --no-setup     just start, skip the first-run install
  bash start.sh --setup-only   install, do not open the server
  bash start.sh --redo-setup   run the install again
  bash start.sh --dir DIR      put the predictors somewhere else
  bash start.sh --dry-run      print what the install would do, change nothing
```

Allow about **20 GB** of disk, or about 4 GB with `--light`. The predictors go in
`~/local_engines` (or beside the clone, if you already had them there).

**Written for Ubuntu, including WSL2**, which is what it is developed and tested on. On WSL,
read §2.1 first — where you put the files matters a great deal for speed. macOS runs the
engine but not DTU's Linux binaries (NetMHCpan, NetMHCIIpan, NetCTLpan, NetChop).

---

## The rest of this page

You do not need it to install. It is here for when something breaks, or when you want to do a
step yourself.

The engine is a small Flask app and **contains no predictors**. Each one is third-party
software with its own licence, installed separately, which the engine finds on disk and runs
as a subprocess **in its own environment** — so torch, TensorFlow and the two IEDB bundles
never have to agree about anything.

| | |
|---|---|
| §2 | WSL2, and the system packages the installer puts in |
| §3 | conda (optional, but the smoother path) |
| §4 | the engine itself, and where everything lives |
| §5 | each open-licence predictor by hand |
| §6 | **NetMHCpan** and the IEDB bundles |
| §7 | **verify every model** |
| §8 | troubleshooting the faults we actually hit |
| §9–12 | optional data, variables, updating, licences |

---

## 2. WSL2 and system packages

### 2.1 WSL2

In **Windows PowerShell as administrator**:

```powershell
wsl --install -d Ubuntu-22.04
```

Then open *Ubuntu* from the Start menu and create your user. Check you are on WSL **2**
(`wsl -l -v` in PowerShell must show `VERSION 2`; `wsl --set-version Ubuntu-22.04 2` fixes
it).

Three WSL things that matter here:

**Keep everything in the Linux filesystem.** Install the engine and the tools under your
Linux home (`~/public_codes`, `~/local_engines`) and **never under `/mnt/c/...`**. Crossing
to the Windows drive is 10–50× slower for the many small reads these tools do, and it does
not preserve the execute permission that DTU's and IEDB's shell programs need. To reach
your results from Windows, go the other way: type `\\wsl$\Ubuntu-22.04\home\<you>` in
Explorer, or run `explorer.exe .` from the folder.

**Give WSL enough RAM.** BigMHC, MHCflurry, ImmuScope and TLimmuno2 each load a neural
network; the default WSL allowance can be too small on a 16 GB machine. Create
`C:\Users\<you>\.wslconfig`:

```ini
[wsl2]
memory=12GB
processors=8
swap=4GB
```

Then `wsl --shutdown` in PowerShell and reopen Ubuntu. The engine's *Performance* profile
(Settings) chooses **Maximum** automatically with 16+ cores and 16+ GB free, so what you
grant here decides how parallel the runs are.

**localhost just works.** WSL2 forwards it, so the server started in Ubuntu is reachable at
<http://localhost:5000> from Windows with no extra setup. It binds `127.0.0.1` only, so it
is not exposed to your network.

### 2.2 System packages

```bash
sudo apt update && sudo apt install -y \
    git curl wget unzip \
    g++ perl tcsh gawk r-base-core
```

| Package | Needed by |
|---|---|
| `git`, `curl`/`wget`, `unzip` | fetching everything |
| `g++` | PRIME (the engine rebuilds its C++ binary for Linux on first use) |
| `perl` | PRIME's script, and CD4episcore's NN-align |
| `tcsh`, `gawk` | the IEDB bundles (NetCTLpan, NetChop, NetMHCIIpan, CD4episcore) |
| `r-base-core` | CD4episcore's NN-align |
| `qemu-user` | only if CD4episcore segfaults on a very new Intel CPU — §8.1(b) |

`uv` is needed by both IEDB bundles (§6.2, §6.3) and is not an apt package:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Reopen the shell afterwards so `uv` is on `PATH`.

---

## 3. conda (Miniforge)

**Miniforge** is the recommended installer: it is conda with `conda-forge` as the default
channel, and it avoids the licensing question around Anaconda's own `defaults` channel for
commercial and large-organisation use. Any conda (Miniconda, Anaconda, Mambaforge) works.

```bash
curl -L -o /tmp/miniforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
bash /tmp/miniforge.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash
exec bash -l                       # reload the shell so `conda` works
conda --version
```

If you would rather not have conda activate its base environment in every shell:

```bash
conda config --set auto_activate_base false
```

**The engine does not need conda to be active.** It locates each tool's environment by
name, so you never have to activate anything before a run — see §7.

> **No conda, or you prefer virtualenvs?** Everything below works with
> `bash install_tools.sh --venv`, which makes a `.venv` inside each tool's folder instead.
> The one place conda is genuinely easier is TLimmuno2, which needs Python ≤ 3.11
> (§5.7). Skip to §4.

---

## 4. The engine itself

```bash
git clone https://github.com/ydv-amit-ydv/public_codes.git ~/public_codes
cd ~/public_codes/epitope_filter_engine
conda env create -f environment.yml          # creates "epitope-engine": python 3.11, flask, pandas
```

```bash
ENGINE_MODELS_DIR=~/local_engines bash start.sh
```

`start.sh` prints the interpreter it chose and picks it in this order:

1. `$ENGINE_PYTHON`, if you set it
2. the conda environment that is **currently active**, if it has flask and pandas
3. the conda environment named **`epitope-engine`**, active or not ← the normal case
4. `./.venv`, which it creates with flask and pandas if nothing else is found

So `conda activate epitope-engine` is optional. Put the models directory in your profile
and even that argument goes away:

```bash
echo 'export ENGINE_MODELS_DIR=~/local_engines' >> ~/.bashrc
```

Open <http://localhost:5000>. It already works: **IEDB class I immunogenicity (Calis 2013)**
is implemented natively in the engine, so it needs no download and is the one model that is
ready before you install anything. Settings lists every other model as not installed, with
the reason.

### Where things live

```
~/local_engines/                 <- ENGINE_MODELS_DIR: searched for the tools
    netmhc/                      NetMHCpan 4.2             §6.1  (DTU licence form)
    bigmhc/                      BigMHC clone              §5.4
    gfeller/
        MixMHCpred/              MixMHCpred 3.0            §5.1
        PRIME/                   PRIME 2.1                 §5.2
        MixMHC2pred/             MixMHC2pred 2.1           §5.3
    iedb/
        tc1/                     NetCTLpan + NetChop       §6.2  (743 MB)
        tc2/                     NetMHCIIpan + CD4episcore §6.3  (1.6 GB)
    ImmuScope/                   ImmuScope-IM              §5.6
    TLimmuno2/                   TLimmuno2                 §5.7
    engine_workspace/            runs, caches, learned timings, host proteome

~/public_codes/epitope_filter_engine/     <- the engine; `git pull` to update
```

**The folder names do not matter.** Discovery searches for each tool by the files inside it
(a `netMHCpan` script beside `Linux_x86_64/`, `src/predict.py` beside `models/`, and so on),
up to 6 levels deep, shallowest match wins. Conda environments are found by name. An
explicit path in `engine_config.json` always overrides discovery.

Your runs and settings live in `engine_workspace/`, **outside** the code folder, so
updating or replacing the engine never touches them.

---

## 5. The open-licence predictors — one command

```bash
cd ~/public_codes/epitope_filter_engine
bash install_tools.sh
```

That clones and sets up seven predictors, each in its own conda environment named after it,
into `~/local_engines`:

| Environment | Tool | Class |
|---|---|---|
| *(none needed)* | MixMHCpred 3.0 | I |
| *(none needed)* | PRIME 2.1 | I |
| *(none needed)* | MixMHC2pred 2.1 | II |
| `bigmhc` | BigMHC EL + IM (torch) | I |
| `mhcflurry` | MHCflurry 2 (torch) | I |
| `immuscope` | ImmuScope-IM (torch) | II |
| `tlimmuno2` | TLimmuno2 (TensorFlow 2.15) | II |

Useful flags: `--dry-run` (print, change nothing), `--only bigmhc,tlimmuno2`,
`--dir /other/path`, `--venv` (virtualenvs instead of conda), `--iedb` (also §6.2 and §6.3).
Re-running skips whatever is already there, so it is safe after a failure.

All seven are free for academic use. The three Gfeller lab tools require a licence from the
Ludwig Institute for for-profit use — see each repository's README.

The rest of this section is what the script does, if you would rather do it by hand or need
to repair one tool.

### 5.1 MixMHCpred 3.0 — class I presentation

```bash
mkdir -p ~/local_engines/gfeller && cd ~/local_engines/gfeller
git clone --depth 1 https://github.com/GfellerLab/MixMHCpred.git
```

Needs only numpy and pandas, which the engine's own environment already has. Nothing to
build, no environment of its own.

### 5.2 PRIME 2.1 — class I immunogenicity

```bash
cd ~/local_engines/gfeller
git clone --depth 1 https://github.com/GfellerLab/PRIME.git
```

PRIME **calls MixMHCpred**, so §5.1 is a prerequisite. The repository ships a binary for
other systems; the engine rebuilds `lib/PRIME.x` for Linux the first time PRIME runs, which
is why `g++` is in §2.2.

### 5.3 MixMHC2pred 2.1 — class II presentation

The repository ships the program but **not** the allele definitions, which are a separate
download:

```bash
cd ~/local_engines/gfeller
git clone --depth 1 https://github.com/GfellerLab/MixMHC2pred.git && cd MixMHC2pred
mkdir -p PWMdef
for sp in Human Mouse; do
  curl -LO https://mixmhc2pred.unil.ch/data/PWMdef/PWMdef_$sp.zip
  unzip -oj PWMdef_$sp.zip -d PWMdef && rm PWMdef_$sp.zip   # -j matters: the zip holds a
done                                                       # PWMdef_Human/ folder, and the
                                                           # engine wants the .txt files flat
```

Without `PWMdef/` the engine reports MixMHC2pred as not ready and names this URL, and a run
that selects it is refused with that reason instead of failing silently.

### 5.4 BigMHC — class I eluted ligand and immunogenicity

```bash
cd ~/local_engines
git clone --depth 1 https://github.com/KarchinLab/bigmhc.git
conda create -y -n bigmhc python=3.11
conda run -n bigmhc pip install torch --index-url https://download.pytorch.org/whl/cpu numpy pandas psutil
```

The CPU build of torch is enough. The engine looks for an interpreter that can import
`torch, numpy, pandas, psutil`, preferring a conda environment whose name contains
`bigmhc`; `BIGMHC_PYTHON` overrides.

BigMHC is the slowest model (about 26 s to start plus 9 ms per peptide-allele pair, for each
of EL and IM), which is why the engine chains it: EL scores only pairs that can still become
binders and IM only pairs already called binders.

### 5.5 MHCflurry 2 — class I affinity, processing, presentation

Apache 2.0, needs Python ≥ 3.10 and torch:

```bash
conda create -y -n mhcflurry python=3.11
conda run -n mhcflurry pip install mhcflurry --extra-index-url https://download.pytorch.org/whl/cpu
conda run -n mhcflurry mhcflurry-downloads fetch models_class1_presentation     # 160 MB
```

MHCflurry has no repository to clone, so the environment is the whole install — another
reason the conda path is tidier here. With `--venv` you must set
`MHCFLURRY_PYTHON=~/local_engines/mhcflurry-env/bin/python`, which the script tells you.

The engine checks that the models are downloaded and says so in Settings if they are not.
For peptides cut from a protein it passes 15 residues of flanking sequence on each side,
which is what the processing predictor expects.

### 5.6 ImmuScope-IM — class II immunogenicity per peptide-allele pair

Code (GPL-3.0) plus weights (CC-BY-4.0) from the authors' Zenodo record 14810445:

```bash
cd ~/local_engines
git clone --depth 1 https://github.com/shenlongchen/ImmuScope.git && cd ImmuScope

curl -L -o w.tgz https://zenodo.org/api/records/14810445/files/ImmuScope-weights.tar.gz/content
tar -xzf w.tgz && rm w.tgz                 # 149.8 MB -> weights/IM, weights/EL, weights/CD4

# only one file is needed out of a 3 GB archive: stream it, never keep the archive
mkdir -p data/raw
curl -sL https://zenodo.org/api/records/14810445/files/ImmuScope-data.tar.gz/content \
  | tar -xz --wildcards '*pseudosequence.2023.dat'
find . -name pseudosequence.2023.dat       # must end at data/raw/ — move it if it landed in data/data/raw/

conda create -y -n immuscope python=3.11
conda run -n immuscope pip install torch --index-url https://download.pytorch.org/whl/cpu \
    numpy pandas h5py scikit-learn ruamel.yaml tqdm click
```

`IMMUSCOPE_PYTHON` overrides. The engine averages the 10 IM networks the way the authors'
test script does.

### 5.7 TLimmuno2 — class II immunogenicity, second network

The models are in the repository. It needs **TensorFlow 2.15 (Keras 2)** — Keras 3 cannot
load these models — and TF 2.15 has no wheel above Python 3.11, which is exactly the case
conda handles best:

```bash
cd ~/local_engines
git clone --depth 1 https://github.com/XSLiuLab/TLimmuno2.git
conda create -y -n tlimmuno2 python=3.11
conda run -n tlimmuno2 pip install "tensorflow-cpu==2.15.1" "numpy<2" pandas
```

`TLIMMUNO2_PYTHON` overrides. With `--venv` on a system whose `python3` is 3.12 or newer,
the script looks for `python3.11`/`3.10`/`3.9` and tells you to install one if there is
none.

### 5.8 IEDB class I immunogenicity (Calis 2013) — already installed

Nothing to download, no environment. A 20-value amino-acid scale times nine position
weights, implemented natively in `MHC1/Immunogenicity_scripts/calis_immunogenicity.py` and
verified two ways: identical to IEDB's own `immunogenicity_predictor.py` on 28,000 random
peptides across all lengths and mask choices, and reproducing all 11,781 peptide scores in
an IEDB website export (maximum difference 0.0).

---

## 6. NetMHCpan and the IEDB bundles

### 6.1 NetMHCpan 4.2 — class I binding affinity and eluted ligand

DTU academic-licence software: you must accept the licence on their site, so this step
cannot be scripted and `install_tools.sh` leaves it out.

1. Go to <https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/> and request the
   **Linux** version. You receive the software tarball and a link to the data tarball.
2. Extract both into `~/local_engines/netmhc/` and follow the `netMHCpan-4.2.readme`
   inside, which tells you to unpack the data file in the same folder.

```bash
mkdir -p ~/local_engines/netmhc && cd ~/local_engines/netmhc
tar -xzf /mnt/c/Users/<you>/Downloads/netMHCpan-4.2b.Linux.tar.gz    # reading from Windows is fine
cd netMHCpan-4.2 && tar -xzf /mnt/c/Users/<you>/Downloads/data.tar.gz
```

(Downloading in Windows and extracting *into* the Linux filesystem is the easy route under
WSL. Reading a tarball off `/mnt/c` once is fine — it is only the *installed* tools that
must not live there.)

**You do not need to edit the script.** DTU's `netMHCpan` ships with a placeholder `NMHOME`;
rather than modify your copy, the engine runs a patched copy from its scratch directory
whenever the placeholder is wrong. Editing it yourself also works.

The engine finds any folder with a file named `netMHCpan` beside a `Linux_x86_64/`
directory, at any depth, so the extra `netMHCpan-4.2/` level above is harmless.

### 6.2 NetCTLpan 1.1 and NetChop 3.1 — class I antigen processing

DTU programs that IEDB redistributes inside its next-generation T Cell Class I standalone
(free for academic and non-profit use; the bundle carries DTU's licence text). The engine
drives the bundle's own `tcell_mhci.py` in the bundle's own `uv` environment.

```bash
mkdir -p ~/local_engines/iedb && cd ~/local_engines/iedb
wget https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/IEDB_NG_TC1-v0.2-beta.tar.gz   # 743 MB
tar -xzf IEDB_NG_TC1-v0.2-beta.tar.gz
mv ng_tc1-v0.2-beta tc1        # NetCTL refuses a path over 77 characters
cd tc1 && uv sync && ./configure
```

Check the `README` in that downloads folder for the current file name if this one is gone.
Needs `tcsh` and `gawk`.

Two things to know before using them:

- **They need the whole protein.** Cleavage depends on the residues around a peptide, so
  these two apply when a protein is cut into fragments, not to a loose peptide list; the
  form switches them off in that mode. Putting them late in a funnel saves nothing.
- **NetCTLpan scores 8–11-mers only.** Other lengths get no score and the page says so.

### 6.3 NetMHCIIpan and CD4episcore — class II binding and immunogenicity

Same licence terms, 1.6 GB. Provides NetMHCIIpan EL + BA in versions 4.3 (default), 4.2 and
*IEDB default*, plus CD4episcore.

```bash
cd ~/local_engines/iedb
wget https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/IEDB_NG_TC2-v0.2.3-beta.tar.gz
tar -xzf IEDB_NG_TC2-v0.2.3-beta.tar.gz
mv ng_tc2-v0.2.3-beta tc2      # keep the path short here too
cd tc2 && uv sync
```

Needs `tcsh`, `gawk`, `perl` and `r-base-core`.

**Now apply §8.1** — in this release CD4episcore does not work without it.
`install_tools.sh --iedb` applies fix (a) for you.

---

## 7. Verify every model

```bash
cd ~/public_codes/epitope_filter_engine && ENGINE_MODELS_DIR=~/local_engines bash start.sh
```

In the browser, **Settings** shows *In use right now* — every path and environment that was
discovered — and a readiness checklist per model with the reason when something is missing.
Then press **Test run** for each key: it predicts a known peptide end to end through the
real tool.

| Key | Covers | Installed in |
|---|---|---|
| `netmhcpan` | NetMHCpan 4.2 | §6.1 |
| `bigmhc` | BigMHC EL + IM | §5.4 |
| `mixmhcpred` | MixMHCpred 3.0 | §5.1 |
| `prime` | PRIME 2.1 | §5.2 |
| `mhcflurry` | MHCflurry 2 | §5.5 |
| `iedb_bundle` | NetCTLpan, NetChop | §6.2 |
| `mixmhc2pred` | MixMHC2pred 2.1 | §5.3 |
| `iedb_bundle_ii` | NetMHCIIpan, CD4episcore | §6.3 |
| `immuscope` | ImmuScope-IM | §5.6 |
| `tlimmuno2` | TLimmuno2 | §5.7 |

The first BigMHC call loads its models from disk and can take minutes on a cold cache;
later calls are quick.

Two numbers to check by hand if you suspect a tool is wrong:

- **GILGFVFTL on HLA-A\*02:01** — what *Test run* uses for class I.
- **PKYVKQNTLKLATGM on DRB1\*01:01 must give CD4episcore 96.32.** If not, §8.1.

And the test suite, which needs **no predictors at all** (it stubs them), so it is also the
way to check the engine itself after a `git pull`:

```bash
conda run -n epitope-engine python test_engine.py      # config, server, API, funnel, workspace
conda run -n epitope-engine python test_consensus.py   # votes, agreement, coverage, priority score
conda run -n epitope-engine python test_host.py        # host-proteome search against brute force
```

All three must print `... passed`. `test_engine.py` skips its NetCTLpan and class II
sections when tcsh/gawk/perl/R are absent and says so. To exercise real BigMHC instead of
the stub: `BIGMHC_PYTHON=<env python> BIGMHC_HOME=<clone> conda run -n epitope-engine python test_engine.py`.

---

## 8. Troubleshooting the faults we actually hit

### 8.1 CD4episcore returns nothing: "loop of ufunc does not support argument 0 of type float"

Two independent causes in IEDB's class II bundle v0.2.3-beta, both seen on a
Core Ultra 9 / WSL machine. Run these from inside the class II bundle (`~/local_engines/iedb/tc2/`).

**(a) The R path is built from bytes**, giving `-R b/usr/bin/R`.
`install_tools.sh --iedb` applies this one automatically.

```bash
sed -i 's|R_path = subprocess.check_output("which R", shell=True).rstrip()$|R_path = subprocess.check_output("which R", shell=True).rstrip().decode()|' src/imscore.py
grep -n "R_path =" src/imscore.py        # must show exactly one .decode()
```

Run that `sed` **once only**.

**(b) NN-align's binaries segfault on very new Intel CPUs.** *The installer does this for you
when qemu is present: it runs one of the binaries, and only if it is killed does it wrap them.
If qemu is missing it says so. In that case:*

```bash
sudo apt install -y qemu-user
bash install_tools.sh --only iedb2      # quick: the bundle is already downloaded
```

*The rest of this entry is what that does, for a manual install.*

**(b) NN-align's binaries segfault on very new Intel CPUs.** The shipped
`src/nnalign/bin/Linux_x86_64/*` are statically linked for kernel 2.6 and crash
("Segmentation fault" before "Uploaded 0 peptides"). Run them under qemu:

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

### 8.1b NetCTL: "will not work with the current executable path ... at most 77 characters"

IEDB's `./configure` prints this when the bundle sits too deep. The installer avoids it by
renaming the extracted folder to `tc1`, and repairs an older install by renaming it for you
(re-running `bash start.sh --redo-setup` is enough). If you still see it, your home directory
is long: install with a shorter path, `bash start.sh --dir /opt/le`.

### 8.1c MixMHC2pred stays "not ready" after the PWMdef download

The zip contains a `PWMdef_Human/` folder, but the engine needs the `.txt` files **directly**
in `PWMdef/`. Unzip with `-j` (junk paths). An install made before this was fixed is repaired
in place by `bash start.sh --redo-setup` — it flattens the folder instead of downloading again.

### 8.1d "Could not find a version that satisfies the requirement pandas"

Seen while building the BigMHC or ImmuScope environment. `--index-url` *replaces* PyPI, so a
single `pip install torch --index-url <pytorch> numpy pandas ...` looks for pandas on the
PyTorch index, where it does not exist. Install torch on its own:

```bash
pip install numpy pandas h5py scikit-learn ruamel.yaml tqdm click
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

The installer does this correctly now; `bash start.sh --redo-setup` fixes an environment left
empty by the old one.

### 8.2 PRIME hangs

PRIME 2.1's C++ program never returns for **HLA-C\*07:06** (it spins at 100 % CPU), so the
engine skips that allele for PRIME with a note. PRIME is also run one allele per call under
a time limit, so any other allele that hangs is given up on instead of stalling the run.

### 8.3 tcsh prints "libtinfo.so.6: no version information available", or perl/R fail silently

A conda environment is active and these programs are loading conda's libraries instead of
the system's. **The engine already starts tcsh, gawk, perl and R with a cleaned
environment** (no `LD_LIBRARY_PATH`, no conda variables, no perl/R paths) — this is what
fixed CD4episcore failing with the ufunc error. If you still see it, `perl -v` and
`R --version` must work in a plain terminal outside conda.

This is the one real friction between conda and these tools, and it is handled; it is not a
reason to avoid conda.

### 8.4 WSL-specific

- **Everything is slow.** The tools are probably installed under `/mnt/c`. Move them into
  the Linux filesystem (§2.1).
- **"Permission denied" or "not executable"** running a DTU/IEDB program: same cause — a
  Windows drive does not carry the execute bit.
- **A run dies with no message, or the machine swaps hard.** WSL ran out of memory; raise
  `memory=` in `.wslconfig` (§2.1), or pick the **Gentle** performance profile.
- **`conda: command not found` after installing it.** Run `exec bash -l`, or reopen Ubuntu.
- **Clock skew after the laptop sleeps** can break TLS on downloads: `sudo hwclock -s`.

### 8.5 An allele gets no score

A model that does not know an allele casts no vote for it and the page names the skipped
alleles: MixMHCpred and PRIME have their own allele lists, BigMHC is human only, and HLA-DP
and mouse alleles have no population frequency data (they count in agreement but not in
coverage). An unknown allele makes NetMHCIIpan return an empty table with no error and
MixMHC2pred stop with a message and no output file, so the engine checks alleles first and
reports which it skipped.

### 8.6 A FASTA with stray characters

A `"` or `*` in a sequence makes NetCTLpan return nothing for the whole protein, so those
are replaced by `X` first, with a warning. Windows containing other stray characters are
skipped and the page says so.

### 8.7 A model fails mid-run

The run logs it (`✗ CD4episcore failed and is skipped: ...`), closes that stage with the
reason, warns on the results page and carries on with the other models. A class that
finished with a failed model is not saved for reuse, so pressing Run again retries it. Only
if *every* selected model fails does the run fail.

### 8.8 Port 5000 is taken

`PORT=5001 bash start.sh`.

---

## 9. Two optional data files

Both optional; the page says so when either is missing.

**Population HLA frequencies** — needed for population-weighted coverage and the priority
score's coverage term. Point `population.weights_csv` (Settings) at an HLA allele weightage
CSV, for example `india_hla_allele_weightage.csv` from the Epitope Workbench, or drop it
beside the engine folder; `ENGINE_POPULATION_CSV` also works. With no file, coverage falls
back to the share of the alleles *you* selected, and the page says so.

**Host proteome** — for the self-similarity check that flags candidates identical to, or one
residue from, a protein the host already makes. Settings → *Host proteome* → **Download**
(about 8 MB, kept in the workspace), or set `host_proteome.fasta` to your own FASTA (plain
or `.gz`). Set `host_proteome.organism` to `mouse` for mouse studies. A run whose rules must
*remove* self-peptides will not start until this is downloaded; *flag only* (the default)
just notes that it was not checked.

---

## 10. Environment variables

All optional; each can be set in `engine_config.json` instead.

| Variable | Effect |
|---|---|
| `ENGINE_MODELS_DIR` | Folder searched for the tools (default: the engine's parent folder) |
| `ENGINE_PYTHON` | Interpreter for the server, overriding the conda search in §4 |
| `ENGINE_CONFIG` | Path of the settings JSON (default: `engine_config.json` beside the engine) |
| `ENGINE_WORKSPACE` | Where runs, caches and timings go (default `<models folder>/engine_workspace`) |
| `ENGINE_POPULATION_CSV` | HLA frequency CSV (§9) |
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

## 11. Updating, removing, headless use

```bash
git -C ~/public_codes pull
conda env update -f ~/public_codes/epitope_filter_engine/environment.yml --prune   # only if it changed
```

Then re-run the §7 test suite. Your runs and settings are in `engine_workspace/`, outside
the code folder, so they survive. To remove everything:

```bash
rm -rf ~/public_codes ~/local_engines
conda env remove -n epitope-engine -y
for e in bigmhc mhcflurry immuscope tlimmuno2; do conda env remove -n $e -y; done
```

Nothing was installed system-wide except the apt packages in §2.2.

**No browser needed:**

```bash
# press "Save settings" in the page to get the JSON, then
conda run -n epitope-engine python run_engine.py epitope_settings.json -o results/
conda run -n epitope-engine python run_engine.py epitope_settings.json --fasta other.fa -o results/
```

The page's API takes the same JSON: `POST /api/pipeline/run`.

---

## 12. Licences

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
