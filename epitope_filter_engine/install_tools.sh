#!/usr/bin/env bash
# Install the open-licence predictors for the Epitope Filter Engine.
#
# Conda-based by default: each predictor gets its own conda environment, named so that
# the engine finds it on its own (discover.py prefers an environment whose name contains
# the tool's name). With --venv, or when conda is absent, a virtualenv is made inside
# each tool's folder instead.
#
# Covers the tools that can be fetched without accepting a licence form:
#   MixMHCpred 3.0, PRIME 2.1, MixMHC2pred 2.1, BigMHC, MHCflurry 2, ImmuScope-IM, TLimmuno2
# and, with --iedb, IEDB's two next-generation standalones (NetCTLpan/NetChop and
# NetMHCIIpan/CD4episcore; about 2.4 GB together).
#
# NetMHCpan 4.2 is NOT installed here: DTU require you to accept their licence on their
# own site, so it cannot be scripted. See INSTALL.md section 6.1.
#
# Every step is skipped when its result is already there, so re-running is safe.
# Tested on Ubuntu under WSL2. INSTALL.md section 2.2 lists the apt packages.
#
# Usage:
#   bash install_tools.sh [--dir DIR] [--iedb] [--only NAME,NAME] [--venv] [--dry-run]
#
#   --dir DIR      where the tool folders go (default ~/local_engines)
#   --iedb         also download the two IEDB bundles (large)
#   --only LIST    install only these, comma separated:
#                  mixmhcpred,prime,mixmhc2pred,bigmhc,mhcflurry,immuscope,tlimmuno2,iedb1,iedb2
#   --venv         use python3 -m venv instead of conda
#   --dry-run      print what would happen, change nothing

set -u -o pipefail

DIR="$HOME/local_engines"
WANT_IEDB=0
ONLY=""
DRY=0
FORCE_VENV=0
PY="${PYTHON:-python3}"
TORCH_CPU="https://download.pytorch.org/whl/cpu"

while [ $# -gt 0 ]; do
  case "$1" in
    --dir)     DIR="${2:?--dir needs a path}"; shift 2 ;;
    --iedb)    WANT_IEDB=1; shift ;;
    --only)    ONLY="${2:?--only needs a list}"; shift 2 ;;
    --venv)    FORCE_VENV=1; shift ;;
    --dry-run) DRY=1; shift ;;
    -h|--help) awk 'NR>1 && /^#/ {sub(/^# ?/,""); print; next} NR>1 {exit}' "$0"; exit 0 ;;
    *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
  esac
done

DIR="${DIR/#\~/$HOME}"
DONE=(); SKIPPED=(); FAILED=()

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
run()  { if [ "$DRY" = 1 ]; then printf '    $ %s\n' "$*"; else eval "$@"; fi; }

# ---------------------------------------------------------------- conda or venv

CONDA=""
if [ "$FORCE_VENV" = 0 ]; then
  if command -v mamba >/dev/null 2>&1; then CONDA="mamba"            # same CLI, much faster solver
  elif command -v conda >/dev/null 2>&1; then CONDA="conda"
  elif [ -n "${CONDA_EXE:-}" ] && [ -x "$CONDA_EXE" ]; then CONDA="$CONDA_EXE"
  fi
fi

env_exists() {                                 # env_exists <name>
  [ -n "$CONDA" ] || return 1
  "$CONDA" env list 2>/dev/null | awk '{print $1}' | grep -qx "$1"
}

# ensure_env <conda name> <venv dir> <python version> <import check> <pip args...>
#
# Creates the environment if missing and installs the packages if the import check
# fails. Returns the interpreter path in REPLY_PY.
ensure_env() {
  local name="$1" vdir="$2" pyver="$3" check="$4"; shift 4
  REPLY_PY=""
  if [ -n "$CONDA" ]; then
    if env_exists "$name"; then
      info "conda env exists: $name"
    else
      run "'$CONDA' create -y -q -n '$name' python='$pyver'" || return 1
    fi
    if [ "$DRY" != 1 ] && "$CONDA" run -n "$name" python -c "import $check" >/dev/null 2>&1; then
      info "packages already installed in $name ($check)"
    else
      run "'$CONDA' run -n '$name' pip install -q $*" || return 1
    fi
    [ "$DRY" = 1 ] || REPLY_PY="$("$CONDA" run -n "$name" python -c 'import sys; print(sys.executable)' 2>/dev/null)"
    info "environment: conda env '$name' (found by the engine automatically)"
  else
    if [ -x "$vdir/bin/python" ]; then
      info "environment exists: $vdir/"
    else
      run "'$PY' -m venv '$vdir'" || return 1
    fi
    if [ "$DRY" != 1 ] && "$vdir/bin/python" -c "import $check" >/dev/null 2>&1; then
      info "packages already installed in $vdir ($check)"
    else
      run "'$vdir/bin/pip' install -q --upgrade pip" || true
      run "'$vdir/bin/pip' install -q $*" || return 1
    fi
    REPLY_PY="$vdir/bin/python"
    info "environment: $vdir/"
  fi
  return 0
}

need() {                                       # need <cmd> <apt package>
  command -v "$1" >/dev/null 2>&1 && return 0
  info "missing: $1   (sudo apt install -y $2)"
  return 1
}

clone() {                                      # clone <url> <target>
  if [ -e "$2/.git" ]; then info "already cloned: $2/"; return 0; fi
  run "git clone --depth 1 '$1' '$2'"
}

fetch() {                                      # fetch <url> <output file>
  run "curl -fsSL --retry 3 -o '$2' '$1'"
}

wanted() {                                     # wanted <name>
  if [ -n "$ONLY" ]; then
    case ",$ONLY," in *",$1,"*) return 0 ;; *) return 1 ;; esac
  fi
  case "$1" in
    iedb1|iedb2) [ "$WANT_IEDB" = 1 ] ;;       # large: only with --iedb
    *) return 0 ;;
  esac
}

step() {                                       # step <name> <label> <function>
  if ! wanted "$1"; then SKIPPED+=("$1"); return 0; fi
  say "$2"
  if "$3"; then DONE+=("$1"); else FAILED+=("$1"); info "FAILED: $1"; fi
}

# ---------------------------------------------------------------- the tools

do_mixmhcpred() {
  # needs only numpy and pandas, which the engine's own environment already has
  clone https://github.com/GfellerLab/MixMHCpred.git "$DIR/gfeller/MixMHCpred"
}

do_prime() {
  need g++ g++ || info "PRIME needs g++ to rebuild lib/PRIME.x on first use"
  clone https://github.com/GfellerLab/PRIME.git "$DIR/gfeller/PRIME" || return 1
  info "PRIME calls MixMHCpred, so keep that one installed too"
}

do_mixmhc2pred() {
  local m="$DIR/gfeller/MixMHC2pred" sp
  clone https://github.com/GfellerLab/MixMHC2pred.git "$m" || return 1
  if [ -d "$m/PWMdef" ] && [ -n "$(ls -A "$m/PWMdef" 2>/dev/null)" ]; then
    info "PWMdef already present"; return 0
  fi
  need unzip unzip || return 1
  # the repository ships the program but not the allele definitions
  run "mkdir -p '$m/PWMdef'" || return 1
  for sp in Human Mouse; do
    fetch "https://mixmhc2pred.unil.ch/data/PWMdef/PWMdef_$sp.zip" "$m/PWMdef_$sp.zip" || return 1
    run "unzip -oq '$m/PWMdef_$sp.zip' -d '$m/PWMdef'" || return 1
    run "rm -f '$m/PWMdef_$sp.zip'"
  done
}

do_bigmhc() {
  clone https://github.com/KarchinLab/bigmhc.git "$DIR/bigmhc" || return 1
  ensure_env bigmhc "$DIR/bigmhc/.venv" 3.11 torch \
    "torch --index-url $TORCH_CPU numpy pandas psutil" || return 1
}

do_mhcflurry() {
  # MHCflurry has no repository to clone: the conda env (or MHCFLURRY_PYTHON) is all there is
  ensure_env mhcflurry "$DIR/mhcflurry-env" 3.11 mhcflurry \
    "mhcflurry --extra-index-url $TORCH_CPU" || return 1
  local py="$REPLY_PY" dl
  if [ "$DRY" = 1 ]; then
    info '$ mhcflurry-downloads fetch models_class1_presentation'
    return 0
  fi
  dl="$(dirname "$py")/mhcflurry-downloads"
  if "$dl" path models_class1_presentation >/dev/null 2>&1; then
    info "models_class1_presentation already downloaded"
  else
    run "'$dl' fetch models_class1_presentation" || return 1   # 160 MB
  fi
  [ -n "$CONDA" ] || info "no clone to find: export MHCFLURRY_PYTHON=$py"
}

do_immuscope() {
  local h="$DIR/ImmuScope" f
  clone https://github.com/shenlongchen/ImmuScope.git "$h" || return 1
  if [ -d "$h/weights/IM" ]; then
    info "weights already present"
  else
    fetch https://zenodo.org/api/records/14810445/files/ImmuScope-weights.tar.gz/content "$h/w.tgz" || return 1
    run "tar -xzf '$h/w.tgz' -C '$h'" || return 1      # 149.8 MB -> weights/IM, EL, CD4
    run "rm -f '$h/w.tgz'"
  fi
  if [ -f "$h/data/raw/pseudosequence.2023.dat" ]; then
    info "pseudosequence table already present"
  else
    # one file out of a 3 GB archive: stream it, never keep the archive
    run "mkdir -p '$h/data/raw'" || return 1
    run "curl -fsSL --retry 3 'https://zenodo.org/api/records/14810445/files/ImmuScope-data.tar.gz/content' | tar -xz -C '$h' --wildcards '*pseudosequence.2023.dat'" || return 1
    if [ "$DRY" != 1 ] && [ ! -f "$h/data/raw/pseudosequence.2023.dat" ]; then
      f="$(find "$h" -name pseudosequence.2023.dat -print -quit 2>/dev/null)"
      [ -n "$f" ] && run "mv '$f' '$h/data/raw/'"
    fi
  fi
  ensure_env immuscope "$h/.venv" 3.11 torch \
    "torch --index-url $TORCH_CPU numpy pandas h5py scikit-learn ruamel.yaml tqdm click" || return 1
}

do_tlimmuno2() {
  local h="$DIR/TLimmuno2" py="" c
  clone https://github.com/XSLiuLab/TLimmuno2.git "$h" || return 1
  # TensorFlow 2.15 (Keras 2) is required: Keras 3 cannot load these models. With conda the
  # environment is simply built on Python 3.11; with venv we must find such an interpreter.
  if [ -n "$CONDA" ]; then
    ensure_env tlimmuno2 "$h/.venv" 3.11 tensorflow \
      '"tensorflow-cpu==2.15.1" "numpy<2" pandas' || return 1
    return 0
  fi
  if "$PY" -c 'import sys; sys.exit(0 if sys.version_info[:2] <= (3,11) else 1)' 2>/dev/null; then
    py="$PY"
  else
    for c in python3.11 python3.10 python3.9; do
      command -v "$c" >/dev/null 2>&1 && { py="$c"; break; }
    done
  fi
  if [ -z "$py" ]; then
    info "tensorflow-cpu 2.15.1 has no wheel for $("$PY" -V 2>&1 | cut -d' ' -f2), and no python3.9-3.11 was found"
    info "install conda (INSTALL.md 3) or: sudo apt install -y python3.11-venv"
    return 1
  fi
  PY="$py" ensure_env tlimmuno2 "$h/.venv" 3.11 tensorflow \
    '"tensorflow-cpu==2.15.1" "numpy<2" pandas'
}

do_iedb1() {
  local d="$DIR/iedb" t=IEDB_NG_TC1-v0.2-beta.tar.gz
  need tcsh tcsh; need gawk gawk
  need uv uv || { info "uv: curl -LsSf https://astral.sh/uv/install.sh | sh"; return 1; }
  [ -d "$d/ng_tc1-v0.2-beta" ] && { info "already extracted"; return 0; }
  run "mkdir -p '$d'" || return 1
  run "curl -fL --retry 3 -o '$d/$t' 'https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/$t'" || return 1
  run "tar -xzf '$d/$t' -C '$d'" || return 1
  # the bundle keeps its own uv environment; the engine drives its tcell_mhci.py through it
  run "cd '$d/ng_tc1-v0.2-beta' && uv sync && ./configure" || return 1
  run "rm -f '$d/$t'"
}

do_iedb2() {
  local d="$DIR/iedb" t=IEDB_NG_TC2-v0.2.3-beta.tar.gz h="$DIR/iedb/ng_tc2-v0.2.3-beta"
  need tcsh tcsh; need gawk gawk; need perl perl; need R r-base-core
  need uv uv || { info "uv: curl -LsSf https://astral.sh/uv/install.sh | sh"; return 1; }
  if [ -d "$h" ]; then
    info "already extracted"
  else
    run "mkdir -p '$d'" || return 1
    run "curl -fL --retry 3 -o '$d/$t' 'https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/$t'" || return 1
    run "tar -xzf '$d/$t' -C '$d'" || return 1
    run "cd '$h' && uv sync" || return 1
    run "rm -f '$d/$t'"
  fi
  # INSTALL.md 8.1(a): src/imscore.py builds the R path from bytes, giving "-R b/usr/bin/R",
  # and CD4episcore then returns nothing. Applied once, guarded by the grep.
  if [ "$DRY" = 1 ]; then
    info '$ sed -i ... src/imscore.py ...   (CD4episcore R-path fix)'
  elif [ ! -f "$h/src/imscore.py" ]; then
    info "src/imscore.py not found: check the bundle version, then INSTALL.md 8.1"
  elif grep -q 'rstrip().decode()' "$h/src/imscore.py"; then
    info "CD4episcore R-path fix already applied"
  else
    sed -i 's|R_path = subprocess.check_output("which R", shell=True).rstrip()$|R_path = subprocess.check_output("which R", shell=True).rstrip().decode()|' "$h/src/imscore.py" \
      && info "applied the CD4episcore R-path fix (INSTALL.md 8.1a)"
  fi
  info "if CD4episcore segfaults on a very new Intel CPU, apply INSTALL.md 8.1(b) (qemu-user)"
}

# ---------------------------------------------------------------- run

say "Epitope Filter Engine: installing predictors into $DIR"
[ "$DRY" = 1 ] && info "(dry run: nothing will be changed)"
if [ -n "$CONDA" ]; then
  info "environments: conda ($CONDA), one per predictor, found by the engine by name"
else
  if [ "$FORCE_VENV" = 1 ]; then
    info "environments: python3 -m venv inside each tool folder (--venv)"
  else
    info "environments: python3 -m venv inside each tool folder"
    info "conda was not found; install Miniforge for the smoother path (INSTALL.md 3)"
  fi
fi
need git git   || exit 1
need curl curl || exit 1
run "mkdir -p '$DIR'" || exit 1

step mixmhcpred  "MixMHCpred 3.0 - class I presentation"            do_mixmhcpred
step prime       "PRIME 2.1 - class I immunogenicity"               do_prime
step mixmhc2pred "MixMHC2pred 2.1 - class II presentation"          do_mixmhc2pred
step bigmhc      "BigMHC - class I EL + IM (torch)"                 do_bigmhc
step mhcflurry   "MHCflurry 2 - class I affinity and processing"    do_mhcflurry
step immuscope   "ImmuScope-IM - class II immunogenicity (torch)"   do_immuscope
step tlimmuno2   "TLimmuno2 - class II immunogenicity (TF 2.15)"    do_tlimmuno2
step iedb1       "IEDB class I bundle - NetCTLpan, NetChop"         do_iedb1
step iedb2       "IEDB class II bundle - NetMHCIIpan, CD4episcore"  do_iedb2

say "Summary"
[ ${#DONE[@]}    -gt 0 ] && info "installed or already present: ${DONE[*]}"
[ ${#SKIPPED[@]} -gt 0 ] && info "not requested:                ${SKIPPED[*]}"
[ ${#FAILED[@]}  -gt 0 ] && info "failed:                       ${FAILED[*]}"
[ "$WANT_IEDB" = 0 ] && [ -z "$ONLY" ] && info "the IEDB bundles were not requested: add --iedb"

cat <<EOF

    NetMHCpan 4.2 is not installed by this script (DTU licence form): INSTALL.md 6.1

    Next:
      cd "$(cd "$(dirname "$0")" && pwd)" && ENGINE_MODELS_DIR="$DIR" bash start.sh
    then open http://127.0.0.1:5000 and press Test run for each model under Settings.
EOF

[ ${#FAILED[@]} -gt 0 ] && exit 1
exit 0
