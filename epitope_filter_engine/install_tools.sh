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
#                  mixmhcpred,prime,mixmhc2pred,bigmhc,mhcflurry,immuscope,tlimmuno2,screen,iedb1,iedb2
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
    [ "$DRY" = 1 ] || REPLY_PY="$("$CONDA" run -n "$name" python -c 'import sys; print(sys.executable)' 2>/dev/null)"
    info "environment: conda env '$name' (found by the engine automatically)"
  else
    if [ -x "$vdir/bin/python" ]; then
      info "environment exists: $vdir/"
    else
      run "'$PY' -m venv '$vdir'" || return 1
      run "'$vdir/bin/pip' install -q --upgrade pip" || true
    fi
    REPLY_PY="$vdir/bin/python"
    info "environment: $vdir/"
  fi
  [ $# -gt 0 ] || return 0
  env_pip "$name" "$vdir" "$check" "$@"
}

env_has() {                                    # env_has <name> <vdir> <module>
  [ "$DRY" = 1 ] && return 1
  if [ -n "$CONDA" ]; then
    "$CONDA" run -n "$1" python -c "import $3" >/dev/null 2>&1
  else
    [ -x "$2/bin/python" ] && "$2/bin/python" -c "import $3" >/dev/null 2>&1
  fi
}

env_pip() {                                    # env_pip <name> <vdir> <import check> <pip args...>
  local name="$1" vdir="$2" check="$3"; shift 3
  if env_has "$name" "$vdir" "$check"; then
    info "already installed: $check"
    return 0
  fi
  if [ -n "$CONDA" ]; then
    run "'$CONDA' run -n '$name' pip install -q $*"
  else
    run "'$vdir/bin/pip' install -q $*"
  fi
}

# CPU torch, installed ON ITS OWN.
#
# --index-url REPLACES PyPI rather than adding to it, so a single
#   pip install torch --index-url <pytorch cpu> numpy pandas h5py ...
# sends pip to the PyTorch index for pandas too and dies with
#   "Could not find a version that satisfies the requirement pandas".
# That is exactly what happened, so torch gets its own call and everything else
# comes from PyPI normally.
ensure_torch() {                               # ensure_torch <name> <vdir>
  env_pip "$1" "$2" torch "torch --index-url $TORCH_CPU"
}

# Echo a python3.9-3.11, or nothing. Several pinned wheels (pandas 1.5.x, TensorFlow 2.15)
# stop at 3.11, and building them from source on a newer interpreter fails.
older_python() {
  local c
  if "$PY" -c 'import sys; sys.exit(0 if sys.version_info[:2] <= (3,11) else 1)' 2>/dev/null; then
    echo "$PY"; return 0
  fi
  for c in python3.11 python3.10 python3.9; do
    command -v "$c" >/dev/null 2>&1 && { echo "$c"; return 0; }
  done
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
  clone https://github.com/GfellerLab/MixMHCpred.git "$DIR/gfeller/MixMHCpred" || return 1
  return 0
}

do_prime() {
  need g++ g++ || info "PRIME needs g++ to rebuild lib/PRIME.x on first use"
  clone https://github.com/GfellerLab/PRIME.git "$DIR/gfeller/PRIME" || return 1
  info "PRIME calls MixMHCpred, so keep that one installed too"
  return 0
}

do_mixmhc2pred() {
  local m="$DIR/gfeller/MixMHC2pred" sp
  clone https://github.com/GfellerLab/MixMHC2pred.git "$m" || return 1
  # The engine needs the definition .txt files DIRECTLY in PWMdef/. The zip holds them in a
  # PWMdef_Human/ folder, so unzip must junk paths (-j); an install made before that fix has
  # them one level too deep and is flattened here instead of being re-downloaded.
  if [ "$DRY" != 1 ] && [ -d "$m/PWMdef" ] && ! ls "$m"/PWMdef/*.txt >/dev/null 2>&1; then
    if find "$m/PWMdef" -mindepth 2 -name '*.txt' -print -quit 2>/dev/null | grep -q .; then
      info "flattening PWMdef (the files were one folder too deep)"
      find "$m/PWMdef" -mindepth 2 -type f -exec mv -f {} "$m/PWMdef/" \; 2>/dev/null
      find "$m/PWMdef" -mindepth 1 -type d -empty -delete 2>/dev/null
    fi
  fi
  if [ "$DRY" != 1 ] && ls "$m"/PWMdef/*.txt >/dev/null 2>&1; then
    info "PWMdef already present ($(ls "$m"/PWMdef/*.txt | wc -l) allele files)"; return 0
  fi
  need unzip unzip || return 1
  # the repository ships the program but not the allele definitions
  run "mkdir -p '$m/PWMdef'" || return 1
  for sp in Human Mouse; do
    fetch "https://mixmhc2pred.unil.ch/data/PWMdef/PWMdef_$sp.zip" "$m/PWMdef_$sp.zip" || return 1
    run "unzip -ojq '$m/PWMdef_$sp.zip' -d '$m/PWMdef'" || return 1
    run "rm -f '$m/PWMdef_$sp.zip'"
  done
  return 0
}

do_bigmhc() {
  clone https://github.com/KarchinLab/bigmhc.git "$DIR/bigmhc" || return 1
  ensure_env bigmhc "$DIR/bigmhc/.venv" 3.11 pandas numpy pandas psutil || return 1
  ensure_torch bigmhc "$DIR/bigmhc/.venv" || return 1
  return 0
}

do_mhcflurry() {
  # MHCflurry has no repository to clone: the conda env (or MHCFLURRY_PYTHON) is all there is
  # --extra-index-url ADDS to PyPI (unlike --index-url), so mhcflurry and its torch
  # dependency can be resolved in one call
  ensure_env mhcflurry "$DIR/mhcflurry-env" 3.11 mhcflurry "mhcflurry --extra-index-url $TORCH_CPU" || return 1
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
  return 0
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
  ensure_env immuscope "$h/.venv" 3.11 h5py numpy pandas h5py scikit-learn ruamel.yaml tqdm click || return 1
  ensure_torch immuscope "$h/.venv" || return 1
  return 0
}

do_tlimmuno2() {
  local h="$DIR/TLimmuno2" py="" c
  clone https://github.com/XSLiuLab/TLimmuno2.git "$h" || return 1
  # TensorFlow 2.15 (Keras 2) is required: Keras 3 cannot load these models. With conda the
  # environment is simply built on Python 3.11; with venv we must find such an interpreter.
  if [ -n "$CONDA" ]; then
    ensure_env tlimmuno2 "$h/.venv" 3.11 tensorflow '"tensorflow-cpu==2.15.1" "numpy<2" pandas' || return 1
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
  PY="$py" ensure_env tlimmuno2 "$h/.venv" 3.11 tensorflow '"tensorflow-cpu==2.15.1" "numpy<2" pandas' || return 1
  return 0
}

# NetCTL refuses to run if its own executable path is longer than 77 characters, and
# "<models>/iedb/ng_tc1-v0.2-beta/method/nctl/nctle/netCTL.1" alone is 57 of them, so a
# normal home directory pushes it over and ./configure warns that NetCTL will not work.
# The engine finds a bundle by its files, not its folder name, so the extracted directory
# is renamed to something short. shorten_bundle also repairs an install done before this.
shorten_bundle() {                             # shorten_bundle <dir> <glob> <short name>
  local d="$1" pat="$2" short="$3" old
  [ -d "$d/$short" ] && return 0
  old="$(find "$d" -maxdepth 1 -name "$pat" -type d -print -quit 2>/dev/null)"
  [ -n "$old" ] || return 1
  run "mv '$old' '$d/$short'" || return 1
  info "renamed $(basename "$old") -> $short (NetCTL needs a path under 78 characters)"
  # a uv/virtualenv records absolute paths in pyvenv.cfg and in every script shebang, so it
  # does not survive its folder being renamed: drop it and let uv sync build it again (seconds)
  if [ "$DRY" != 1 ] && [ -d "$d/$short/.venv" ]; then
    rm -rf "$d/$short/.venv"
    info "removed the bundle's old .venv so it is rebuilt for the new path"
  fi
}

do_screen() {
  # Safety and function screens: ToxinPred3 (toxicity), AlgPred2 (allergenicity),
  # IFNepitope2 (IFN-gamma induction).
  #
  # ToxinPred2 is deliberately NOT installed. It is the previous version of ToxinPred3,
  # trained on whole proteins, and running it over epitopes misfires: it called GILGFVFTL,
  # a benign and very well characterised influenza epitope, a toxin. ToxinPred3 accepts
  # proteins as well as peptides, so one model covers both.
  #
  # The version pins are not cosmetic, they are what makes these run at all:
  #   * scikit-learn 1.2.x — the shipped models are pickled trees from before scikit-learn
  #     1.3 added missing_go_to_left to the node dtype; a newer one refuses to unpickle them.
  #     It is installed LAST because the other packages pull a newer one in.
  #   * pandas 1.5.3 — ToxinPred3 writes a temporary FASTA with to_csv(sep="\n"), which
  #     pandas 2+ rejects as a bad delimiter.
  #   * IFNepitope2 hard-pins pandas==1.4.0 and scikit-learn==1.2.1, neither of which has a
  #     wheel here, so it goes in with --no-deps; 1.5.3 and 1.2.2 run it fine.
  local ipy="$PY"
  if [ -z "$CONDA" ]; then
    ipy="$(older_python)"
    if [ -z "$ipy" ]; then
      info "these tools pin pandas 1.5.x and scikit-learn 1.2.x, which have no wheel for"
      info "$("$PY" -V 2>&1). Install an older interpreter (sudo apt install -y python3.11-venv)"
      info "and run again, or install conda, which the installer uses to make a 3.11 environment."
      return 1
    fi
  fi
  PY="$ipy" ensure_env vaxscreen "$DIR/vaxscreen-env" 3.11 pandas '"numpy<2" "pandas==1.5.3" joblib tqdm' || return 1   # tqdm: IFNepitope2 needs it and --no-deps skips it
  local py="$REPLY_PY"
  [ "$DRY" = 1 ] || [ -n "$py" ] || return 1
  env_pip vaxscreen "$DIR/vaxscreen-env" toxinpred3 toxinpred3 algpred2 || return 1
  if [ "$DRY" = 1 ]; then
    info '$ pip install --no-deps ifnepitope2'
  elif ! "$py" -c "import importlib.util,sys; sys.exit(0 if importlib.util.find_spec('ifnepitope2') else 1)" 2>/dev/null; then
    run "'$py' -m pip install -q --no-deps ifnepitope2" || info "IFNepitope2 did not install; the others still work"
  else
    info "already installed: ifnepitope2"
  fi
  # last word on scikit-learn, whatever the packages above pulled in
  run "'${py:-python3}' -m pip install -q 'scikit-learn==1.2.2'" || return 1
  info "HemoPI2 (haemolysis) is not installed by default: it needs torch. Add it with"
  info "  '${py:-<screening python>}' -m pip install hemopi2 torch --index-url $TORCH_CPU"
  return 0
}

do_iedb1() {
  local d="$DIR/iedb" t=IEDB_NG_TC1-v0.2-beta.tar.gz h="$DIR/iedb/tc1"
  need tcsh tcsh; need gawk gawk
  need uv uv || { info "uv: curl -LsSf https://astral.sh/uv/install.sh | sh"; return 1; }
  if [ ! -d "$h" ]; then
    if shorten_bundle "$d" 'ng_tc1-*' tc1; then             # repair an earlier install
      run "cd '$h' && uv sync && ./configure" || return 1
      return 0
    fi
    run "mkdir -p '$d'" || return 1
    run "curl -fL --retry 3 -o '$d/$t' 'https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/$t'" || return 1
    run "tar -xzf '$d/$t' -C '$d'" || return 1
    shorten_bundle "$d" 'ng_tc1-*' tc1 || return 1
    run "rm -f '$d/$t'"
  else
    info "already extracted"
  fi
  # the bundle keeps its own uv environment; the engine drives its tcell_mhci.py through it
  if [ "$DRY" = 1 ] || [ ! -x "$h/.venv/bin/python" ]; then
    run "cd '$h' && uv sync && ./configure" || return 1
  else
    info "bundle environment already built"
  fi
  if [ "$DRY" != 1 ] && [ ${#h} -gt 20 ]; then
    local probe="$h/method/nctl/nctle/netCTL.1"
    [ ${#probe} -gt 77 ] && info "WARNING: $probe is ${#probe} characters; NetCTL needs 77 or fewer. Use --dir with a shorter path."
  fi
  return 0
}

do_iedb2() {
  local d="$DIR/iedb" t=IEDB_NG_TC2-v0.2.3-beta.tar.gz h="$DIR/iedb/tc2"
  need tcsh tcsh; need gawk gawk; need perl perl; need R r-base-core
  need uv uv || { info "uv: curl -LsSf https://astral.sh/uv/install.sh | sh"; return 1; }
  if [ -d "$h" ]; then
    info "already extracted"
  elif shorten_bundle "$d" 'ng_tc2-*' tc2; then             # repair an earlier install
    :
  else
    run "mkdir -p '$d'" || return 1
    run "curl -fL --retry 3 -o '$d/$t' 'https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/$t'" || return 1
    run "tar -xzf '$d/$t' -C '$d'" || return 1
    shorten_bundle "$d" 'ng_tc2-*' tc2 || return 1
    run "rm -f '$d/$t'"
  fi
  if [ "$DRY" = 1 ] || [ ! -x "$h/.venv/bin/python" ]; then
    run "cd '$h' && uv sync" || return 1
  else
    info "bundle environment already built"
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
  return 0
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
step screen      "Safety screens - toxicity, allergenicity, IFN-gamma" do_screen
step iedb1       "IEDB class I bundle - NetCTLpan, NetChop"         do_iedb1
step iedb2       "IEDB class II bundle - NetMHCIIpan, CD4episcore"  do_iedb2

say "Summary"
[ ${#DONE[@]}    -gt 0 ] && info "installed or already present: ${DONE[*]}"
[ ${#SKIPPED[@]} -gt 0 ] && info "not requested:                ${SKIPPED[*]}"
[ ${#FAILED[@]}  -gt 0 ] && info "failed:                       ${FAILED[*]}"
[ "$WANT_IEDB" = 0 ] && [ -z "$ONLY" ] && info "the IEDB bundles were not requested: add --iedb"

if [ -z "${ENGINE_FROM_START:-}" ]; then                 # start.sh says this part itself
cat <<EOF

    NetMHCpan 4.2 is not installed by this script (DTU licence form): INSTALL.md 6.1

    Next:
      cd "$(cd "$(dirname "$0")" && pwd)" && ENGINE_MODELS_DIR="$DIR" bash start.sh
    then open http://127.0.0.1:5000 and press Test run for each model under Settings.
EOF
fi

[ ${#FAILED[@]} -gt 0 ] && exit 1
exit 0
