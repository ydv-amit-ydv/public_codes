#!/usr/bin/env bash
# Install the open-licence predictors for the Epitope Filter Engine.
#
# Covers the tools that can be fetched without accepting a licence form:
#   MixMHCpred 3.0, PRIME 2.1, MixMHC2pred 2.1, BigMHC, MHCflurry 2, ImmuScope-IM, TLimmuno2
# and, with --iedb, IEDB's two next-generation standalones (NetCTLpan/NetChop and
# NetMHCIIpan/CD4episcore; about 2.4 GB together).
#
# NetMHCpan 4.2 is NOT installed here: DTU require you to accept their licence on their
# own site, so it cannot be scripted. See INSTALL.md section 4.1.
#
# Every step is skipped when its result is already there, so re-running is safe.
# Nothing is installed system-wide; INSTALL.md section 2 lists the apt packages.
#
# Usage:
#   bash install_tools.sh [--dir DIR] [--iedb] [--only NAME,NAME] [--dry-run]
#
#   --dir DIR      where the tools go (default ~/local_engines)
#   --iedb         also download the two IEDB bundles (large)
#   --only LIST    install only these, comma separated:
#                  mixmhcpred,prime,mixmhc2pred,bigmhc,mhcflurry,immuscope,tlimmuno2,iedb1,iedb2
#   --dry-run      print what would happen, change nothing

set -u -o pipefail

DIR="$HOME/local_engines"
WANT_IEDB=0
ONLY=""
DRY=0
PY="${PYTHON:-python3}"
TORCH_CPU="https://download.pytorch.org/whl/cpu"

while [ $# -gt 0 ]; do
  case "$1" in
    --dir)     DIR="${2:?--dir needs a path}"; shift 2 ;;
    --iedb)    WANT_IEDB=1; shift ;;
    --only)    ONLY="${2:?--only needs a list}"; shift 2 ;;
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

wanted() {                                     # wanted <name>
  if [ -n "$ONLY" ]; then
    case ",$ONLY," in *",$1,"*) return 0 ;; *) return 1 ;; esac
  fi
  case "$1" in
    iedb1|iedb2) [ "$WANT_IEDB" = 1 ] ;;       # large: only with --iedb
    *) return 0 ;;
  esac
}

need() {                                       # need <cmd> <apt package>
  command -v "$1" >/dev/null 2>&1 && return 0
  info "missing: $1   (sudo apt install -y $2)"
  return 1
}

mkvenv() {                                     # mkvenv <dir> [interpreter]
  local d="$1" py="${2:-$PY}"
  if [ -x "$d/bin/python" ]; then info "environment exists: $d/"; return 0; fi
  run "'$py' -m venv '$d'"
}

pipin() {                                      # pipin <dir> <stamp> <pip args...>
  local d="$1" stamp="$2"; shift 2
  if [ -f "$d/.installed-$stamp" ]; then info "packages already installed: $stamp"; return 0; fi
  run "'$d/bin/pip' install -q $*" || return 1
  [ "$DRY" = 1 ] || : > "$d/.installed-$stamp"
}

clone() {                                      # clone <url> <target>
  if [ -e "$2/.git" ]; then info "already cloned: $2/"; return 0; fi
  run "git clone --depth 1 '$1' '$2'"
}

fetch() {                                      # fetch <url> <output file>
  run "curl -fsSL --retry 3 -o '$2' '$1'"
}

step() {                                       # step <name> <label> <function>
  if ! wanted "$1"; then SKIPPED+=("$1"); return 0; fi
  say "$2"
  if "$3"; then DONE+=("$1"); else FAILED+=("$1"); info "FAILED: $1"; fi
}

# ---------------------------------------------------------------- the tools

do_mixmhcpred() {
  clone https://github.com/GfellerLab/MixMHCpred.git "$DIR/gfeller/MixMHCpred"
}

do_prime() {
  need g++ g++ || info "PRIME needs g++ to rebuild lib/PRIME.x on first use"
  clone https://github.com/GfellerLab/PRIME.git "$DIR/gfeller/PRIME" || return 1
  info "PRIME calls MixMHCpred, so keep that installed too"
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
  mkvenv "$DIR/bigmhc/.venv" || return 1
  pipin  "$DIR/bigmhc/.venv" torch "torch --index-url $TORCH_CPU" || return 1
  pipin  "$DIR/bigmhc/.venv" deps numpy pandas psutil
}

do_mhcflurry() {
  local e="$DIR/mhcflurry-env"
  mkvenv "$e" || return 1
  pipin  "$e" mhcflurry "mhcflurry --extra-index-url $TORCH_CPU" || return 1
  if [ "$DRY" = 1 ]; then
    info '$ mhcflurry-downloads fetch models_class1_presentation'
  elif "$e/bin/mhcflurry-downloads" path models_class1_presentation >/dev/null 2>&1; then
    info "models_class1_presentation already downloaded"
  else
    run "'$e/bin/mhcflurry-downloads' fetch models_class1_presentation" || return 1
  fi
  info "found automatically; MHCFLURRY_PYTHON=$e/bin/python overrides"
}

do_immuscope() {
  local h="$DIR/ImmuScope" f
  clone https://github.com/shenlongchen/ImmuScope.git "$h" || return 1
  if [ -d "$h/weights/IM" ]; then
    info "weights already present"
  else
    fetch https://zenodo.org/api/records/14810445/files/ImmuScope-weights.tar.gz/content "$h/w.tgz" || return 1
    run "tar -xzf '$h/w.tgz' -C '$h'" || return 1
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
  mkvenv "$h/.venv" || return 1
  pipin  "$h/.venv" torch "torch --index-url $TORCH_CPU" || return 1
  pipin  "$h/.venv" deps numpy pandas h5py scikit-learn ruamel.yaml tqdm click
}

do_tlimmuno2() {
  local h="$DIR/TLimmuno2" py="" c
  clone https://github.com/XSLiuLab/TLimmuno2.git "$h" || return 1
  # TensorFlow 2.15 (Keras 2) is required: Keras 3 cannot load these models. TF 2.15 has
  # no wheel above Python 3.11, so prefer an older interpreter when one is installed.
  if "$PY" -c 'import sys; sys.exit(0 if sys.version_info[:2] <= (3,11) else 1)' 2>/dev/null; then
    py="$PY"
  else
    for c in python3.11 python3.10 python3.9; do
      command -v "$c" >/dev/null 2>&1 && { py="$c"; break; }
    done
  fi
  if [ -z "$py" ]; then
    info "no python3.9-3.11 found, and tensorflow-cpu 2.15.1 has no wheel for $("$PY" -V 2>&1 | cut -d' ' -f2)"
    info "install one (sudo apt install -y python3.11-venv), then: bash $0 --only tlimmuno2"
    return 1
  fi
  mkvenv "$h/.venv" "$py" || return 1
  pipin  "$h/.venv" tf '"tensorflow-cpu==2.15.1" "numpy<2" pandas'
}

do_iedb1() {
  local d="$DIR/iedb" t=IEDB_NG_TC1-v0.2-beta.tar.gz
  need tcsh tcsh; need gawk gawk
  need uv uv || { info "uv: curl -LsSf https://astral.sh/uv/install.sh | sh"; return 1; }
  [ -d "$d/ng_tc1-v0.2-beta" ] && { info "already extracted"; return 0; }
  run "mkdir -p '$d'" || return 1
  run "curl -fL --retry 3 -o '$d/$t' 'https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/$t'" || return 1
  run "tar -xzf '$d/$t' -C '$d'" || return 1
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
  # INSTALL.md 7.1(a): src/imscore.py builds the R path from bytes, giving "-R b/usr/bin/R",
  # and CD4episcore then returns nothing. Applied once, guarded by the grep.
  if [ "$DRY" = 1 ]; then
    info '$ sed -i ... src/imscore.py ...   (CD4episcore R-path fix)'
  elif [ ! -f "$h/src/imscore.py" ]; then
    info "src/imscore.py not found: check the bundle version, then INSTALL.md 7.1"
  elif grep -q 'rstrip().decode()' "$h/src/imscore.py"; then
    info "CD4episcore R-path fix already applied"
  else
    sed -i 's|R_path = subprocess.check_output("which R", shell=True).rstrip()$|R_path = subprocess.check_output("which R", shell=True).rstrip().decode()|' "$h/src/imscore.py" \
      && info "applied the CD4episcore R-path fix (INSTALL.md 7.1a)"
  fi
  info "if CD4episcore still segfaults on a very new Intel CPU, apply INSTALL.md 7.1(b) (qemu-user)"
}

# ---------------------------------------------------------------- run

say "Epitope Filter Engine: installing predictors into $DIR"
[ "$DRY" = 1 ] && info "(dry run: nothing will be changed)"
need git git  || exit 1
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

    NetMHCpan 4.2 is not installed by this script (DTU licence form): INSTALL.md 4.1

    Next:
      cd "$(cd "$(dirname "$0")" && pwd)" && ENGINE_MODELS_DIR="$DIR" bash start.sh
    then open http://127.0.0.1:5000 and press Test run for each model under Settings.
EOF

[ ${#FAILED[@]} -gt 0 ] && exit 1
exit 0
