#!/bin/bash
# Clone, run this, done.
#
#     git clone https://github.com/ydv-amit-ydv/public_codes.git
#     cd public_codes/immunoforge
#     bash start.sh
#
# The first run sets everything up — the engine's own environment, and every predictor that
# can be downloaded without accepting a licence form — and then opens the server. Every run
# after that just opens the server, in a second.
#
# NetMHCpan 4.2 is the one exception: DTU require you to accept their licence on their own
# site, so it cannot be scripted. The engine runs fine without it, and tells you at the end
# where to put it. See INSTALL.md section 6.1.
#
# Options (you need none of them):
#   --no-setup     just start; skip the first-run install
#   --setup-only   install, do not start the server
#   --light        skip the two big IEDB bundles (~2.4 GB) on the first run
#   --redo-setup   run the installer again even if it ran before
#   --dir DIR      where the predictors live (default ~/local_engines)
#   --dry-run      print what the first-run install would do, change nothing
#
# Environment: ENGINE_PYTHON, ENGINE_MODELS_DIR, PORT (see INSTALL.md section 10).
set -e
cd "$(dirname "$0")"

SETUP=auto; START=1; IEDB=--iedb; DRY=""; MODELS="${ENGINE_MODELS_DIR:-}"
while [ $# -gt 0 ]; do
  case "$1" in
    --no-setup)   SETUP=no; shift ;;
    --setup-only) START=0; SETUP=force; shift ;;
    --redo-setup) SETUP=force; shift ;;
    --light)      IEDB=""; shift ;;
    --dir)        MODELS="${2:?--dir needs a path}"; shift 2 ;;
    --dry-run)    DRY=--dry-run; SETUP=force; START=0; shift ;;
    -h|--help)    awk 'NR>1 && /^#/ {sub(/^# ?/,""); print; next} NR>1 {exit}' "$0"; exit 0 ;;
    *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
  esac
done

say() { printf '\n\033[1m==> %s\033[0m\n' "$*"; }

# uv and other user-level tools live here. Export it for every run, not just the first:
# the server runs the installer too (Settings -> Repair), and without this the IEDB
# bundles fail with "missing: uv" even though uv is installed.
[ -d "$HOME/.local/bin" ] && case ":$PATH:" in *":$HOME/.local/bin:"*) ;; *) export PATH="$HOME/.local/bin:$PATH" ;; esac

# ---------------------------------------------------------------- where the predictors live
#
# Discovery looks beside this folder unless ENGINE_MODELS_DIR says otherwise. A fresh clone
# has nothing beside it, so default to ~/local_engines — but if tools ARE already sitting
# next to the clone (the older layout), leave discovery alone and keep using them.
if [ -z "$MODELS" ]; then
  here="$(cd .. && pwd)"
  if ls -d "$here"/*/ 2>/dev/null | grep -qiE '/(netmhc|bigmhc|gfeller|iedb|immuscope|tlimmuno2)[^/]*/$' \
     || find "$here" -maxdepth 3 -name netMHCpan -print -quit 2>/dev/null | grep -q .; then
    MODELS="$here"                                   # tools are already beside the clone
  else
    MODELS="$HOME/local_engines"
  fi
fi
MODELS="${MODELS/#\~/$HOME}"
export ENGINE_MODELS_DIR="$MODELS"
STAMP="$MODELS/.immunoforge-setup-done"
[ ! -f "$STAMP" ] && [ -f "$MODELS/.engine-setup-done" ] && STAMP="$MODELS/.engine-setup-done"

# ---------------------------------------------------------------- system packages
#
# A handful of the predictors are shell programs: the IEDB bundles need tcsh and gawk,
# CD4episcore needs perl and R, PRIME is built with g++. They are the one part that needs
# root, so they are asked for once, plainly, rather than failing later inside a tool.
system_packages() {
  command -v apt-get >/dev/null 2>&1 || return 0           # not Debian/Ubuntu: INSTALL.md 2.2
  local want=() pkg cmd
  # qemu-user is for INSTALL.md 8.1(b): NN-align's old static binaries crash on recent Intel
  # CPUs, which makes CD4episcore fail with an unrelated-looking ufunc error. The installer
  # wraps them automatically, but only if qemu is here.
  for pair in "git:git" "curl:curl" "unzip:unzip" "g++:g++" "perl:perl" "tcsh:tcsh" "gawk:gawk" \
              "R:r-base-core" "qemu-x86_64:qemu-user"; do
    cmd="${pair%%:*}"; pkg="${pair##*:}"
    command -v "$cmd" >/dev/null 2>&1 || want+=("$pkg")
  done
  python3 -c "import venv" >/dev/null 2>&1 || want+=("python3-venv")
  [ ${#want[@]} -eq 0 ] && return 0

  say "System packages needed by some predictors"
  echo "    ${want[*]}"
  if [ "$(id -u)" = 0 ]; then
    apt-get update -qq && apt-get install -y -qq "${want[@]}" && return 0
  elif command -v sudo >/dev/null 2>&1; then
    if [ -t 0 ]; then
      printf '    Install them now with sudo? [Y/n] '
      read -r a </dev/tty || a=n
      case "$a" in [Nn]*) ;; *) sudo apt-get update -qq && sudo apt-get install -y "${want[@]}" && return 0 ;; esac
    elif sudo -n true 2>/dev/null; then                    # passwordless sudo, nobody to ask
      sudo apt-get update -qq && sudo apt-get install -y -qq "${want[@]}" && return 0
    fi
  fi
  echo "    Not installed. Run this yourself, then start again:"
  echo "      sudo apt install -y ${want[*]}"
  echo "    Carrying on: the models that need them will show as not ready in Settings."
  return 0
}

# uv runs the two IEDB bundles' own environments. It installs into the user's home, no root.
ensure_uv() {
  command -v uv >/dev/null 2>&1 && return 0
  [ -x "$HOME/.local/bin/uv" ] && { export PATH="$HOME/.local/bin:$PATH"; return 0; }
  [ -z "$IEDB" ] && return 0                               # --light: the bundles are skipped anyway
  say "Installing uv (runs the IEDB bundles' environments)"
  if curl -LsSf https://astral.sh/uv/install.sh | sh >/dev/null 2>&1; then
    export PATH="$HOME/.local/bin:$PATH"
    command -v uv >/dev/null 2>&1 && { echo "    ok"; return 0; }
  fi
  echo "    could not install uv; the two IEDB bundles will be skipped."
  echo "    Install it yourself and re-run with --redo-setup:"
  echo "      curl -LsSf https://astral.sh/uv/install.sh | sh"
  return 0
}

# ---------------------------------------------------------------- the engine's own interpreter
has_deps() { [ -x "$1" ] && "$1" -c "import flask, pandas" >/dev/null 2>&1; }

conda_bin() {
  if command -v mamba >/dev/null 2>&1; then echo mamba
  elif command -v conda >/dev/null 2>&1; then echo conda
  elif [ -n "${CONDA_EXE:-}" ] && [ -x "$CONDA_EXE" ]; then echo "$CONDA_EXE"
  fi
}

conda_env_python() {
  local root="" base
  [ -n "${CONDA_EXE:-}" ] && root="$(dirname "$(dirname "$CONDA_EXE")")"
  for base in "${CONDA_PREFIX_1:-}" "${CONDA_PREFIX:-}" "$root" \
              "$HOME/miniforge3" "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" "/opt/conda"; do
    [ -n "$base" ] || continue
    base="${base%%/envs/*}"
    # immunoforge is the current name; epitope-engine is what an install made
    # before the project was renamed, and is still used if it is what exists.
    for env in immunoforge epitope-engine; do
      [ -x "$base/envs/$env/bin/python" ] && { echo "$base/envs/$env/bin/python"; return; }
    done
  done
}

PY=""
if [ -n "${ENGINE_PYTHON:-}" ]; then                      # 1. explicit override
  PY="${ENGINE_PYTHON/#\~/$HOME}"
  has_deps "$PY" || { echo "ENGINE_PYTHON=$PY cannot import flask and pandas"; exit 1; }
elif [ -n "${CONDA_PREFIX:-}" ] && has_deps "$CONDA_PREFIX/bin/python"; then
  PY="$CONDA_PREFIX/bin/python"                           # 2. the active conda environment
else
  cand="$(conda_env_python || true)"                      # 3. the conda env epitope-engine
  if [ -n "$cand" ] && has_deps "$cand"; then
    PY="$cand"
  elif [ -n "$(conda_bin)" ] && [ "$SETUP" != no ] && [ -z "$DRY" ]; then
    say "Creating the conda environment 'immunoforge' (python, flask, pandas)"
    "$(conda_bin)" env create -q -f environment.yml >/dev/null 2>&1 || true
    cand="$(conda_env_python || true)"
    [ -n "$cand" ] && has_deps "$cand" && PY="$cand"
  fi
fi
if [ -z "$PY" ]; then                                     # 4. a local virtualenv
  if ! has_deps .venv/bin/python; then
    say "Creating .venv (flask + pandas)"
    python3 -m venv .venv || {
      echo "python3 -m venv failed. Install it (sudo apt install -y python3-venv) or install conda."
      exit 1
    }
    .venv/bin/pip install -q -r requirements.txt
  fi
  PY=.venv/bin/python
fi

# ---------------------------------------------------------------- first run: install the predictors
if [ "$SETUP" = force ] || { [ "$SETUP" = auto ] && [ ! -f "$STAMP" ]; }; then
  say "First run: installing the predictors into $MODELS"
  echo "    This downloads several GB and takes a while. It happens once."
  echo "    Later runs start in a second. Press Ctrl-C to skip (then: bash start.sh --no-setup)."
  [ -z "$DRY" ] && { system_packages; ensure_uv; }
  if ENGINE_FROM_START=1 bash install_tools.sh --dir "$MODELS" $IEDB $DRY; then
    [ -z "$DRY" ] && date -u +"%Y-%m-%dT%H:%M:%SZ" > "$STAMP"
  else
    echo
    echo "    Some predictors did not install. The engine still runs with the rest;"
    echo "    Settings shows what is missing. To try again: bash start.sh --redo-setup"
    [ -z "$DRY" ] && date -u +"%Y-%m-%dT%H:%M:%SZ" > "$STAMP"
  fi
  if [ -z "$DRY" ] && ! find "$MODELS" -maxdepth 4 -name netMHCpan -print -quit 2>/dev/null | grep -q .; then
    say "One thing left for you, by hand"
    echo "    NetMHCpan 4.2 needs DTU's licence form, so it cannot be downloaded by a script:"
    echo "      https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/  (choose Linux)"
    echo "    Extract it anywhere under $MODELS and the engine finds it. Everything else is ready."
  fi
fi

[ "$START" = 1 ] || exit 0

say "Starting on http://127.0.0.1:${PORT:-5000}"
echo "    predictors: $ENGINE_MODELS_DIR"
echo "    interpreter: $PY"
exec "$PY" server.py
