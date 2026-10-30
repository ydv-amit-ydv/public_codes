#!/usr/bin/env bash
# circopt - one-shot setup and launch for WSL / Linux.
#
# Works out what state the machine is in and does the right thing:
#   * no clone yet            -> clones it
#   * a clone in the wrong    -> rescues run results, removes the redundant nested copy,
#     place (nested)             keeps the real one, updates it
#   * a clone already fine    -> just updates it
#   * conda env missing       -> creates it and installs the dependencies
# then runs the self-test and starts the web server.
#
# It never deletes anything that holds work you have not pushed: uncommitted changes,
# unpushed commits or a nested copy that is actually tracked all stop the delete and get
# moved to a backup instead.
#
#   bash setup_and_run.sh                 # ask before changing anything
#   bash setup_and_run.sh --yes           # no questions
#   bash setup_and_run.sh --no-server     # set everything up but do not start the server
#   bash setup_and_run.sh --port 8080     # serve on another port
#   bash setup_and_run.sh --dry-run       # say what it would do, change nothing

set -uo pipefail

REPO_URL="https://github.com/ydv-amit-ydv/public_codes.git"
REPO="$HOME/public_codes"
PKG_DIR="$REPO/circrna_opt"
RUNS_DIR="$HOME/circopt_runs"
ENV_NAME="circopt"
PORT=5000
ASSUME_YES=0
START_SERVER=1
DRY_RUN=0
STAMP="$(date +%Y%m%d-%H%M%S)"
BACKUP="$HOME/circopt_backup-$STAMP"

while [ $# -gt 0 ]; do
  case "$1" in
    --yes|-y) ASSUME_YES=1 ;;
    --no-server) START_SERVER=0 ;;
    --dry-run) DRY_RUN=1 ;;
    --port) PORT="${2:?--port needs a number}"; shift ;;
    --env) ENV_NAME="${2:?--env needs a name}"; shift ;;
    --repo) REPO="${2:?--repo needs a path}"; PKG_DIR="$REPO/circrna_opt"; shift ;;
    -h|--help) sed -n '2,24p' "$0"; exit 0 ;;
    *) echo "unknown option: $1 (try --help)" >&2; exit 2 ;;
  esac
  shift
done

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
info() { printf '    %s\n' "$*"; }
warn() { printf '\033[33m    ! %s\033[0m\n' "$*"; }
die()  { printf '\n\033[31m!! %s\033[0m\n' "$*" >&2; exit 1; }
run()  { if [ "$DRY_RUN" = 1 ]; then info "[dry-run] $*"; else "$@"; fi; }

confirm() {
  [ "$ASSUME_YES" = 1 ] && return 0
  [ "$DRY_RUN" = 1 ] && return 0
  printf '    %s [y/N] ' "$1"
  read -r reply </dev/tty || return 1
  case "$reply" in [yY]*) return 0 ;; *) return 1 ;; esac
}

# never operate from inside a directory we may delete
cd "$HOME" || die "cannot cd to $HOME"
command -v git >/dev/null 2>&1 || die "git is not installed:  sudo apt install git"

# ---------------------------------------------------------------- 1. the repo
say "Checking $REPO"
is_repo() { [ -d "$1/.git" ] && git -C "$1" rev-parse --is-inside-work-tree >/dev/null 2>&1; }

if [ ! -e "$REPO" ]; then
  info "nothing there yet - cloning"
  run git clone "$REPO_URL" "$REPO" || die "clone failed"
elif is_repo "$REPO"; then
  remote="$(git -C "$REPO" remote get-url origin 2>/dev/null || true)"
  case "$remote" in
    *ydv-amit-ydv/public_codes*) info "already a clone of public_codes" ;;
    "") warn "this git repo has no 'origin' remote - leaving it alone" ;;
    *) die "$REPO is a clone of $remote, not public_codes. Move it aside first." ;;
  esac
elif [ -z "$(ls -A "$REPO" 2>/dev/null)" ]; then
  info "empty directory - cloning into it"
  run rmdir "$REPO" && run git clone "$REPO_URL" "$REPO" || die "clone failed"
else
  die "$REPO exists, is not a git repo and is not empty. Move it aside, then re-run."
fi

# ------------------------------------------------- 2. the redundant nested clone
NESTED="$REPO/public_codes"
if [ -e "$NESTED" ]; then
  say "Found a nested copy at $NESTED"

  # rescue run results first, whatever happens next
  if [ -d "$NESTED/circrna_opt/runs" ] && [ -n "$(ls -A "$NESTED/circrna_opt/runs" 2>/dev/null)" ]; then
    info "rescuing run results -> $RUNS_DIR"
    run mkdir -p "$RUNS_DIR"
    run cp -rn "$NESTED/circrna_opt/runs/." "$RUNS_DIR/" 2>/dev/null || true
    [ "$DRY_RUN" = 1 ] || info "kept: $(ls -A "$RUNS_DIR" | wc -l) run(s) in $RUNS_DIR"
  fi

  keep_reason=""
  if git -C "$REPO" ls-files --error-unmatch public_codes >/dev/null 2>&1; then
    keep_reason="it is tracked by the outer repo (deleting would stage a mass deletion)"
  elif is_repo "$NESTED"; then
    if [ -n "$(git -C "$NESTED" status --porcelain 2>/dev/null)" ]; then
      keep_reason="it has uncommitted changes"
    elif git -C "$NESTED" rev-parse '@{u}' >/dev/null 2>&1 \
         && [ -n "$(git -C "$NESTED" log '@{u}..HEAD' --oneline 2>/dev/null)" ]; then
      keep_reason="it has commits that were never pushed"
    fi
  fi

  if [ -n "$keep_reason" ]; then
    warn "not deleting: $keep_reason"
    info "moving it to $BACKUP instead"
    run mkdir -p "$BACKUP"
    run mv "$NESTED" "$BACKUP/public_codes" || die "could not move the nested copy"
  else
    info "it is a clean, redundant second clone - nothing unsaved in it"
    if confirm "delete $NESTED ?"; then
      run rm -rf "$NESTED"
      info "deleted"
    else
      info "left alone (the outer clone is still used from here on)"
    fi
  fi
fi

# ----------------------------------------------------------------- 3. update
say "Updating the repo"
if [ -n "$(git -C "$REPO" status --porcelain 2>/dev/null)" ]; then
  warn "you have local changes here - skipping 'git pull' so nothing is clobbered:"
  git -C "$REPO" status --short | sed 's/^/      /'
else
  run git -C "$REPO" pull --ff-only || warn "pull failed - carrying on with what is on disk"
fi
[ -d "$PKG_DIR/circopt" ] || die "$PKG_DIR/circopt is missing - is the clone complete?"
info "at commit $(git -C "$REPO" log --oneline -1 2>/dev/null || echo '?')"

# ------------------------------------------------------------ 4. python env
say "Python environment"
PY=""
if command -v conda >/dev/null 2>&1; then
  CONDA_BASE="$(conda info --base 2>/dev/null)"
  # shellcheck disable=SC1091
  [ -f "$CONDA_BASE/etc/profile.d/conda.sh" ] && . "$CONDA_BASE/etc/profile.d/conda.sh"
  if conda env list | awk '{print $1}' | grep -qx "$ENV_NAME"; then
    info "conda env '$ENV_NAME' exists"
  else
    info "creating conda env '$ENV_NAME' (python 3.11)"
    run conda create -y -n "$ENV_NAME" python=3.11 >/dev/null || die "could not create the env"
  fi
  if [ "$DRY_RUN" = 0 ]; then
    conda activate "$ENV_NAME" || die "could not activate '$ENV_NAME'"
  fi
  PY="python"
else
  warn "conda not found - using the python already on PATH"
  PY="$(command -v python3 || command -v python)" || die "no python found"
fi

if [ "$DRY_RUN" = 0 ]; then
  if ! $PY -c "import RNA, flask, numpy, pandas, matplotlib, python_codon_tables" >/dev/null 2>&1; then
    info "installing dependencies"
    $PY -m pip install -q -r "$PKG_DIR/requirements.txt" || die "pip install failed"
  fi
  $PY -c "import RNA; print('    ViennaRNA', RNA.__version__)" || die "ViennaRNA will not import"
  if command -v derna >/dev/null 2>&1; then
    info "DERNA found (optional seed/comparison arm)"
  else
    info "DERNA not installed - optional. To add it:"
    info "    conda install -y -c conda-forge -c bioconda derna"
  fi
fi

# -------------------------------------------------------------- 5. self-test
say "Self-test"
if [ "$DRY_RUN" = 0 ]; then
  ( cd "$PKG_DIR" && $PY -m circopt.cli selftest ) || die "self-test failed"
else
  info "[dry-run] would run: python -m circopt.cli selftest"
fi

# ----------------------------------------------------------------- 6. serve
if [ "$DRY_RUN" = 0 ]; then
  mkdir -p "$RUNS_DIR" 2>/dev/null || true
fi
if [ "$START_SERVER" = 0 ] || [ "$DRY_RUN" = 1 ]; then
  say "Ready"
  info "start the server with:"
  info "    cd $PKG_DIR && conda activate $ENV_NAME"
  info "    CIRCOPT_RUNS=$RUNS_DIR python -m circopt.app --port $PORT"
  exit 0
fi

if command -v ss >/dev/null 2>&1 && ss -ltn 2>/dev/null | grep -q ":$PORT "; then
  die "port $PORT is already in use - stop the old server (Ctrl+C) or pass --port 5001"
fi

say "Starting the server on http://localhost:$PORT"
info "results are written to $RUNS_DIR (outside the repo, so 'git pull' never touches them)"
info "press Ctrl+C to stop"
cd "$PKG_DIR" || die "cannot cd to $PKG_DIR"
exec env CIRCOPT_RUNS="$RUNS_DIR" $PY -m circopt.app --port "$PORT"
