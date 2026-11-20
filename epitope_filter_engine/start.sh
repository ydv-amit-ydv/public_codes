#!/bin/bash
# Start the Epitope Filter Engine.
#
# Interpreter, first match wins:
#   1. $ENGINE_PYTHON                       explicit override
#   2. the active conda environment         if it has flask and pandas
#   3. the conda env "epitope-engine"       created by: conda env create -f environment.yml
#   4. ./.venv                              made here on first run (flask + pandas)
#
# The predictors are NOT run by this interpreter: each has its own environment and
# is found automatically beside this folder (or under $ENGINE_MODELS_DIR), or set in
# the Settings tab. See INSTALL.md.
set -e
cd "$(dirname "$0")"

has_deps() { [ -x "$1" ] && "$1" -c "import flask, pandas" >/dev/null 2>&1; }

PY=""

# 1. explicit override
if [ -n "${ENGINE_PYTHON:-}" ]; then
  PY="${ENGINE_PYTHON/#\~/$HOME}"
  has_deps "$PY" || { echo "ENGINE_PYTHON=$PY cannot import flask and pandas"; exit 1; }

# 2. the conda environment that is already active
elif [ -n "${CONDA_PREFIX:-}" ] && has_deps "$CONDA_PREFIX/bin/python"; then
  PY="$CONDA_PREFIX/bin/python"

# 3. a conda environment named epitope-engine, active or not
else
  CONDA_ROOT=""
  if [ -n "${CONDA_EXE:-}" ]; then CONDA_ROOT="$(dirname "$(dirname "$CONDA_EXE")")"; fi
  for base in "${CONDA_PREFIX_1:-}" "${CONDA_PREFIX:-}" "$CONDA_ROOT" \
              "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/miniforge3" "$HOME/mambaforge" "/opt/conda"; do
    [ -n "$base" ] || continue
    base="${base%%/envs/*}"
    if has_deps "$base/envs/epitope-engine/bin/python"; then
      PY="$base/envs/epitope-engine/bin/python"
      break
    fi
  done
fi

# 4. a local virtualenv, created on first run
if [ -z "$PY" ]; then
  if ! has_deps .venv/bin/python; then
    echo "Creating .venv (flask + pandas)."
    echo "Prefer conda? Ctrl-C, then: conda env create -f environment.yml"
    python3 -m venv .venv || {
      echo "python3 -m venv failed: sudo apt install -y python3-venv"
      echo "or use conda: conda env create -f environment.yml"
      exit 1
    }
    .venv/bin/pip install -q -r requirements.txt
  fi
  PY=.venv/bin/python
fi

echo "Engine interpreter: $PY"
exec "$PY" engine_server.py
