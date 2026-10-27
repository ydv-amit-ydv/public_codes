#!/bin/bash
# Start the Epitope Filter Engine. Needs only python3 (with the venv module).
# Creates ./.venv on first run (flask + pandas); BigMHC and NetMHCpan are found
# automatically in folders beside this one, or set in the Settings tab.
set -e
cd "$(dirname "$0")"
if [ ! -x .venv/bin/python ]; then
  python3 -m venv .venv || { echo "python3 -m venv failed: sudo apt install python3-venv"; exit 1; }
  .venv/bin/pip install -q -r requirements.txt
fi
exec .venv/bin/python engine_server.py
