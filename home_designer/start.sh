#!/usr/bin/env bash
# Start Home Designer: http://127.0.0.1:8765  (designs are saved in ./designs)
cd "$(dirname "$0")"
PY=$(command -v python3 || command -v python)
if [ -z "$PY" ]; then echo "Python 3 is needed (https://www.python.org/downloads/)"; exit 1; fi
exec "$PY" server.py --open "$@"
