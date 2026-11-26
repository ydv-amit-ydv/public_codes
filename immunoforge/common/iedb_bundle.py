"""Running IEDB's next-generation standalones (T Cell Class I / II) from the engine.

Both are `uv` projects with their own Python 3.11 environment and a `src/tcell_mhc*.py` entry point that takes
a JSON job. Things learned the hard way, all handled here:
  * standard input must be closed, or the programs wait on it forever (run_cmd does that);
  * the programs call tcsh/gawk (and, for CD4episcore, R and perl) from PATH;
  * every record is sent as ">1", ">2"... and results are mapped back by number, so odd FASTA names cannot matter;
  * non-amino-acid characters in a sequence are replaced by X (a stray quote or '*' otherwise makes
    NetCTLpan return nothing for the whole protein).
"""
import json, os, shutil

import pandas as pd

from . import config, envclean, gfeller
from .jobs import run_cmd

ToolError = gfeller.ToolError
AA = set("ACDEFGHIKLMNPQRSTVWY")


def settings(model_key, entry, label, hint, tools=("tcsh", "gawk")):
    m = config.model(model_key)
    home, py = m.get("home") or "", m.get("python") or ""
    if not os.path.isfile(os.path.join(home, "src", entry)):
        raise ToolError(f"{label} not found under {config.search_root()}. {hint}")
    if not py:
        raise ToolError(f"The bundle at {home} has no .venv: run `uv sync` inside it.")
    for t in tools:
        if not shutil.which(t, path="/usr/bin:/bin"):
            raise ToolError(f"{t} is missing (the bundle's programs need it): sudo apt install {t}")
    return {"home": home, "python": py, "entry": entry}


def env(s):
    return envclean.clean(os.path.dirname(s["python"]))


def fasta(records):
    """records: list of sequences -> numbered FASTA with X for anything that is not an amino acid."""
    return "\n".join(f">{i}\n" + "".join(c if c in AA else "X" for c in seq) for i, seq in enumerate(records, 1))


def run(s, spec, workdir, job=None, task=None):
    path = os.path.join(workdir, "job.json")
    with open(path, "w") as f:
        json.dump(spec, f)
    return run_cmd([s["python"], os.path.join("src", s["entry"]), "-j", path, "-f", "tsv"],
                   env=env(s), job=job, cwd=s["home"], task=task)


def explain(r):
    """Why a bundle run produced no table: the lines of its output that name a problem, and the usual causes."""
    text = (r.stdout or "") + "\n" + (r.stderr or "")
    keys = ("error", "fail", "not found", "cannot", "can't", "unable", "no such", "denied", "traceback", "exception", "killed", "rint method")
    lines = [l.strip() for l in text.splitlines() if any(k in l.lower() for k in keys) and "FutureWarning" not in l and "libtinfo" not in l][-6:]
    hint = ""
    if "rint method" in text or "no callable" in text:
        hint = (" The tool's helper programs (perl NN-align, R) gave no output, so IEDB's own code then failed on an empty table. "
                "Check that `perl -v` and `R --version` work in a plain terminal, and start the engine from a shell without conda active.")
    return f"{' | '.join(lines)[:600]}.{hint}"


def table(stdout, first_col):
    """First tab-separated table whose header starts with `first_col` (logs may surround it)."""
    lines = stdout.splitlines()
    for i, l in enumerate(lines):
        if l.startswith(first_col):
            rows = []
            for l2 in lines[i + 1:]:
                if not l2.strip() or "\t" not in l2:
                    break
                rows.append(l2.split("\t"))
            return pd.DataFrame(rows, columns=l.split("\t")[:len(rows[0])] if rows else l.split("\t"))
    return pd.DataFrame()
