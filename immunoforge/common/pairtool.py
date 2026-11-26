"""Shared plumbing for the class II immunogenicity tools that live in their own Python environment (ImmuScope: torch,
TLimmuno2: TensorFlow). Each is a driver script run by that environment's interpreter on a CSV of peptide-allele pairs;
this module maps the allele names, uses the prediction cache, runs the driver under the watchdog and reads the scores back.
The server never imports torch or TensorFlow.
"""
import os, re, tempfile

import pandas as pd

from . import alleles as allele_lib, config, predcache
from .jobs import run_cmd

MAX_LEN = 21               # both networks encode at most 21 residues
CHUNK = 400_000            # pairs per driver start (a start costs seconds: model loading)


class PairToolError(RuntimeError):
    pass


def settings(name, label, marker, hint):
    m = config.model(name)
    if not m.get("enabled", True):
        raise PairToolError(f"{label} is switched off in Settings.")
    home, py = m.get("home") or "", m.get("python") or ""
    if not home or not os.path.isfile(os.path.join(home, *marker)):
        raise PairToolError(f"{label} not found under {config.search_root()} ({hint}) or set models.{name}.home in Settings.")
    if not py or not os.path.isfile(py):
        raise PairToolError(f"No Python environment for {label}: {hint} or set models.{name}.python in Settings.")
    return {"home": home, "python": py}


def predict_pairs(tool, label, col, s, driver, driver_col, pairs, extra=(), job=None, stage=None, sig=()):
    """pairs: DataFrame[Peptide, Allele] (canonical class II names). -> (DataFrame[Peptide, Allele, col], notes)."""
    notes = []
    grid = pairs[["Peptide", "Allele"]].drop_duplicates().reset_index(drop=True)
    grid["pseudo"] = [allele_lib.pseudo_name(a) for a in grid["Allele"]]
    bad = sorted(set(grid.loc[grid["pseudo"].isna(), "Allele"]))
    if bad:
        notes.append(f"{label} has no pseudosequence for: " + ", ".join(bad[:8]) + (" ..." if len(bad) > 8 else "")
                     + " (DR chains and DQ/DP alpha-beta pairs only).")
    long_ = grid["Peptide"].str.len() > MAX_LEN
    if long_.any():
        notes.append(f"{label} reads at most {MAX_LEN} residues: {int(long_.sum()):,} longer pair(s) are not scored.")
    ok = grid[grid["pseudo"].notna() & ~long_].reset_index(drop=True)
    if ok.empty:
        return pd.DataFrame(columns=["Peptide", "Allele", col]), notes
    todo = ok.rename(columns={"Peptide": "peptide", "Allele": "allele"})
    parts = (len(todo) + CHUNK - 1) // CHUNK
    done = [0]
    if job and stage:
        job.update(stage, done=0, total=len(todo), msg=f"{len(todo):,} pairs")

    def compute(sub):
        out = []
        for i in range(0, len(sub), CHUNK):
            if job:
                job.check()
            part = sub.iloc[i:i + CHUNK]
            with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
                inp, res = os.path.join(wd, "in.csv"), os.path.join(wd, "out.csv")
                pd.DataFrame({"Peptide": part["peptide"].values, "Allele": part["pseudo"].values}).to_csv(inp, index=False)
                cmd = [s["python"], driver, "--home", s["home"], "--pairs", inp, "--out", res, *extra]
                pat = re.compile(r"(\d+)/(\d+)")

                def on_line(line, base=done[0]):
                    m = pat.search(line)
                    if m and job and stage and int(m.group(2)):
                        job.update(stage, done=base + round(len(part) * int(m.group(1)) / int(m.group(2))), msg=line.strip()[:80])
                r = run_cmd(cmd, env={k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}, job=job, merge=True, on_line=on_line,
                            task={"model": label, "what": f"{len(part):,} pairs", "n": len(part),
                                  "sample": [f"{a} {p}" for a, p in zip(part["allele"][:2], part["peptide"][:2])]})
                if r.returncode or not os.path.exists(res):
                    raise PairToolError(f"{label} failed (exit {r.returncode}):\n{(r.stdout or '')[-1500:]}")
                got = pd.read_csv(res)
                if len(got) != len(part) or not (got["Peptide"].astype(str).values == part["peptide"].values).all():
                    raise PairToolError(f"{label} returned rows that do not match the input; refusing to align by position")
                out.append(got[driver_col].astype(float).values)
            done[0] += len(part)
        import numpy as np
        return np.concatenate(out)

    ok["score"] = predcache.fetch_scores(tool, (tool,) + tuple(sig), todo, compute, allele_lib.canonical_ii, job)
    ok[col] = ok.pop("score")
    if job and stage:
        job.update(stage, done=len(todo), msg=f"{len(todo):,} pairs scored")
    res = ok[["Peptide", "Allele", col]]
    if res[col].isna().all():
        raise PairToolError(f"{label} scored nothing.")
    return res, notes
