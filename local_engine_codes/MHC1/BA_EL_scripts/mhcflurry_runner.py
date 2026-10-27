"""MHCflurry 2.x (O'Donnell et al.; open source): binding affinity + antigen processing + presentation.

Output per peptide x allele (all from one model set, `models_class1_presentation`):
  MHCF_Affinity_nM / MHCF_Affinity_%Rank   binding affinity and its percentile (lower = better)
  MHCF_Processing                          antigen-processing score, 0-1 (higher = better)
  MHCF_Presentation / MHCF_Presentation_%Rank   combined presentation score and percentile
When the peptides were cut from a protein, 15 residues of flanking sequence on each side are passed
along, which is what the processing predictor is designed for.

MHCflurry needs torch (and Python 3.10+), so it runs in its own environment like BigMHC; the engine
starts its `mhcflurry-predict` command as a subprocess.
"""
import os, sys, tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, gfeller                       # noqa: E402
from common.jobs import run_cmd                           # noqa: E402

CHUNK = 400_000
FLANK = 15
ToolError = gfeller.ToolError
COLS = {"mhcflurry_affinity": "MHCF_Affinity_nM", "mhcflurry_affinity_percentile": "MHCF_Affinity_%Rank",
        "mhcflurry_processing_score": "MHCF_Processing", "mhcflurry_presentation_score": "MHCF_Presentation",
        "mhcflurry_presentation_percentile": "MHCF_Presentation_%Rank"}


def settings():
    m = config.model("mhcflurry")
    py = m.get("python")
    if not py:
        raise ToolError("No Python with mhcflurry found. Create an environment for it (see README) "
                        "or set models.mhcflurry.python in Settings.")
    cmd = m.get("command")
    if not cmd or cmd == "auto":
        cmd = os.path.join(os.path.dirname(py), "mhcflurry-predict")
    if not os.path.isfile(cmd):
        raise ToolError(f"mhcflurry-predict not found: {cmd} (pip install mhcflurry in that environment)")
    return {"python": py, "command": cmd}


def flanks(sequences, rows):
    """n_flank / c_flank for rows with Protein, Start, Peptide from the named sequences."""
    seq = dict(sequences)
    n, c = [], []
    for prot, start, pep in zip(rows["Protein"], rows["Start"], rows["Peptide"]):
        s = seq.get(prot, "")
        i = int(start) - 1
        n.append(s[max(0, i - FLANK):i] if s[i:i + len(pep)] == pep else "")
        c.append(s[i + len(pep):i + len(pep) + FLANK] if s[i:i + len(pep)] == pep else "")
    return n, c


def predict(rows, alleles, sequences=None, job=None, stage="mhcflurry"):
    """rows: DataFrame[Protein, Start, Peptide]; sequences: [(name, seq)] to supply flanks (optional).
    -> (DataFrame[Protein, Start, Peptide, Allele, MHCF_*], notes)."""
    s = settings()
    rows = rows[["Protein", "Start", "Peptide"]].drop_duplicates().reset_index(drop=True)
    if sequences:
        rows["n_flank"], rows["c_flank"] = flanks(sequences, rows)
    grid = rows.merge(pd.DataFrame({"allele": list(alleles)}), how="cross")
    parts = list(gfeller.chunks(list(range(len(grid))), CHUNK))
    frames = []
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV")}
    env["PATH"] = os.path.dirname(s["python"]) + ":/usr/bin:/bin"
    if job:
        job.update(stage, done=0, total=len(parts), msg=f"{len(grid):,} peptide-allele pairs")
    with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
        for i, idx in enumerate(parts, 1):
            part = grid.iloc[idx]
            inp, out = os.path.join(wd, f"in{i}.csv"), os.path.join(wd, f"out{i}.csv")
            cols = ["allele", "Peptide"] + (["n_flank", "c_flank"] if sequences else [])
            part[cols].rename(columns={"Peptide": "peptide"}).to_csv(inp, index=False)
            r = run_cmd([s["command"], inp, "--out", out, "--no-throw"], env=env, job=job)
            gfeller.require_output(out, r, "MHCflurry")
            t = pd.read_csv(out)
            res = part[["Protein", "Start", "Peptide", "allele"]].reset_index(drop=True)
            for src, dst in COLS.items():
                res[dst] = t[src].values
            frames.append(res.rename(columns={"allele": "Allele"}))
            if job:
                job.update(stage, done=i, msg=f"chunk {i}/{len(parts)}")
    out = pd.concat(frames, ignore_index=True)
    notes = []
    dead = [a for a, g in out.groupby("Allele") if g["MHCF_Presentation"].isna().all()]
    if dead:
        notes.append("MHCflurry does not support: " + ", ".join(dead))
        out = out[~out["Allele"].isin(dead)].reset_index(drop=True)
    return out, notes
