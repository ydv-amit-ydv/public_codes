"""MixMHCpred 3.0 (Gfeller lab): pan-allele MHC-I presentation predictor.

Output per peptide x allele: MixMHC_Score (higher = better) and MixMHC_%Rank (lower = better).
Runs `bash MixMHCpred` from the clone found by common/discover.py (or models.mixmhcpred in
Settings); it only needs numpy and pandas, so the server's own environment is enough.
"""
import os, sys, tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, gfeller                       # noqa: E402
from common.jobs import run_cmd                           # noqa: E402

CHUNK = 200_000
ToolError = gfeller.ToolError


def settings():
    m = config.model("mixmhcpred")
    if not m.get("home"):
        raise ToolError(f"MixMHCpred not found under {config.search_root()} "
                        "(clone https://github.com/GfellerLab/MixMHCpred there, or set models.mixmhcpred.home).")
    script = os.path.join(m["home"], "MixMHCpred")
    if not os.path.isfile(script):
        raise ToolError(f"MixMHCpred script not found: {script}")
    if not m.get("python"):
        raise ToolError("No Python with numpy and pandas found for MixMHCpred (set models.mixmhcpred.python).")
    return {"home": m["home"], "script": script, "python": m["python"]}


def supported(home):
    return gfeller.read_allele_list(os.path.join(home, "lib", "alleles_list.txt"), skip_header=True)


def predict(peptides, alleles, job=None, stage="mixmhcpred"):
    """-> (DataFrame[Peptide, Allele, MixMHC_Score, MixMHC_%Rank], notes)."""
    s = settings()
    names, notes = gfeller.split_alleles(alleles, supported(s["home"]), "MixMHCpred")
    if not names:
        return pd.DataFrame(columns=["Peptide", "Allele", "MixMHC_Score", "MixMHC_%Rank"]), notes
    peps = sorted(set(peptides))
    parts = list(gfeller.chunks(peps, CHUNK))
    frames, env = [], gfeller.tool_env(s["python"])
    if job:
        job.update(stage, done=0, total=len(parts), msg=f"{len(peps):,} peptides x {len(names)} alleles")
    with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
        for i, part in enumerate(parts, 1):
            inp, out = os.path.join(wd, f"in{i}.txt"), os.path.join(wd, f"out{i}.txt")
            with open(inp, "w") as f:
                f.write("\n".join(part) + "\n")
            r = run_cmd(["bash", s["script"], "-i", inp, "-o", out, "-a", ",".join(names)], env=env, job=job)
            gfeller.require_output(out, r, "MixMHCpred")
            t = gfeller.read_table(out)
            for g, canon in names.items():
                frames.append(pd.DataFrame({"Peptide": t["Peptide"], "Allele": canon,
                                            "MixMHC_Score": t[f"Score_{g}"], "MixMHC_%Rank": t[f"%Rank_{g}"]}))
            if job:
                job.update(stage, done=i, msg=f"chunk {i}/{len(parts)}")
    return pd.concat(frames, ignore_index=True), notes
