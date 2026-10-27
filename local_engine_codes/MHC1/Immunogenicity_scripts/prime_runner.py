"""PRIME 2.1 (Gfeller lab): MHC-I immunogenicity = binding (via MixMHCpred 3.0) + TCR recognition propensity.

Output per peptide x allele: PRIME_%Rank (lower = better; the value to rank peptides by),
PRIME_Score, and PRIME_%RankBinding (MixMHCpred's binding rank used inside PRIME).
Needs the PRIME clone and a MixMHCpred 3.0 clone; both are found beside this folder.
On first use the bundled C++ program is rebuilt for Linux, as the README instructs (needs g++).
"""
import os, shutil, subprocess, sys, tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, gfeller                       # noqa: E402
from common.jobs import run_cmd                           # noqa: E402

CHUNK = 200_000
ToolError = gfeller.ToolError


def settings():
    m = config.model("prime")
    if not m.get("home"):
        raise ToolError(f"PRIME not found under {config.search_root()} "
                        "(clone https://github.com/GfellerLab/PRIME there, or set models.prime.home).")
    mix = config.model("mixmhcpred")
    if not mix.get("home"):
        raise ToolError("PRIME needs MixMHCpred 3.0 (https://github.com/GfellerLab/MixMHCpred); "
                        "clone it beside PRIME or set models.mixmhcpred.home.")
    if not m.get("python"):
        raise ToolError("No Python with numpy and pandas found for PRIME (set models.prime.python).")
    return {"home": m["home"], "script": os.path.join(m["home"], "PRIME"), "python": m["python"],
            "mix": os.path.join(mix["home"], "MixMHCpred")}


def supported(home):
    return gfeller.read_allele_list(os.path.join(home, "lib", "alleles.txt"))


def ensure_built(home):
    """The repo ships a macOS/Windows binary; rebuild lib/PRIME.x here unless it is already a Linux ELF."""
    exe = os.path.join(home, "lib", "PRIME.x")
    try:
        with open(exe, "rb") as f:
            if f.read(4) == b"\x7fELF":
                return exe
    except OSError:
        pass
    if not shutil.which("g++"):
        raise ToolError("PRIME's C++ program has to be compiled for Linux and g++ is missing: sudo apt install g++")
    r = subprocess.run(["g++", "-O3", "PRIME.cc", "-o", "PRIME.x"], cwd=os.path.join(home, "lib"),
                       capture_output=True, text=True)
    if r.returncode:
        raise ToolError("Compiling PRIME.cc failed:\n" + r.stderr[-1500:])
    return exe


def predict(peptides, alleles, job=None, stage="prime"):
    """-> (DataFrame[Peptide, Allele, PRIME_%Rank, PRIME_Score, PRIME_%RankBinding], notes)."""
    s = settings()
    ensure_built(s["home"])
    os.makedirs(os.path.join(s["home"], "temp"), exist_ok=True)       # PRIME writes scratch files here
    names, notes = gfeller.split_alleles(alleles, supported(s["home"]), "PRIME")
    cols = ["Peptide", "Allele", "PRIME_%Rank", "PRIME_Score", "PRIME_%RankBinding"]
    if not names:
        return pd.DataFrame(columns=cols), notes
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
            r = run_cmd(["bash", s["script"], "-i", inp, "-o", out, "-a", ",".join(names), "-mix", s["mix"]],
                        env=env, job=job)
            gfeller.require_output(out, r, "PRIME")
            t = gfeller.read_table(out)
            for g, canon in names.items():
                frames.append(pd.DataFrame({"Peptide": t["Peptide"], "Allele": canon,
                                            "PRIME_%Rank": t[f"%Rank_{g}"], "PRIME_Score": t[f"Score_{g}"],
                                            "PRIME_%RankBinding": t[f"%RankBinding_{g}"]}))
            if job:
                job.update(stage, done=i, msg=f"chunk {i}/{len(parts)}")
    return pd.concat(frames, ignore_index=True), notes
