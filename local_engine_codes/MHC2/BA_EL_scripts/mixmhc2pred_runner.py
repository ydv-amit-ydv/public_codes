"""MixMHC2pred 2.1 (Gfeller lab): pan-allele MHC-II ligand predictor (HLA-DR, -DQ, -DP and mouse I-A/I-E).

Output per peptide x allele: MIX2_%Rank (lower = better). For peptides cut from a protein the 12-residue
context (3 residues either side of the peptide and 3 at each of its ends) is supplied, which is how the
model is meant to be used for epitope search; loose peptides are scored with --no_context.

The git repository has the Linux binary (`MixMHC2pred_unix`) but not the allele definitions: download
PWMdef_Human.zip (and PWMdef_Mouse.zip) from https://mixmhc2pred.unil.ch/PWMdef into <home>/PWMdef.
"""
import os, sys, tempfile
from functools import lru_cache

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import alleles as allele_lib, config, gfeller   # noqa: E402
from common.jobs import run_cmd                              # noqa: E402

CHUNK = 100_000
ToolError = gfeller.ToolError
HINT = ("Download PWMdef_Human.zip and PWMdef_Mouse.zip from https://mixmhc2pred.unil.ch/PWMdef (the 'PWM definition "
        "files' table) and unzip them into the PWMdef folder inside the MixMHC2pred folder.")


def settings():
    m = config.model("mixmhc2pred")
    if not m.get("home"):
        raise ToolError(f"MixMHC2pred not found under {config.search_root()} "
                        "(clone https://github.com/GfellerLab/MixMHC2pred there, or set models.mixmhc2pred.home).")
    exe = os.path.join(m["home"], "MixMHC2pred_unix")
    if not os.path.isfile(exe):
        raise ToolError(f"MixMHC2pred_unix not found in {m['home']}")
    pwm = os.path.join(m["home"], "PWMdef")
    if not os.path.isdir(pwm) or not any(f.endswith(".txt") for f in os.listdir(pwm)):
        raise ToolError("MixMHC2pred has no allele definitions yet (PWMdef is missing or empty). " + HINT)
    return {"home": m["home"], "exe": exe}


@lru_cache(maxsize=2)
def _alleles(home, stamp):
    """{OfficialName: AlleleName} from the Alleles_map_*.txt files shipped with the allele definitions."""
    out = {}
    pwm = os.path.join(home, "PWMdef")
    for fn in os.listdir(pwm):
        if fn.startswith("Alleles_map_") and fn.endswith(".txt"):
            with open(os.path.join(pwm, fn)) as f:
                for line in f:
                    if line.startswith("#") or line.startswith("OfficialName"):
                        continue
                    t = line.split("\t")
                    if len(t) >= 2:
                        out[t[0]] = t[1]
    return out


def allele_map(home):
    return _alleles(home, os.path.getmtime(os.path.join(home, "PWMdef")))


def resolve(alleles, home):
    """-> ({MixMHC2pred name: canonical allele}, notes)."""
    amap = allele_map(home)
    names, bad = {}, []
    for a in alleles:
        c = allele_lib.canonical_ii(a)
        n = amap.get(c)
        if n is None:
            g = allele_lib.mix2_name(c)
            n = g if g in set(amap.values()) else None
        if n:
            names[n] = c
        else:
            bad.append(a)
    return names, ([f"MixMHC2pred does not support: " + ", ".join(bad)] if bad else [])


def context(sequences, rows):
    """12-residue context per row: 3 upstream + first 3 + last 3 of the peptide + 3 downstream ('-' past the ends)."""
    seq = dict(sequences)
    out = []
    for prot, start, pep in zip(rows["Protein"], rows["Start"], rows["Peptide"]):
        s = seq.get(prot, "")
        i, j = int(start) - 1, int(start) - 1 + len(pep)
        if s[i:j] != pep:
            out.append("X" * 12)
            continue
        up = ("---" + s[max(0, i - 3):i])[-3:]
        down = (s[j:j + 3] + "---")[:3]
        out.append(up + pep[:3] + pep[-3:] + down)
    return out


def predict(rows, alleles, sequences=None, job=None, stage="mixmhc2pred"):
    """rows: DataFrame[Protein, Start, Peptide]. -> (DataFrame[Protein, Start, Peptide, Allele, MIX2_%Rank], notes)."""
    s = settings()
    names, notes = resolve(alleles, s["home"])
    cols = ["Protein", "Start", "Peptide", "Allele", "MIX2_%Rank"]
    if not names:
        return pd.DataFrame(columns=cols), notes
    rows = rows[["Protein", "Start", "Peptide"]].drop_duplicates().reset_index(drop=True)
    if sequences:
        rows["Context"] = context(sequences, rows)
    parts = list(gfeller.chunks(list(range(len(rows))), CHUNK))
    frames = []
    if job:
        job.update(stage, done=0, total=len(parts), msg=f"{len(rows):,} peptides x {len(names)} alleles")
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
        for i, idx in enumerate(parts, 1):
            part = rows.iloc[idx]
            inp, out = os.path.join(wd, f"in{i}.txt"), os.path.join(wd, f"out{i}.txt")
            with open(inp, "w") as f:
                for _, r in part.iterrows():
                    f.write(r["Peptide"] + ("\t" + r["Context"] if sequences else "") + "\n")
            cmd = [s["exe"], "-i", inp, "-o", out, "-a", *names]
            if not sequences:
                cmd.append("--no_context")
            r = run_cmd(cmd, env=env, job=job)
            gfeller.require_output(out, r, "MixMHC2pred")
            t = gfeller.read_table(out)
            if list(t["Peptide"]) != list(part["Peptide"]):
                raise ToolError("MixMHC2pred returned the peptides in a different order than they were sent")
            for n, canon in names.items():
                frames.append(pd.DataFrame({"Protein": part["Protein"].values, "Start": part["Start"].values,
                                            "Peptide": part["Peptide"].values, "Allele": canon,
                                            "MIX2_%Rank": t[f"%Rank_{n}"].values}))
            if job:
                job.update(stage, done=i, msg=f"chunk {i}/{len(parts)}")
    return pd.concat(frames, ignore_index=True), notes
