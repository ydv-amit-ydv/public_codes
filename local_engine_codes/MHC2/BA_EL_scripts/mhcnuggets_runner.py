"""MHCnuggets (Bhattacharya et al.): LSTM binding-affinity predictor for MHC class II, an independent third voter next to NetMHCIIpan and MixMHC2pred.

Output per peptide x allele: MHCnuggets_IC50_nM and MHCnuggets_%Rank (its percentile among human-proteome peptides, 0-100, lower = better).
It covers HLA-DR here (rare DR alleles use the closest trained one); DQ, DP and mouse alleles are skipped with a note. It needs its own
environment (TensorFlow): see the README. The server never imports it; mhcnuggets_driver.py runs inside that environment.
"""
import json
import os
import sys
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import alleles as allele_lib, config, envclean, gfeller, predcache   # noqa: E402
from common.jobs import run_cmd, THREADS                                          # noqa: E402

ToolError = gfeller.ToolError
DRIVER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mhcnuggets_driver.py")
HINT = ("conda create -n mhcnuggets python=3.10 -y && conda activate mhcnuggets && pip install mhcnuggets   "
        "(then set models.mhcnuggets.python in Settings if it is not found)")
COLS = ["Peptide", "Allele", "MHCnuggets_IC50_nM", "MHCnuggets_%Rank"]


def settings():
    m = config.model("mhcnuggets")
    if not m.get("python"):
        raise ToolError("MHCnuggets not found (no Python that can import mhcnuggets). Install it: " + HINT)
    return {"python": m["python"]}


def mn_name(a):
    """HLA-DRB1*15:01 -> HLA-DRB115:01; HLA-DQA1*05:01/DQB1*02:01 -> HLA-DQA105:01-DQB102:01 (MHCnuggets' own spelling)."""
    return a.replace("*", "").replace("/", "-")


def supported(alleles):
    ok = [a for a in alleles if a.startswith(("HLA-DRB", "HLA-DQA", "HLA-DPA"))]
    return ok, ([f"MHCnuggets skipped (human DR, DQ and DP only): {', '.join(a for a in alleles if a not in ok)}"] if len(ok) < len(alleles) else [])


def _run(s, jobs, job, unknown):
    """jobs: {allele: [peptides]} -> DataFrame[Peptide, Allele, ic50, rank], in a few processes at once. Alleles MHCnuggets was not trained for go to `unknown`."""
    perf = config.performance()
    k = max(1, min(perf["shards"], len(jobs)))
    names = list(jobs)
    groups = [names[i::k] for i in range(k)]
    share = max(1, perf["threads"] // k)
    env = envclean.clean(os.path.dirname(s["python"]))
    env["PATH"] = os.path.dirname(s["python"]) + ":" + env["PATH"]
    done, lock = [0], threading.Lock()
    if job:
        job.update("mhcnuggets", done=0, total=len(names), msg=f"{sum(len(v) for v in jobs.values()):,} peptide-allele pairs, {k} processes")

    def one(group):
        THREADS.set(share)
        with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
            spec, out = os.path.join(wd, "jobs.json"), os.path.join(wd, "out.csv")
            with open(spec, "w") as f:
                json.dump({"class": "II", "jobs": {mn_name(a): jobs[a] for a in group}}, f)
            n = sum(len(jobs[a]) for a in group)
            r = run_cmd([s["python"], DRIVER, spec, out], env=env, job=job,
                        task={"model": "MHCnuggets", "what": f"{len(group)} alleles ({group[0]}...)", "n": n, "sample": jobs[group[0]][:3]})
            if not os.path.exists(out):
                raise ToolError("MHCnuggets wrote no output:\n" + ((r.stderr or r.stdout) or "")[-1200:])
            df = pd.read_csv(out)
            with lock:
                unknown.extend(json.load(open(out + ".unknown")))
        with lock:
            done[0] += len(group)
            if job:
                job.update("mhcnuggets", done=done[0], msg=f"{done[0]}/{len(names)} alleles")
        return df
    with ThreadPoolExecutor(max_workers=k) as ex:
        frames = list(ex.map(one, groups))
    df = pd.concat(frames, ignore_index=True)
    back = {mn_name(a): a for a in jobs}
    df["allele"] = df["allele"].map(back)
    return df.rename(columns={"peptide": "Peptide", "allele": "Allele", "ic50": "MHCnuggets_IC50_nM"}).assign(
        **{"MHCnuggets_%Rank": lambda d: d["rank"] * 100}).drop(columns="rank")[COLS]


def predict(peptides, alleles, job=None, stage="mhcnuggets"):
    """-> (DataFrame[Peptide, Allele, MHCnuggets_IC50_nM, MHCnuggets_%Rank], notes). Cached per peptide and allele."""
    s = settings()
    ok, notes = supported([allele_lib.canonical_ii(a) for a in alleles])
    if not ok:
        return pd.DataFrame(columns=COLS), notes
    peps = sorted(set(peptides))

    unknown = []

    def compute(missing, als):
        return _run(s, {a: missing for a in als}, job, unknown)
    sig = ("mhcnuggets", s["python"])
    got = predcache.fetch("mhcnuggets", sig, peps, ok, compute, allele_lib.canonical_ii, job, stage)
    if unknown:
        known_names = {mn_name(a): a for a in ok}
        notes.append("MHCnuggets was not trained for: " + ", ".join(sorted({known_names.get(u, u) for u in unknown})))
    return (got if not got.empty else pd.DataFrame(columns=COLS)), notes
