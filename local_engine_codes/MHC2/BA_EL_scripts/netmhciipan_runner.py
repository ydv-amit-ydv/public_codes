"""NetMHCIIpan (DTU) through IEDB's next-generation T Cell Class II standalone.

  EL: score (higher = better) and %Rank (lower = better, <= 2 is the usual binder cut-off)
  BA: predicted IC50 in nM and %Rank
Versions: "4.3" (newest), "4.2", or "default" = the unsuffixed IEDB method, which is what IEDB's own website
returns (it reproduces the scores in IEDB web exports). The bundle's version-suffixed names exist for 4.2 and 4.3 only.

Alleles are sent several to a job (the bundle labels rows correctly for class II); jobs run in parallel.
An allele the program does not know produces no rows and no error, so missing alleles are reported afterwards.
"""
import os, sys, tempfile, threading
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import alleles as allele_lib, config, gfeller, iedb_bundle as bundle   # noqa: E402

ToolError = gfeller.ToolError
VERSIONS = ("4.3", "4.2", "default")
HINT = ("Download IEDB_NG_TC2-*.tar.gz from https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/ and "
        "extract it there, then run `uv sync` inside it (see README).")
PEPTIDES_PER_JOB = 2000
COLS = ["Protein", "Start", "Peptide", "Allele", "NII_Core", "NII_Score_EL", "NII_%Rank_EL", "NII_IC50_nM", "NII_%Rank_BA"]


def settings():
    return bundle.settings("iedb_bundle_ii", "tcell_mhcii.py", "IEDB T Cell Class II bundle", HINT)


def method(kind, version):
    return f"netmhciipan_{kind}" + ("" if version == "default" else f"-{version}")


def predict(rows, alleles, lengths, sequences=None, kinds=("el", "ba"), version="4.3", job=None, stage="netmhciipan"):
    """rows: DataFrame[Protein, Start, Peptide] (the windows wanted). sequences: [(name, seq)] when peptides were cut
    from a protein; None for loose peptides. -> (DataFrame keyed by Protein/Start/Peptide/Allele (scan) or
    Peptide/Allele (loose peptides) with NII_* columns, notes)."""
    s = settings()
    if version not in VERSIONS:
        raise ValueError(f"NetMHCIIpan version must be one of {VERSIONS}")
    alleles = [allele_lib.canonical_ii(a) for a in alleles]
    cfg = config.model("iedb_bundle_ii")
    workers = config.performance()["workers"]
    per_call = max(1, min(config.model("netmhcpan").get("alleles_per_call") or 3, -(-len(alleles) // workers)) or 1)
    chunks = list(gfeller.chunks(alleles, per_call))
    lengths = sorted(set(int(x) for x in lengths))
    names = [n for n, _ in sequences] if sequences else None
    peptides = sorted(set(rows["Peptide"]))
    if sequences:
        tasks = [(k, c, None) for k in kinds for c in chunks]
    else:
        # The class II standalone (v0.2.3-beta) fails on `peptide_list` input, even in its own README example, so each
        # loose peptide is sent as its own one-window sequence, grouped by length.
        by_len = {}
        for p in peptides:
            by_len.setdefault(len(p), []).append(p)
        groups = [(L, part) for L, ps in sorted(by_len.items()) for part in gfeller.chunks(ps, PEPTIDES_PER_JOB)]
        tasks = [(k, c, g) for k in kinds for c in chunks for g in groups]
    lock, done, frames = threading.Lock(), [0], {"el": [], "ba": []}
    if job:
        job.update(stage, done=0, total=len(tasks), msg=f"{len(tasks)} runs, {workers} in parallel")

    def one(task):
        kind, als, group = task
        if job:
            job.check()
        spec = {"alleles": ",".join(als),
                "predictors": [{"type": "binding", "method": method(kind, version)}]}
        if sequences:
            spec.update({"input_sequence_text": bundle.fasta([sq for _, sq in sequences]),
                         "peptide_length_range": [lengths[0], lengths[-1]]})
        else:
            spec.update({"input_sequence_text": bundle.fasta(group[1]), "peptide_length_range": [group[0], group[0]]})
        with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
            r = bundle.run(s, spec, wd, job)
        with lock:
            done[0] += 1
            if job:
                job.update(stage, done=done[0], msg=f"{done[0]}/{len(tasks)} runs finished")
        if r.returncode:
            raise ToolError(f"NetMHCIIpan failed:\n{(r.stderr or r.stdout)[-1500:]}")
        t = bundle.table(r.stdout, "allele")
        if group is not None and not t.empty:                      # restore the real peptide for each one-window record
            t["peptide"] = [group[1][int(i) - 1] for i in t["seq_num"]]
        return kind, t

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for kind, t in ex.map(one, tasks):
            if not t.empty:
                frames[kind].append(t)

    def shape(parts, kind):
        if not parts:
            return pd.DataFrame()
        t = pd.concat(parts, ignore_index=True)
        out = pd.DataFrame({"Allele": t["allele"], "Peptide": t["peptide"], "NII_Core": t["core"]})
        if sequences:
            out.insert(0, "Start", t["start"].astype(int).values)
            out.insert(0, "Protein", [names[int(i) - 1] for i in t["seq_num"]])
        if kind == "el":
            out["NII_Score_EL"] = t["score"].astype(float).values
            out["NII_%Rank_EL"] = t["percentile"].astype(float).values
        else:
            out["NII_IC50_nM"] = t["ic50"].astype(float).values
            out["NII_%Rank_BA"] = t["percentile"].astype(float).values
        return out

    keys = (["Protein", "Start"] if sequences else []) + ["Peptide", "Allele"]
    el, ba = shape(frames["el"], "el"), shape(frames["ba"], "ba")
    if el.empty and ba.empty:
        return pd.DataFrame(columns=COLS), [f"NetMHCIIpan returned nothing for: {', '.join(alleles)}"]
    if el.empty or ba.empty:
        out = el if not el.empty else ba
    else:
        out = el.merge(ba.drop(columns=["NII_Core"]), on=keys, how="outer")
    got = set(out["Allele"])
    missing = [a for a in alleles if a not in got]
    notes = [f"NetMHCIIpan does not support: {', '.join(missing)}"] if missing else []
    # keep only the windows that were asked for (the bundle returns every length between the smallest and largest)
    want = rows[["Protein", "Start", "Peptide"]].drop_duplicates() if sequences else pd.DataFrame({"Peptide": peptides})
    return out.merge(want, on=[k for k in keys if k != "Allele"], how="inner"), notes
