"""CD4episcore (Dhanda et al. 2018), IEDB's CD4 T-cell immunogenicity predictor, from the T Cell Class II standalone.

It scores the peptide alone (no HLA typing): CD4_Episcore is a percentile where LOWER means more immunogenic
(IEDB's convention; the Epitope Workbench ranks it the same way). It runs NN-align (Perl) and R, so both must be
installed (sudo apt install r-base-core perl). The class II tool only takes a protein as input for this method, so
every peptide is sent as its own one-window sequence, which reproduces IEDB's website scores exactly.
"""
import os, shutil, sys, tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, gfeller, iedb_bundle as bundle   # noqa: E402

ToolError = gfeller.ToolError
CHUNK = 1500
HINT = ("Download IEDB_NG_TC2-*.tar.gz from https://downloads.iedb.org/nextgen-tools/tcell_mhcii/LATEST/ and "
        "extract it there, then run `uv sync` inside it (see README).")


def settings():
    s = bundle.settings("iedb_bundle_ii", "tcell_mhcii.py", "IEDB T Cell Class II bundle", HINT, tools=("tcsh", "gawk", "perl", "R"))
    return s


def predict(peptides, job=None, stage="cd4episcore"):
    """-> DataFrame[Peptide, CD4_Episcore]."""
    s = settings()
    by_len = defaultdict(list)
    for p in sorted(set(peptides)):
        by_len[len(p)].append(p)
    tasks = [(L, part) for L, ps in sorted(by_len.items()) for part in gfeller.chunks(ps, CHUNK)]
    done = [0]
    if job:
        job.update(stage, done=0, total=len(tasks), msg=f"{sum(len(p) for _, p in tasks):,} peptides in {len(tasks)} runs")

    def one(task):
        L, part = task
        if job:
            job.check()
        spec = {"input_sequence_text": bundle.fasta(part), "peptide_length_range": [L, L], "alleles": "DRB1*01:01",
                "predictors": [{"type": "immunogenicity", "method": "cd4episcore"}]}
        with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
            r = bundle.run(s, spec, wd, job, task={"model": "CD4episcore", "what": f"{L}-mers", "n": len(part), "sample": part[:3]})
        done[0] += 1
        if job:
            job.update(stage, done=done[0], msg=f"{done[0]}/{len(tasks)} runs finished")
        t = bundle.table(r.stdout, "seq_num")
        if t.empty:
            raise ToolError(f"CD4episcore returned nothing for {len(part)} {L}-mers:\n{(r.stderr or r.stdout)[-800:]}")
        return pd.DataFrame({"Peptide": [part[int(i) - 1] for i in t["seq_num"]], "CD4_Episcore": t["score"].astype(float).values})

    with ThreadPoolExecutor(max_workers=max(1, min(len(tasks), config.performance()["workers"]))) as ex:
        frames = list(ex.map(one, tasks))
    return pd.concat(frames, ignore_index=True).drop_duplicates("Peptide") if frames else pd.DataFrame(columns=["Peptide", "CD4_Episcore"])
