"""Antigen-processing predictors from IEDB's next-generation T Cell Class I standalone:

  NetCTLpan 1.1   per 8-11-mer and allele: MHC binding, TAP transport, C-terminal cleavage, a combined
                  score (higher = better) and its %Rank (lower = better; IEDB suggests <= 1 for epitopes)
  NetChop 3.1     per residue: probability of proteasomal cleavage after it. For a peptide the value that
                  matters is the cleavage at its C-terminus (NetChop_Cterm, higher = better)

Both are DTU programs that IEDB redistributes inside its bundle (free for academic and non-profit use);
the engine runs the bundle's own `src/tcell_mhci.py` with its own Python environment. Both need the
protein the peptide came from (cleavage depends on the residues around it), so they apply when a protein
is cut into fragments, not to loose peptides.
"""
import json, os, re, shutil, sys, tempfile, threading
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, gfeller                       # noqa: E402
from common.jobs import run_cmd                           # noqa: E402

ToolError = gfeller.ToolError
AA = set("ACDEFGHIKLMNPQRSTVWY")
CTLPAN_LENGTHS = (8, 11)
_INCOMPATIBLE = re.compile(r"incompatible with netctlpan: ([^.]*)\.")


def settings():
    m = config.model("iedb_bundle")
    if not m.get("home"):
        raise ToolError(f"IEDB T Cell Class I bundle not found under {config.search_root()}. Download "
                        "IEDB_NG_TC1-*.tar.gz from https://downloads.iedb.org/nextgen-tools/tcell_mhci/LATEST/ "
                        "and extract it there (see README).")
    if not m.get("python"):
        raise ToolError(f"The bundle at {m['home']} has no .venv: run `uv sync` then `./configure` inside it.")
    for tool in ("tcsh", "gawk"):
        if not shutil.which(tool, path="/usr/bin:/bin"):
            raise ToolError(f"{tool} is missing (the bundle's programs need it): sudo apt install tcsh gawk")
    return {"home": m["home"], "python": m["python"]}


def _env(s):
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME", "VIRTUAL_ENV")}
    env["PATH"] = os.path.dirname(s["python"]) + ":/usr/bin:/bin"
    return env


def _fasta(sequences):
    """Numbered records; anything that is not one of the 20 amino acids becomes X (keeps positions,
    and a stray quote or '*' would otherwise make NetCTLpan return nothing for the whole protein)."""
    return "\n".join(f">{i}\n" + "".join(c if c in AA else "X" for c in seq) for i, (_, seq) in enumerate(sequences, 1))


def _run(s, spec, job, workdir, task=None):
    path = os.path.join(workdir, "job.json")
    with open(path, "w") as f:
        json.dump(spec, f)
    return run_cmd([s["python"], "src/tcell_mhci.py", "-j", path, "-f", "tsv"], env=_env(s), job=job, cwd=s["home"], task=task)


def _table(stdout, first_col):
    lines = stdout.splitlines()
    for i, l in enumerate(lines):
        if l.startswith(first_col):
            return pd.read_csv(pd.io.common.StringIO("\n".join(lines[i:])), sep="\t")
    return pd.DataFrame()


def netctlpan(sequences, lengths, alleles, job=None, stage="netctlpan"):
    """-> (DataFrame[Protein, Start, Peptide, Allele, CTLpan_*], notes).
    One bundle job per allele: with several alleles in one job the bundle labels every row with the whole
    allele list, so the alleles could only be told apart by guessing the row order."""
    s = settings()
    use = sorted(l for l in set(lengths) if CTLPAN_LENGTHS[0] <= l <= CTLPAN_LENGTHS[1])
    notes = []
    if set(lengths) - set(use):
        notes.append("NetCTLpan scores 8-11-mers only; other lengths get no NetCTLpan score")
    cols = ["Protein", "Start", "Peptide", "Allele", "CTLpan_MHC", "CTLpan_TAP", "CTLpan_Cleavage", "CTLpan_Combined", "CTLpan_%Rank"]
    alleles = list(alleles)
    if not use or not alleles:
        return pd.DataFrame(columns=cols), notes
    names = [n for n, _ in sequences]
    fasta = _fasta(sequences)
    workers = config.performance()["workers"]
    lock, done, unsupported, errors = threading.Lock(), [0], [], []
    if job:
        job.update(stage, done=0, total=len(alleles), msg=f"{len(alleles)} alleles on {len(sequences)} proteins")

    def one(allele):
        if job:
            job.check()
        spec = {"input_sequence_text": fasta, "peptide_length_range": [use[0], use[-1]],
                "alleles": allele.replace("H-2-", "H2-"),            # the bundle spells mouse alleles H2-Kb
                "predictors": [{"type": "processing", "method": "netctlpan", "tap_weight": 0.05,
                                "cleavage_weight": 0.225, "epitope_threshold": 1}]}
        with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
            r = _run(s, spec, job, wd, task={"model": "NetCTLpan", "what": f"{allele} \u00b7 {len(sequences)} protein(s)", "n": None})
        with lock:
            done[0] += 1
            if job:
                job.update(stage, done=done[0], msg=f"{done[0]}/{len(alleles)} alleles")
        if r.returncode:
            with lock:
                (unsupported if _INCOMPATIBLE.search(r.stderr + r.stdout) else errors).append((allele, (r.stderr or r.stdout)[-400:]))
            return None
        t = _table(r.stdout, "sequence_number")
        if t.empty:
            return None
        return pd.DataFrame({"Protein": [names[i - 1] for i in t["sequence_number"]], "Start": t["start"],
                             "Peptide": t["peptide"], "Allele": allele, "CTLpan_MHC": t["mhc_prediction"],
                             "CTLpan_TAP": t["tap_prediction_score"], "CTLpan_Cleavage": t["cleavage_prediction_score"],
                             "CTLpan_Combined": t["combined_prediction_score"], "CTLpan_%Rank": t["percentile_rank"]})

    with ThreadPoolExecutor(max_workers=workers) as ex:
        frames = [f for f in ex.map(one, alleles) if f is not None]
    if errors:
        raise ToolError(f"NetCTLpan failed for {errors[0][0]}:\n{errors[0][1]}")
    if unsupported:
        notes.append("NetCTLpan does not support: " + ", ".join(a for a, _ in unsupported))
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=cols)), notes


def netchop_cterm(sequences, rows, job=None, stage="netchop"):
    """-> DataFrame[Protein, Start, Peptide, NetChop_Cterm]: cleavage probability after each peptide's last residue."""
    s = settings()
    names = [n for n, _ in sequences]
    with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
        spec = {"input_sequence_text": _fasta(sequences), "peptide_length_range": [9, 9], "alleles": "HLA-A*02:01",
                "predictors": [{"type": "processing", "method": "netchop", "threshold": 0.5, "network_method": "c_term"}]}
        if job:
            job.update(stage, done=0, total=1, msg=f"{len(sequences)} proteins")
        r = _run(s, spec, job, wd, task={"model": "NetChop", "what": f"{len(sequences)} protein(s)", "n": len(rows), "sample": list(rows["Peptide"][:3])})
    if r.returncode:
        raise ToolError(f"NetChop failed:\n{(r.stderr or r.stdout)[-1500:]}")
    t = _table(r.stdout, "sequence_number")
    if t.empty:
        return pd.DataFrame(columns=["Protein", "Start", "Peptide", "NetChop_Cterm"])
    score = {(names[a - 1], p): v for a, p, v in zip(t["sequence_number"], t["position"], t["prediction_score"])}
    rows = rows[["Protein", "Start", "Peptide"]].drop_duplicates().copy()
    rows["NetChop_Cterm"] = [score.get((p, int(st) + len(pep) - 1)) for p, st, pep in zip(rows["Protein"], rows["Start"], rows["Peptide"])]
    if job:
        job.update(stage, done=1, msg=f"{len(rows):,} scores")
    return rows
