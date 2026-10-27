"""Run BigMHC (EL and/or IM) on peptides or a scanned protein.

BigMHC (Albert et al., Nat Mach Intell 2023) is a pan-allele MHC-I model:
  EL = antigen presentation (eluted ligand), IM = immunogenicity.

BigMHC needs torch, so it lives in its own conda env. This module never imports
torch: it starts BigMHC's predict.py with that env's interpreter, taken from
engine_config.json (models.bigmhc.python) or --python. It can therefore be
imported by the Flask server (Python 3.12) and still drive a Python 3.9 env.

Example
  python3 bigmhc_runner.py --fasta orf1.fa --lengths 9 10 --alleles HLA-A02:01 HLA-B07:02 -o out.csv
"""
import argparse, os, re, sys, tempfile

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import alleles as allele_lib, config   # noqa: E402
from common.jobs import run_cmd                    # noqa: E402

AA = set("ACDEFGHIKLMNPQRSTVWY")
COLS = ["seq #", "Protein", "start", "end", "peptide length", "peptide"]


class BigMHCError(RuntimeError):
    pass


def parse_fasta_text(text):
    """FASTA or a bare sequence -> [(name, seq)]."""
    seqs, name, buf = [], None, []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if buf:
                seqs.append((name or f"seq{len(seqs)+1}", "".join(buf)))
            parts = line[1:].split()
            name, buf = (parts[0] if parts else None), []
        elif line:
            buf.append(line.upper())
    if buf:
        seqs.append((name or f"seq{len(seqs)+1}", "".join(buf)))
    return seqs


def windows(seqs, lengths):
    rows = []
    for i, (name, seq) in enumerate(seqs, 1):
        for L in lengths:
            for s in range(len(seq) - L + 1):
                pep = seq[s:s + L]
                if not set(pep) - AA:
                    rows.append((i, name, s + 1, s + L, L, pep))
    return pd.DataFrame(rows, columns=COLS)


def given_peptides(peps):
    rows = [(i, "peptide", 1, len(p), len(p), p.upper()) for i, p in enumerate(peps, 1)]
    return pd.DataFrame(rows, columns=COLS)


def settings(python=None, home=None, device=None, jobs=None):
    """Explicit arguments win over engine_config.json (which may say 'auto' = discover)."""
    m = config.model("bigmhc")
    s = {"python": python or m.get("python"),
         "home": os.path.expanduser(home) if home else m.get("home"),
         "device": device or m.get("device") or "cpu",
         "jobs": m.get("jobs", 0) if jobs is None else jobs}
    root = config.search_root()
    if not s["python"]:
        raise BigMHCError("No Python with torch/numpy/pandas/psutil found. Create the env "
                          "(conda create -n bigmhc python=3.9; pip install torch numpy pandas psutil) "
                          "or set models.bigmhc.python in Settings.")
    if not s["home"]:
        raise BigMHCError(f"BigMHC clone not found under {root} (need a folder with src/predict.py "
                          "and models/). Clone it there or set models.bigmhc.home in Settings.")
    if not os.path.isfile(s["python"]):
        raise BigMHCError(f"BigMHC python not found: {s['python']} (set models.bigmhc.python)")
    if not os.path.isfile(os.path.join(s["home"], "src", "predict.py")):
        raise BigMHCError(f"BigMHC not found at {s['home']} (set models.bigmhc.home)")
    return s


def _clean_env():
    return {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}


_LOADED = re.compile(r"loading model .*\.\.\.done")
_BATCH = re.compile(r"batch (\d+)/(\d+)")


def run_bigmhc(s, model, pairs, workdir, job=None, stage=None):
    """pairs: DataFrame with 'allele' and 'peptide'. Returns scores aligned to pairs.
    With a job, BigMHC's own progress lines ("loading model ...done", "batch i/N") drive the bar."""
    inp = os.path.join(workdir, f"in_{model}.csv")
    pairs[["allele", "peptide"]].to_csv(inp, header=False, index=False)    # col0 allele, col1 peptide
    out = inp + ".prd"
    cmd = [s["python"], os.path.join(s["home"], "src", "predict.py"),
           f"-i={inp}", f"-m={model}", "-a=0", "-p=1", "-c=0",
           f"-d={s['device']}", f"-j={s['jobs']}", f"-o={out}"]
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "PYTHONHOME")}
    try:
        n_models = len([d for d in os.listdir(os.path.join(s["home"], "models")) if d.startswith("bat")]) or 7
    except OSError:
        n_models = 7
    state = {"loaded": 0, "total": None}

    def on_line(line):
        # one 0-100 scale per stage, so the bar never jumps back: loading the networks is the
        # first 40%, scoring the batches the rest (BigMHC only reveals the batch count late)
        if job is None or stage is None:
            return
        if _LOADED.search(line):
            state["loaded"] += 1
            job.update(stage, done=round(40 * state["loaded"] / n_models), total=100,
                       msg=f"loading network {state['loaded']}/{n_models}")
        m = _BATCH.search(line)
        if m:
            i, n = int(m.group(1)), int(m.group(2))
            job.update(stage, done=40 + round(60 * (i - 1) / n), total=100, msg=f"scoring batch {i}/{n}")

    if job and stage:
        job.update(stage, done=0, total=100, msg=f"starting BigMHC {model.upper()} on {len(pairs):,} pairs")
    r = run_cmd(cmd, env=env, job=job, merge=True, on_line=on_line)
    if r.returncode or not os.path.exists(out):
        raise BigMHCError(f"BigMHC {model} failed:\n{' '.join(cmd)}\n{r.stdout[-3000:]}")
    prd = pd.read_csv(out)
    score_cols = [c for c in prd.columns if c.lower().startswith("bigmhc")]
    if not score_cols:
        raise BigMHCError(f"No BigMHC score column in {out}; columns are {list(prd.columns)}")
    if len(prd) != len(pairs):
        raise BigMHCError(f"BigMHC returned {len(prd)} rows for {len(pairs)} inputs")
    pep_col = next((c for c in prd.columns if c.lower() in ("pep", "peptide")), None)
    if pep_col is not None and not (prd[pep_col].astype(str).values == pairs["peptide"].values).all():
        raise BigMHCError("BigMHC output order differs from input; refusing to align by position")
    if job and stage:
        job.finish_stage(stage, f"{len(pairs):,} pairs scored")
    return prd[score_cols[0]].values


def predict_pairs(pairs, models=("el", "im"), job=None, **overrides):
    """pairs: DataFrame[peptide, allele] (any combination). Adds BigMHC_EL / BigMHC_IM columns."""
    s = settings(**overrides)
    grid = pairs.copy()
    grid["allele"] = [allele_lib.canonical(a) for a in grid["allele"]]
    with tempfile.TemporaryDirectory(dir=config.tmp_dir()) as wd:
        for m in models:
            grid[f"BigMHC_{m.upper()}"] = run_bigmhc(s, m, grid, wd, job=job, stage=f"bigmhc_{m}")
    return grid


def predict(peps, alleles, models=("el", "im"), job=None, **overrides):
    """peps: DataFrame from windows()/given_peptides(). Returns peptide x allele rows
    with BigMHC_EL / BigMHC_IM columns."""
    peps = peps.copy()
    peps["peptide index"] = peps["peptide"].factorize()[0] + 1
    grid = peps.merge(pd.DataFrame({"allele": [allele_lib.canonical(a) for a in alleles]}), how="cross")
    return predict_pairs(grid, models, job=job, **overrides)


def workbench_im(grid):
    """IM scores in the column layout of IEDB/HEV/Raw/MHC_1_Immunogenicity_*.csv."""
    wb = grid.rename(columns={"BigMHC_IM": "immunogenicity score"})
    return wb[["seq #", "peptide", "start", "end", "peptide length", "allele",
               "peptide index", "immunogenicity score"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--fasta")
    src.add_argument("--peptides", nargs="+")
    ap.add_argument("--lengths", nargs="+", type=int, default=[9], help="k-mer sizes for --fasta")
    ap.add_argument("--alleles", nargs="+", required=True)
    ap.add_argument("--models", nargs="+", default=["el", "im"], choices=["el", "im"])
    ap.add_argument("--python", help="BigMHC env interpreter (default: engine_config.json)")
    ap.add_argument("--bigmhc-home", help="BigMHC clone (default: engine_config.json)")
    ap.add_argument("--device", help="cpu, all, or GPU indices e.g. 0,1")
    ap.add_argument("-o", "--out", default="bigmhc_results.csv")
    ap.add_argument("--workbench-csv", help="also write IM scores in the Epitope Workbench schema "
                    "(BigMHC IM and IEDB's score differ in scale; keep them in separate folders)")
    a = ap.parse_args()

    if a.fasta:
        peps = windows(parse_fasta_text(open(a.fasta).read()), a.lengths)
    else:
        peps = given_peptides(a.peptides)
    if peps.empty:
        sys.exit("No valid peptides.")
    try:
        grid = predict(peps, a.alleles, a.models, python=a.python, home=a.bigmhc_home, device=a.device)
    except (BigMHCError, config.ConfigError) as e:
        sys.exit(str(e))
    grid.to_csv(a.out, index=False)
    print(f"wrote {a.out} ({len(grid)} rows)", file=sys.stderr)
    if a.workbench_csv:
        if "BigMHC_IM" not in grid:
            sys.exit("--workbench-csv needs the im model")
        workbench_im(grid).to_csv(a.workbench_csv, index=False)
        print(f"wrote {a.workbench_csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
