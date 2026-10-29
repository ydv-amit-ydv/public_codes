"""ImmuScope-IM (Chen et al. 2025, GPL-3.0): deep-learning CD4 T-cell immunogenicity of a peptide-class II pair, an ensemble
of 10 networks (instance-level probability, 0-1, HIGHER = more immunogenic). Scores DR, DQ and DP; the pseudosequence table
has ~5,600 alleles. Peptides up to 21 residues.

Install (once): clone https://github.com/<ImmuScope repository> beside this engine, download ImmuScope-weights.tar.gz and the
pseudosequence table from the authors' Zenodo record (doi 10.5281/zenodo.14810445) into it, then create the environment:
    ImmuScope/weights/IM/ImmuScope-IM-0..9.pt        data/raw/pseudosequence.2023.dat
    python -m venv .venv && .venv/bin/pip install torch numpy pandas h5py scikit-learn ruamel.yaml tqdm click
The driver (immuscope_driver.py) builds the tensors itself, because the repository only ships training and test scripts.
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, pairtool   # noqa: E402

PairToolError = pairtool.PairToolError
COL = "ImmuScope_IM"
DRIVER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "immuscope_driver.py")
HINT = ("clone ImmuScope there with weights/IM/ImmuScope-IM-*.pt and data/raw/pseudosequence.2023.dat (see README, 'ImmuScope'), "
        "and make a Python environment with torch")


def settings():
    s = pairtool.settings("immuscope", "ImmuScope", ("ImmuScope", "models", "ImmuScope.py"), HINT)
    missing = [p for p in (("weights", "IM", "ImmuScope-IM-0.pt"), ("data", "raw", "pseudosequence.2023.dat")) if not os.path.isfile(os.path.join(s["home"], *p))]
    if missing:
        raise PairToolError("ImmuScope is missing " + ", ".join("/".join(m) for m in missing) + f" in {s['home']} (download from Zenodo 14810445).")
    return s


def predict(pairs, job=None, stage="immuscope"):
    """pairs: DataFrame[Peptide, Allele] -> (DataFrame[Peptide, Allele, ImmuScope_IM], notes)."""
    s = settings()
    return pairtool.predict_pairs("immuscope", "ImmuScope", COL, s, DRIVER, COL, pairs, extra=("--threads", str(config.performance()["threads"])),
                                  job=job, stage=stage, sig=(s["home"],))
