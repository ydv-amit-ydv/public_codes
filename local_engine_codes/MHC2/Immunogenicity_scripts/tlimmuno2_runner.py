"""TLimmuno2 (Zhang et al. 2024): deep-learning CD4 T-cell immunogenicity of a peptide-class II pair. A binding-affinity
network feeds an immunogenicity network (TensorFlow 2.15 / Keras 2; Keras 3 cannot load these SavedModels). The score is the
raw network output, 0-1, HIGHER = more immunogenic; the repository's own percentile against 90,001 background peptides is not
used (it is slow and the engine ranks within the candidates anyway). Peptides up to 21 residues, DR plus DQ/DP pairs.

Install (once): clone the TLimmuno2 repository beside this engine (its models are in the repository), then
    python -m venv .venv && .venv/bin/pip install "tensorflow-cpu==2.15.1" "numpy<2" pandas
"""
import os, sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from common import config, pairtool   # noqa: E402

PairToolError = pairtool.PairToolError
COL = "TLimmuno2"
DRIVER = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tlimmuno2_driver.py")
HINT = "clone TLimmuno2 there and make a Python environment with tensorflow 2.15 (see README, 'TLimmuno2')"


def settings():
    return pairtool.settings("tlimmuno2", "TLimmuno2", ("Python", "TLimmuno2.py"), HINT)


def predict(pairs, job=None, stage="tlimmuno2"):
    s = settings()
    return pairtool.predict_pairs("tlimmuno2", "TLimmuno2", COL, s, DRIVER, COL, pairs, extra=("--threads", str(config.performance()["threads"])),
                                  job=job, stage=stage, sig=(s["home"],))
