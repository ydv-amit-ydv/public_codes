"""TLimmuno2 prediction driver. Runs INSIDE TLimmuno2's own Python environment (TensorFlow 2.15 / Keras 2), started by
tlimmuno2_runner.py. Same encoding and models as Python/TLimmuno2.py of the repository (BLOSUM62 peptide of up to 21 aa,
34-residue MHC pseudosequence, binding-affinity network feeding the immunogenicity network), but scores a pair list in
batches and skips the background-rank step: the output is the raw immunogenicity score (0-1, higher = more immunogenic).

  python tlimmuno2_driver.py --home ~/TLimmuno2 --pairs in.csv --out out.csv [--threads 8]

in.csv: Peptide,Allele (allele = a name from Python/data/pseudosequence.2016.all.X.dat, e.g. DRB1_0101)
out.csv: Peptide,Allele,TLimmuno2   (empty when the peptide is longer than 21 aa or the allele is unknown)
"""
import argparse, csv, os, sys

ALPHA = list("ARNDCQEGHILKMFPSTWYVX")
PEP_LEN, MHC_LEN = 21, 34


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--home", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--batch", type=int, default=4096)
    a = ap.parse_args()
    home = os.path.abspath(os.path.expanduser(a.home))
    os.chdir(home)                                   # the repository loads its models by relative path
    os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")      # CPU only: these networks are small
    import numpy as np
    import pandas as pd
    import tensorflow as tf
    if a.threads:
        tf.config.threading.set_intra_op_parallelism_threads(a.threads)
        tf.config.threading.set_inter_op_parallelism_threads(max(1, a.threads // 4))

    blosum = pd.read_csv("Python/data/BLOSUM62.csv", comment="#")[ALPHA].loc[ALPHA]
    table = {c: blosum[c].to_numpy(float) for c in ALPHA}     # column c = BLOSUM62 scores of residue c against the alphabet
    pseudo = dict(l.split()[:2] for l in open("Python/data/pseudosequence.2016.all.X.dat") if l.strip())

    def enc(seq, L):
        seq = seq + "X" * (L - len(seq))
        return np.stack([table[c if c in table else "X"] for c in seq])

    with open(a.pairs, newline="") as fh:
        rows = list(csv.DictReader(fh))
    ok = [i for i, r in enumerate(rows) if r["Allele"] in pseudo and 5 <= len(r["Peptide"]) <= PEP_LEN]
    BA = tf.keras.models.load_model("Python/model/BAmodel", compile=False)
    BAm = tf.keras.models.Model(inputs=BA.input, outputs=BA.layers[-2].output)
    IMM = tf.keras.models.load_model("Python/model/TLimmuno2", compile=False)
    res = {}
    for s in range(0, len(ok), a.batch):
        part = ok[s:s + a.batch]
        P = np.stack([enc(rows[i]["Peptide"], PEP_LEN) for i in part])
        M = np.stack([enc(pseudo[rows[i]["Allele"]], MHC_LEN) for i in part])
        ba = BAm.predict([P, M], verbose=0, batch_size=1024)
        out = IMM.predict([P, M, ba], verbose=0, batch_size=1024).reshape(-1)
        res.update({i: float(v) for i, v in zip(part, out)})
        print(f"{min(s + a.batch, len(ok))}/{len(ok)}", flush=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Peptide", "Allele", "TLimmuno2"])
        for i, r in enumerate(rows):
            w.writerow([r["Peptide"], r["Allele"], f"{res[i]:.6f}" if i in res else ""])


if __name__ == "__main__":
    main()
