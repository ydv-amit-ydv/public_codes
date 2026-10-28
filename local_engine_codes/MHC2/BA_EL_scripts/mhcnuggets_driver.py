"""Runs inside the MHCnuggets environment (TensorFlow): python mhcnuggets_driver.py jobs.json out.csv

jobs.json = {"class": "II", "jobs": {allele: [peptides]}}. out.csv gets peptide, allele, ic50 (nM) and the rank of that IC50 among human-proteome
peptides (MHCnuggets' own percentile, a fraction 0-1). Kept separate so the Flask server never imports TensorFlow.
"""
import json
import os
import pickle
import sys
import tempfile

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "3")

import pandas as pd                                              # noqa: E402
import mhcnuggets.src.predict as P                               # noqa: E402

try:                                                             # MHCnuggets 2.4 calls Adam(lr=...), which Keras 3 no longer accepts
    import keras
    P.Adam = lambda lr=None, **k: keras.optimizers.Adam(learning_rate=lr or 0.001)
except ImportError:
    pass


def main():
    spec = json.load(open(sys.argv[1]))
    frames, unknown = [], []
    # MHCnuggets silently substitutes the "closest" allele when a name is not one it was trained for (and wrongly picks DRB1*01:01 for names
    # it cannot parse), so only alleles it really has are run and the rest are reported
    known = set(pickle.load(open(os.path.join(P.MHCNUGGETS_HOME, "data/production/examples_per_allele.pkl"), "rb")))
    with tempfile.TemporaryDirectory() as wd:
        for allele, peptides in spec["jobs"].items():
            if allele not in known:
                unknown.append(allele)
                continue
            inp, out = os.path.join(wd, "peps.txt"), os.path.join(wd, "out")      # no extension: the rank file is then out_ranks
            with open(inp, "w") as f:
                f.write("\n".join(peptides) + "\n")
            P.predict(class_=spec["class"], peptides_path=inp, mhc=allele, output=out, rank_output=True)
            ranks = pd.read_csv(out + "_ranks")
            ranks["allele"] = allele
            frames.append(ranks)
    json.dump(unknown, open(sys.argv[2] + ".unknown", "w"))
    out = pd.concat(frames, ignore_index=True).rename(columns={"human_proteome_rank": "rank"}) if frames else pd.DataFrame(columns=["peptide", "ic50", "rank", "allele"])
    out.to_csv(sys.argv[2], index=False)


if __name__ == "__main__":
    main()
