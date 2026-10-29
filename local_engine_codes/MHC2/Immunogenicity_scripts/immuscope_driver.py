"""ImmuScope-IM prediction driver. Runs INSIDE ImmuScope's own Python environment (torch), started by immuscope_runner.py.

The ImmuScope repository only ships training and test scripts, so this builds the same tensors the repository's data
pipeline builds (ImmuScope/utils/data_utils.get_peptide_embedding and SinInstanceBag) for a list of peptide-allele pairs,
and averages the instance probability of the 10 ensemble networks exactly as main_immunogenicity_test.py does.

  python immuscope_driver.py --home ~/ImmuScope --pairs in.csv --out out.csv [--threads 8] [--batch 2048]

in.csv: Peptide,Allele (allele = a name from data/raw/pseudosequence.2023.dat, e.g. DRB1_0101 or HLA-DQA10501-DQB10201)
out.csv: Peptide,Allele,ImmuScope_IM   (empty when the peptide is longer than 21 aa or the allele is unknown)
"""
import argparse, csv, os, sys

PEPTIDE_LEN, PEPTIDE_PAD, MHC_LEN = 21, 3, 34
MODEL_CNF = dict(bag_size=10, emb_size=16, conv_size=[9], conv_num=[64], conv_off=[3], dropout=0.25, peptide_pad=3)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--home", required=True)
    ap.add_argument("--weights", default=None, help="folder with ImmuScope-IM-0..9.pt (default <home>/weights/IM)")
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--batch", type=int, default=2048)
    a = ap.parse_args()
    home = os.path.abspath(os.path.expanduser(a.home))
    sys.path.insert(0, home)
    import numpy as np
    import torch
    if a.threads:
        torch.set_num_threads(a.threads)
    from ImmuScope.models.ImmuScope import ImmuScope
    from ImmuScope.utils.data_utils import ACIDS, get_mhc_name_seq

    wdir = a.weights or os.path.join(home, "weights", "IM")
    files = [os.path.join(wdir, f"ImmuScope-IM-{i}.pt") for i in range(10)]
    files = [f for f in files if os.path.isfile(f)]
    if not files:
        sys.exit(f"no ImmuScope-IM-*.pt weights in {wdir}")
    mhc_seq = get_mhc_name_seq(os.path.join(home, "data", "raw", "pseudosequence.2023.dat"))

    with open(a.pairs, newline="") as fh:
        rows = list(csv.DictReader(fh))
    ok = [i for i, r in enumerate(rows) if r["Allele"] in mhc_seq and 5 <= len(r["Peptide"]) <= PEPTIDE_LEN]
    pep = np.zeros((len(ok), 1, PEPTIDE_LEN + 2 * PEPTIDE_PAD), dtype=np.int64)
    mhc_emb = {k: np.array([ACIDS.index(x if x in ACIDS else "-") for x in v], dtype=np.int64) for k, v in mhc_seq.items()}
    mhc = np.zeros((len(ok), 1, MHC_LEN), dtype=np.int64)
    for j, i in enumerate(ok):
        p = [ACIDS.index(x if x in ACIDS else "-") for x in rows[i]["Peptide"]]
        pep[j, 0, PEPTIDE_PAD:PEPTIDE_PAD + len(p)] = p
        mhc[j, 0] = mhc_emb[rows[i]["Allele"]]

    score = np.zeros(len(ok))
    n = len(files)
    for k, f in enumerate(files):
        net = ImmuScope(**MODEL_CNF)
        net.load_state_dict(torch.load(f, map_location="cpu")["model_state_dict"], strict=False)
        net.eval()
        with torch.no_grad():
            for s in range(0, len(ok), a.batch):
                x = (torch.from_numpy(pep[s:s + a.batch]), torch.from_numpy(mhc[s:s + a.batch]))
                _, inst, _, _, _ = net(x)
                score[s:s + a.batch] += inst.reshape(-1).numpy() / n
        print(f"model {k + 1}/{n} done", flush=True)

    res = {i: score[j] for j, i in enumerate(ok)}
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Peptide", "Allele", "ImmuScope_IM"])
        for i, r in enumerate(rows):
            w.writerow([r["Peptide"], r["Allele"], f"{res[i]:.6f}" if i in res else ""])


if __name__ == "__main__":
    main()
