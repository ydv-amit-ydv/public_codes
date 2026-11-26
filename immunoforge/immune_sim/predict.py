#!/usr/bin/env python3
"""
Antigen -> predicted epitopes -> C-ImmSim parameters.

  python3 predict.py P0DTC2                    # UniProt accession (SARS-CoV-2 spike)
  python3 predict.py myprotein.fasta           # local FASTA
  python3 predict.py P0DTC2 --alleles HLA-A*02:01,HLA-DRB1*01:01

Steps
  1. fetch/read the antigen sequence
  2. tile into 9mers (class I) and 15mers (class II)
  3. predict IC50 for every peptide x allele with the trained models
  4. rank epitopes and write results/<name>_epitopes.csv
  5. calibrate C-ImmSim bit-string parameters from the predictions and
     write results/<name>_cimmsim.txt
"""
import argparse
import pickle
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

from model import AAS, CANON_LEN, encode, inv_log50k

HERE = Path(__file__).parent
MODELS = HERE / "models"
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

NBITSTR = 12            # C-ImmSim bit-string length (2^12 = 4096 receptors)
MINMATCH = 9            # C-ImmSim default binding threshold
STRONG_NM = 50.0        # strong binder
WEAK_NM = 500.0         # weak binder / standard cut-off


def fetch_uniprot(acc):
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    with urllib.request.urlopen(url, timeout=60) as r:
        txt = r.read().decode()
    return parse_fasta(txt)


def parse_fasta(txt):
    header, seq = "", []
    for line in txt.splitlines():
        line = line.strip()
        if line.startswith(">"):
            if seq:
                break
            header = line[1:]
        elif line:
            seq.append(line)
    s = re.sub(f"[^{AAS}]", "", "".join(seq).upper())
    if not s:
        raise ValueError("no amino acid sequence found")
    return header, s


def load_models(wanted=None):
    out = {}
    for p in sorted(MODELS.glob("HLA*.pkl")):
        with open(p, "rb") as f:
            d = pickle.load(f)
        if wanted and d["allele"] not in wanted:
            continue
        out[d["allele"]] = d
    if not out:
        raise ValueError("no models found - run model.py first")
    return out


def tile(seq, k):
    return [(i + 1, seq[i:i + k]) for i in range(len(seq) - k + 1)]


def scan(seq, models):
    rows = []
    for allele, d in models.items():
        L = d["length"]
        peps = tile(seq, L)
        if not peps:
            continue
        pos = [p[0] for p in peps]
        s = [p[1] for p in peps]
        score = d["model"].predict(encode(s, L))
        nm = inv_log50k(score)
        rows.append(pd.DataFrame({
            "allele": allele, "mhc_class": d["mhc_class"], "start": pos,
            "peptide": s, "length": L, "score": score, "pred_IC50_nM": nm,
        }))
    df = pd.concat(rows, ignore_index=True)
    df["percentile_rank"] = (df.groupby("allele").score
                               .rank(pct=True, ascending=False) * 100).round(2)
    df["call"] = np.where(df.pred_IC50_nM < STRONG_NM, "strong",
                 np.where(df.pred_IC50_nM < WEAK_NM, "weak", "non-binder"))
    return df.sort_values("pred_IC50_nM").reset_index(drop=True)


def score_peptides(peptides, allele, models):
    """Predicted IC50 (nM) for given peptides against one allele.

    scan() tiles a protein; this scores peptides you already have, which is
    what a substitution test needs - the same model, the same encoding, one
    residue different.
    """
    d = models.get(allele)
    if d is None:
        raise ValueError(f"no model for {allele}")
    L = d["length"]
    peps = [str(p) for p in peptides]
    bad = [p for p in peps if len(p) != L]
    if bad:
        raise ValueError(f"{allele} scores {L}-mers; got {len(bad[0])}")
    return inv_log50k(d["model"].predict(encode(peps, L)))


# --- mapping predictions onto C-ImmSim's bit-string world -------------------

def affinity_to_match(nm):
    """
    C-ImmSim binds when two bit strings share >= minmatch complementary bits,
    and affinity rises with the match count m (minmatch..NBITSTR).
    Map predicted IC50 onto that integer ladder on a log scale:
        <=50 nM   -> 12 (perfect complement)
        ~500 nM   -> 10
        >=50000nM -> below threshold, no binding
    """
    nm = float(np.clip(nm, 1.0, 5e4))
    frac = 1.0 - np.log(nm / 1.0) / np.log(5e4 / 1.0)      # 1.0 strong .. 0.0 weak
    m = MINMATCH - 1 + frac * (NBITSTR - MINMATCH + 1)
    return int(np.clip(round(m), 0, NBITSTR))


def make_bitstring(mhc_bits, n_match, rng):
    """
    Build a peptide bit string with exactly n_match complementary bits
    against the given MHC string.
    """
    comp = [1 - b for b in mhc_bits]              # perfect complement
    out = comp[:]
    flip = rng.choice(NBITSTR, size=NBITSTR - n_match, replace=False)
    for i in flip:
        out[i] = 1 - out[i]                       # break that bit's match
    return out


def bits_to_int(bits):
    v = 0
    for b in bits:
        v = (v << 1) | b
    return v


def cimmsim_params(df, seq_name, n_epitopes=4, seed=1111):
    rng = np.random.default_rng(seed)
    mhc = {}
    lines = []

    best = (df[df.call != "non-binder"]
              .sort_values("pred_IC50_nM")
              .groupby("allele").head(n_epitopes))
    if best.empty:
        best = df.groupby("allele").head(1)

    for allele in sorted(best.allele.unique()):
        mhc[allele] = list(rng.integers(0, 2, NBITSTR))

    rows = []
    for _, r in best.iterrows():
        m = affinity_to_match(r.pred_IC50_nM)
        bits = make_bitstring(mhc[r.allele], m, rng)
        rows.append({
            "allele": r.allele, "mhc_class": r.mhc_class, "start": int(r.start),
            "peptide": r.peptide, "pred_IC50_nM": round(float(r.pred_IC50_nM), 1),
            "call": r.call, "match_bits": m,
            "bitstring": "".join(map(str, bits)),
            "bitstring_int": bits_to_int(bits),
            "mhc_bitstring_int": bits_to_int(mhc[r.allele]),
        })
    tab = pd.DataFrame(rows).sort_values(["mhc_class", "allele", "pred_IC50_nM"])

    classI = tab[tab.mhc_class == "I"]
    classII = tab[tab.mhc_class == "II"]
    epitope = (classI.iloc[0] if len(classI) else tab.iloc[0])
    peptide = (classII.iloc[0] if len(classII) else tab.iloc[0])

    lines.append(f"# C-ImmSim parameters derived from predicted epitopes of {seq_name}")
    lines.append(f"# generated from IEDB-trained binding models; NBITSTR={NBITSTR}")
    lines.append("")
    lines.append("<1111> RS                  /* random seed */")
    lines.append(f"<{NBITSTR}> NBITSTR            /* bit string length */")
    lines.append(f"<{MINMATCH}> minmatch           /* min complementary bits to bind */")
    lines.append("<1> #microL              /* simulated volume */")
    lines.append("<540> num_step           /* 8 h per step; 540 = ~180 days */")
    lines.append("")
    lines.append("# --- MHC bit strings, numbered per class ---")
    for cls, tag in (("I", "MHCI"), ("II", "MHCII")):
        alleles = sorted(a for a in mhc
                         if tab[tab.allele == a].mhc_class.iloc[0] == cls)
        if not alleles:
            continue
        lines.append(f"<{len(alleles)}> Num{tag}")
        for n, allele in enumerate(alleles, start=1):
            lines.append(f"<{bits_to_int(mhc[allele])}> {tag}{n}   "
                         f"/* {allele} = {''.join(map(str, mhc[allele]))} */")
        lines.append("")
    lines.append("")
    lines.append("# --- antigen: prime + boost ---")
    lines.append("<2> #inj Ag")
    lines.append("<1> Inside               /* peptides per antigen */")
    lines.append("<1> Outside              /* epitopes per antigen */")
    lines.append("<100> Adjuvant           /* LPS / danger signal */")
    lines.append("")
    for t in (2, 86):
        lines.append(f"<{t}> Time to inj ag")
        lines.append("<1000> Num Ag to inj")
        lines.append(f"<{peptide.bitstring_int}> peptide   "
                     f"/* {peptide.peptide} @{peptide.start} {peptide.allele} "
                     f"{peptide.pred_IC50_nM} nM, {peptide.match_bits}/12 match */")
        lines.append(f"<{epitope.bitstring_int}> epitope   "
                     f"/* {epitope.peptide} @{epitope.start} {epitope.allele} "
                     f"{epitope.pred_IC50_nM} nM, {epitope.match_bits}/12 match */")
        lines.append("")

    return tab, "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", help="UniProt accession or FASTA file")
    ap.add_argument("--alleles", default=None, help="comma-separated subset")
    ap.add_argument("--top", type=int, default=25, help="epitopes to print")
    args = ap.parse_args()

    p = Path(args.target)
    if p.exists():
        header, seq = parse_fasta(p.read_text())
        name = p.stem
    else:
        header, seq = fetch_uniprot(args.target)
        name = args.target

    print(f"Antigen : {header[:90]}")
    print(f"Length  : {len(seq)} aa\n")

    wanted = set(args.alleles.split(",")) if args.alleles else None
    try:
        models = load_models(wanted)
    except ValueError as e:                 # a library raises; the CLI exits
        sys.exit(str(e))
    print(f"Alleles : {', '.join(sorted(models))}\n")

    df = scan(seq, models)
    out_csv = RESULTS / f"{name}_epitopes.csv"
    df.to_csv(out_csv, index=False)

    print(f"Scanned {len(df):,} peptide-allele pairs")
    print(f"  strong binders (<{STRONG_NM:.0f} nM): {(df.call=='strong').sum():,}")
    print(f"  weak binders  (<{WEAK_NM:.0f} nM): {(df.call=='weak').sum():,}\n")
    print(f"--- top {args.top} predicted epitopes ---")
    cols = ["allele", "mhc_class", "start", "peptide", "pred_IC50_nM",
            "percentile_rank", "call"]
    print(df[cols].head(args.top).to_string(index=False,
          float_format=lambda v: f"{v:.1f}"))

    tab, text = cimmsim_params(df, header[:60])
    out_par = RESULTS / f"{name}_cimmsim.txt"
    out_par.write_text(text + "\n")
    tab.to_csv(RESULTS / f"{name}_cimmsim_epitopes.csv", index=False)

    print(f"\n--- C-ImmSim bit-string mapping ---")
    print(tab[["allele", "mhc_class", "peptide", "pred_IC50_nM",
               "match_bits", "bitstring", "bitstring_int"]]
          .to_string(index=False, float_format=lambda v: f"{v:.1f}"))
    print(f"\nwrote {out_csv}")
    print(f"wrote {out_par}")


if __name__ == "__main__":
    main()
