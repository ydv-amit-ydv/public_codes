#!/usr/bin/env python3
"""
Harvest MHC mass-spectrometry ELUTED LIGANDS from IEDB, as an independent
validation set for the binding models.

Why this is worth having. The models in `model.py` / `model_pan.py` are
trained on in-vitro binding affinity: a purified MHC molecule, a competitive
radioligand or fluorescence assay, an IC50 in nM. Eluted ligands are a
different measurement entirely - peptides actually recovered from the surface
of a live cell by immunoprecipitation and mass spectrometry. Different labs,
different technology, and different biology, because natural presentation also
requires proteasomal cleavage, TAP transport and peptide loading, none of
which an in-vitro affinity assay sees.

So this is not another held-out split of the same data. It is a test of
whether a model fitted to groove affinity predicts what a cell really
presents. That is the question a vaccine designer is actually asking.

NEGATIVES. Mass spectrometry reports only what it detected, so there are no
measured negatives. The standard benchmark design (NetMHCpan, MHCflurry) is
used here: sample decoy peptides of the same length from the SAME source
proteins the ligands came from, and exclude any peptide observed as a ligand
for any allele. The caveat is real and stated in the output - a decoy is
"probably not presented", not "measured as not presented" - so the absolute
AUC depends on the decoy ratio, and PPV in the top 1% is reported alongside
because it is far less sensitive to it.

    python3 harvest_eluted.py
    python3 harvest_eluted.py --decoy-ratio 19 --max-records 40000
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import urllib.request
import urllib.parse

HERE = Path(__file__).parent
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

IEDB = "https://query-api.iedb.org"
UNIPROT = "https://rest.uniprot.org/uniprotkb/accessions"
PAGE = 5000
UA = "immune-sim-research/1.0"          # generic: no personal detail in headers

COLS = ("linear_sequence,linear_sequence_length,mhc_allele_name,mhc_class,"
        "parent_source_antigen_iri,source_organism_name,pubmed_id,"
        "qualitative_measure,elution_id")


def get(url, tries=4, timeout=120):
    last = None
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:                       # transient proxy / 5xx
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"{url[:110]} failed: {last}")


def harvest_ligands(max_records=60000):
    """Positive eluted ligands, paged. PostgREST needs an explicit order."""
    rows, offset = [], 0
    while len(rows) < max_records:
        q = urllib.parse.urlencode({
            "select": COLS,
            "assay_names": "like.*mass spectrometry*",
            "qualitative_measure": "eq.Positive",
            "linear_sequence_length": "gte.8",
            "order": "elution_id.asc",
            "limit": PAGE, "offset": offset})
        batch = json.loads(get(f"{IEDB}/mhc_search?{q}"))
        if not batch:
            break
        rows.extend(batch)
        offset += PAGE
        print(f"    {len(rows)} ligands", flush=True)
        if len(batch) < PAGE:
            break

    d = pd.DataFrame(rows)
    d = d.rename(columns={"linear_sequence": "peptide",
                          "linear_sequence_length": "length",
                          "mhc_allele_name": "allele"})
    d = d.dropna(subset=["peptide", "allele", "mhc_class"])
    d["peptide"] = d.peptide.str.upper().str.strip()
    # drop anything with a non-standard residue or a modification marker
    ok = d.peptide.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+")
    d = d[ok.fillna(False)]
    # one row per (peptide, allele): the same ligand is reported many times
    d = d.drop_duplicates(subset=["peptide", "allele"])
    d["accession"] = (d.parent_source_antigen_iri.fillna("")
                       .str.replace("UNIPROT:", "", regex=False)
                       .str.split("-").str[0])
    return d.reset_index(drop=True)


def fetch_proteome(accessions, batch=250, cap=4000):
    """
    Sequences of the parent proteins the ligands came from. Decoys are drawn
    from these, so the decoy set has the same amino-acid composition and the
    same protein abundance bias as the positives - which is the point. Drawing
    decoys from an unrelated proteome would make the task easier than it is.
    """
    accs = [a for a in dict.fromkeys(accessions) if a][:cap]
    seqs = {}
    for i in range(0, len(accs), batch):
        chunk = accs[i:i + batch]
        # The documented batch endpoint. The `stream` endpoint rejects a
        # long "accession:A OR accession:B ..." query with HTTP 400.
        q = urllib.parse.urlencode({"accessions": ",".join(chunk),
                                    "format": "fasta"})
        try:
            txt = get(f"{UNIPROT}?{q}", tries=3, timeout=180).decode()
        except Exception as e:
            print(f"    batch {i // batch} failed ({e}); continuing")
            continue
        name, buf = None, []
        for line in txt.splitlines():
            if line.startswith(">"):
                if name and buf:
                    seqs[name] = "".join(buf)
                parts = line[1:].split("|")
                name = parts[1] if len(parts) > 1 else line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name and buf:
            seqs[name] = "".join(buf)
        print(f"    {len(seqs)} sequences", flush=True)
    return seqs


def build_decoys(lig, seqs, ratio=19, seed=0):
    """
    `ratio` decoys per positive, per (allele, length). Decoys are sampled
    uniformly from the parent proteins and screened against every observed
    ligand sequence, whatever allele it came from.
    """
    rng = np.random.default_rng(seed)
    observed = set(lig.peptide)
    pool = [s for s in seqs.values() if len(s) >= 30]
    if not pool:
        return pd.DataFrame()
    lens = np.array([len(s) for s in pool], dtype=float)
    pick = lens / lens.sum()          # long proteins contribute more windows

    out = []
    for (allele, L), grp in lig.groupby(["allele", "length"]):
        L = int(L)
        want = int(len(grp) * ratio)
        mclass = grp.mhc_class.iloc[0]
        seen, guard = set(), 0
        while len(seen) < want and guard < want * 60:
            guard += 1
            s = pool[int(rng.choice(len(pool), p=pick))]
            if len(s) <= L:
                continue
            i = int(rng.integers(0, len(s) - L))
            pep = s[i:i + L]
            if (pep in observed or pep in seen
                    or not set(pep) <= set("ACDEFGHIKLMNPQRSTVWY")):
                continue
            seen.add(pep)
        out.append(pd.DataFrame({"peptide": sorted(seen), "allele": allele,
                                 "length": L, "mhc_class": mclass,
                                 "label": 0}))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--decoy-ratio", type=int, default=19,
                    help="decoys per positive (19 -> 5%% prevalence)")
    ap.add_argument("--max-records", type=int, default=60000)
    ap.add_argument("--min-per-allele", type=int, default=120)
    args = ap.parse_args()

    print("Harvesting MHC mass-spectrometry eluted ligands from IEDB")
    lig = harvest_ligands(args.max_records)
    print(f"  {len(lig)} unique (peptide, allele) ligands, "
          f"{lig.allele.nunique()} alleles, "
          f"{lig.pubmed_id.nunique()} publications")

    keep = lig.allele.value_counts()
    keep = keep[keep >= args.min_per_allele].index
    lig = lig[lig.allele.isin(keep)]
    print(f"  {len(lig)} ligands across {len(keep)} alleles with "
          f">={args.min_per_allele} each")

    print("\nFetching parent protein sequences for decoy sampling")
    seqs = fetch_proteome(lig.accession.values)
    if len(seqs) < 50:
        sys.exit("too few parent sequences fetched to build a decoy set")
    print(f"  {len(seqs)} proteins, "
          f"{sum(len(s) for s in seqs.values()):,} residues")

    print(f"\nSampling {args.decoy_ratio} decoys per positive")
    pos = lig[["peptide", "allele", "length", "mhc_class"]].copy()
    pos["label"] = 1
    neg = build_decoys(lig, seqs, args.decoy_ratio)
    df = pd.concat([pos, neg], ignore_index=True)
    df = df.drop_duplicates(subset=["peptide", "allele"])

    out = DATA / "eluted_ligands.csv"
    df.to_csv(out, index=False)
    meta = {
        "n_positives": int((df.label == 1).sum()),
        "n_decoys": int((df.label == 0).sum()),
        "prevalence": round(float(df.label.mean()), 4),
        "alleles": int(df.allele.nunique()),
        "publications": int(lig.pubmed_id.nunique()),
        "decoy_ratio": args.decoy_ratio,
        "decoy_source": "same parent proteins as the ligands, same lengths",
        "caveat": ("mass spectrometry reports only detections, so decoys are "
                   "unobserved rather than measured negatives; absolute AUC "
                   "depends on the decoy ratio, PPV@1% much less so"),
    }
    (DATA / "eluted_ligands_meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\n  {meta['n_positives']} positives + {meta['n_decoys']} decoys "
          f"-> {out}")
    print(f"  prevalence {meta['prevalence']:.3f}, "
          f"{meta['alleles']} alleles, {meta['publications']} publications")
    print("\n  Independent of the training data in BOTH assay technology "
          "(mass spectrometry\n  from live cells, not purified-MHC affinity) "
          "and publication.")


if __name__ == "__main__":
    main()
