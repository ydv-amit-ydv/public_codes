"""Epitope regions: candidate peptides of both classes grouped by where they sit in the protein.

A CTL (class I) peptide and an HTL (class II) peptide that overlap are one region, labelled CTL+HTL; those are the
usual first choice for a vaccine (one stretch of sequence that both T-cell arms can recognise). The grouping is anchored
on the best peptide of each region, as the overlap collapse in consensus.py: a region cannot creep along the protein.
"""
import numpy as np
import pandas as pd

from common import peptide_props

TYPE_NAMES = {(True, True): "CTL+HTL", (True, False): "CTL only", (False, True): "HTL only"}
DUAL_BONUS = 10.0          # added to the region score of CTL+HTL regions; documented in the report and the README


def _items(pep, cls):
    if pep is None or pep.empty:
        return pd.DataFrame()
    d = pep[pep["Candidate"].astype(bool)].copy()
    d["Class"] = cls
    return d


def build(pep1, pep2, seqs, overlap=0.7):
    """pep1 / pep2: peptide tables from consensus (or None). seqs: {protein: sequence}. -> regions DataFrame."""
    parts = [x for x in (_items(pep1, 1), _items(pep2, 2)) if not x.empty]
    cols = ["Region_Rank", "Region_ID", "Type", "Protein", "Start", "End", "Length", "Sequence", "Region_Score",
            "CTL_Peptides", "CTL_Best_Peptide", "CTL_Best_Allele", "CTL_Coverage", "CTL_Priority_Score",
            "HTL_Peptides", "HTL_Best_Peptide", "HTL_Best_Allele", "HTL_Coverage", "HTL_Priority_Score",
            "In_Panel", "Synthesis_Notes"]
    if not parts:
        return pd.DataFrame(columns=cols)
    it = pd.concat(parts, ignore_index=True).sort_values("Priority_Score", ascending=False, kind="stable").reset_index(drop=True)
    start, end, prot = it["Start"].to_numpy(int), it["End"].to_numpy(int), it["Protein"].to_numpy()
    free = np.ones(len(it), bool)
    groups = []
    for i in range(len(it)):
        if not free[i]:
            continue
        a0, a1 = start[i], end[i]
        ov = np.minimum(end, a1) - np.maximum(start, a0) + 1
        shorter = np.minimum(end - start + 1, a1 - a0 + 1)
        m = free & (prot == prot[i]) & (ov > 0) & (ov / shorter >= overlap)
        m[i] = True
        free &= ~m
        groups.append(np.flatnonzero(m))

    rows = []
    for idx in groups:
        g = it.iloc[idx]
        r = {"Protein": g["Protein"].iloc[0], "Start": int(g["Start"].min()), "End": int(g["End"].max())}
        r["Length"] = r["End"] - r["Start"] + 1
        has = {c: bool((g["Class"] == c).any()) for c in (1, 2)}
        r["Type"] = TYPE_NAMES[(has[1], has[2])]
        sc = []
        for c, tag in ((1, "CTL"), (2, "HTL")):
            gc = g[g["Class"] == c]
            r[f"{tag}_Peptides"] = int(len(gc))
            if gc.empty:
                continue
            best = gc.iloc[0]
            r[f"{tag}_Best_Peptide"], r[f"{tag}_Best_Allele"] = best["Peptide"], best.get("Best_Allele", "")
            r[f"{tag}_Coverage"] = float(gc["Coverage_Candidate"].max())
            r[f"{tag}_Priority_Score"] = float(best["Priority_Score"])
            sc.append(float(best["Priority_Score"]))
        r["Region_Score"] = round(float(np.mean(sc)) + (DUAL_BONUS if len(sc) == 2 else 0.0), 2)
        r["In_Panel"] = bool(g["In_Panel"].astype(bool).any()) if "In_Panel" in g else False
        seq = (seqs or {}).get(r["Protein"], "")
        r["Sequence"] = seq[r["Start"] - 1:r["End"]] if seq else ""
        rows.append(r)
    out = pd.DataFrame(rows)
    out["_dual"] = (out["Type"] == "CTL+HTL").astype(int)
    out = out.sort_values(["_dual", "Region_Score"], ascending=False, kind="stable").drop(columns="_dual").reset_index(drop=True)
    out.insert(0, "Region_Rank", range(1, len(out) + 1))
    out.insert(1, "Region_ID", [f"R{i:03d}" for i in out["Region_Rank"]])
    seqs_ok = out["Sequence"].str.len().between(5, 60) & ~out["Sequence"].str.contains(r"[^ACDEFGHIKLMNPQRSTVWY]")
    notes = pd.Series("", index=out.index)
    if seqs_ok.any():
        d = peptide_props.describe(out.loc[seqs_ok, "Sequence"].tolist())
        notes[seqs_ok] = d["Synthesis_Notes"].to_numpy()
    out["Synthesis_Notes"] = notes
    for c in cols:
        if c not in out:
            out[c] = np.nan
    return out[cols]


def label_peptides(pep, regions, cls):
    """Adds Region_ID / Region_Type to a peptide table (by overlap with a region of the same protein)."""
    if pep is None or pep.empty or regions.empty:
        return pep
    pep = pep.copy()
    rid, rtype = [""] * len(pep), [""] * len(pep)
    by_prot = {p: g for p, g in regions.groupby("Protein")}
    for k, (p, s, e) in enumerate(zip(pep["Protein"], pep["Start"], pep["End"])):
        g = by_prot.get(p)
        if g is None:
            continue
        hit = g[(g["Start"] <= e) & (g["End"] >= s)]
        if len(hit):
            rid[k], rtype[k] = hit["Region_ID"].iloc[0], hit["Type"].iloc[0]
    pep["Region_ID"], pep["Region_Type"] = rid, rtype
    return pep
