"""One request: an antigen in, a simulated immune response out.

The antigen can come from three places, in descending order of how much you
should trust the epitopes:

  engine    a finished engine run. Eleven published predictors, fused, with
            population-weighted HLA coverage. This is the best input there
            is, and it is the whole reason the simulator lives in the engine.
  sequence  a pasted protein, scanned by immune_sim's own IEDB-trained
            models. Convenient, and clearly weaker: on a pathogen it has
            never seen those models reach 0.873 ROC-AUC for class I and 0.748
            for class II, where class II is barely better than counting amino
            acids.
  demo      the committed SARS-CoV-2 spike epitope set, for trying the
            simulator without running anything first.

Both halves run in one call because joining them by hand through a CSV on
disk was the thing that made this two tools instead of one.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
SIM = os.path.join(HERE, "immune_sim")
for p in (SIM, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LEN = 4000


def _results_dir():
    d = os.path.join(SIM, "results")
    os.makedirs(d, exist_ok=True)
    return d


def _clean_sequence(raw, name):
    """FASTA or bare residues -> an amino acid string, or a clear refusal."""
    header = ""
    if ">" in raw:
        lines = raw.splitlines()
        header = next((l for l in lines if l.startswith(">")), "")
        raw = "\n".join(l for l in lines if not l.startswith(">"))
    seq = "".join(c for c in raw.upper() if c.isalpha())
    if len(seq) < 9:
        raise ValueError("Paste at least 9 amino acid residues.")
    if len(seq) > MAX_LEN:
        raise ValueError(f"That sequence is {len(seq):,} residues; the limit "
                         f"is {MAX_LEN:,}. Paste one protein, not a proteome.")
    bad = sorted(set(seq) - AA)
    if bad:
        raise ValueError("These are not amino acids: " + ", ".join(bad)
                         + ". Is this a DNA or RNA sequence?")
    # A, C, G, T and U are ALL valid amino acid letters (Ala, Cys, Gly, Thr,
    # Sec), so a pasted gene passes the check above and would be silently
    # scanned as a protein. Catch it by composition instead.
    if len(seq) >= 24 and set(seq) <= set("ACGTU"):
        kind = "RNA" if "U" in seq and "T" not in seq else "DNA"
        raise ValueError(
            f"This looks like a {kind} sequence, not a protein: every one of "
            f"its {len(seq)} residues is A, C, G, T or U. Those are also "
            f"valid amino acid letters, so nothing here can tell them apart "
            f"for you. Translate it to protein first, then paste that.")
    if not name and header:
        name = header.lstrip(">").split("|")[-1].strip()[:40]
    return seq, name


def _from_engine(run_id, name):
    """Epitopes from a finished engine run, through engine_adapter."""
    from common import config
    import engine_adapter
    d = config.run_dir(run_id)
    if not os.path.isdir(d):
        raise ValueError(f"No engine run {run_id!r} in the workspace.")
    try:
        ep = engine_adapter.load(d)
    except Exception as e:                                   # noqa: BLE001
        raise ValueError(
            f"Could not read epitopes from that run: {e}. A run has to have "
            f"finished and produced a pair table before it can be simulated.")
    if ep is None or not len(ep):
        raise ValueError("That engine run produced no epitopes to simulate.")
    return ep, (name or f"engine run {run_id[:8]}")


def _from_sequence(seq, alleles, name):
    """Epitopes from immune_sim's own models."""
    import predict as PR
    models = PR.load_models(alleles or None)
    if not models:
        raise ValueError(
            "The simulation models are not installed yet. Settings -> Immune "
            "simulation -> Install trains them from the data committed in "
            "this repository; it needs no network.")
    pred = PR.scan(seq, models)
    binders = pred[pred.call != "non-binder"].sort_values("pred_IC50_nM")
    if binders.empty:
        raise ValueError(
            "No binders were predicted for this sequence and these alleles, "
            "so there is nothing to simulate. Try more alleles or a longer "
            "sequence.")
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
    binders.to_csv(os.path.join(_results_dir(), f"{safe}_epitopes.csv"),
                   index=False)
    return binders, name


def design(payload):
    import platforms as PF
    from dynamics import Sim, load_epitopes

    source = (payload.get("antigen_source") or "sequence").strip()
    name = (payload.get("name") or "").strip()[:40]
    seq = ""
    table = None
    provenance = ""

    if source == "engine":
        ep_raw, name = _from_engine(str(payload.get("run_id") or ""), name)
        table = ep_raw
        provenance = ("epitopes from an engine run: eleven predictors fused, "
                      "with population-weighted HLA coverage")
    elif source == "demo":
        name = name or "demo (SARS-CoV-2 spike)"
        f = os.path.join(SIM, "results", "demo_epitopes.csv")
        if not os.path.exists(f):
            raise ValueError("results/demo_epitopes.csv is missing.")
        table = pd.read_csv(f)
        provenance = "the committed demo epitope set"
    else:
        seq, name = _clean_sequence(payload.get("sequence") or "", name)
        name = name or "query"
        table, name = _from_sequence(seq, payload.get("alleles"), name)
        provenance = ("epitopes from immune_sim's own IEDB-trained models. "
                      "An engine run would be a better input: on an unseen "
                      "pathogen these reach 0.873 ROC-AUC for class I and "
                      "0.748 for class II")

    ep = load_epitopes(table)

    # --- parameters -------------------------------------------------------
    import json as _json
    pfile = os.path.join(SIM, "models", "calibrated_wetlab.json")
    params, pname = None, "literature defaults (NOT validated)"
    if payload.get("calibration", "wetlab") != "none" and os.path.exists(pfile):
        with open(pfile) as fh:
            params = _json.load(fh)["params"]
        pname = "calibrated_wetlab.json"
    if payload.get("replicating"):
        params = dict(params or {}, ag_growth=0.55)

    days = float(payload.get("days", 240) or 240)
    n_steps = int(max(30, min(days * 3, 3000)))
    res_h = float(payload.get("resolution_h", 8) or 8)
    inj_days = payload.get("inject_days") or [0, 90]
    adj = payload.get("adjuvant") or None
    inj = sorted({(max(2, int(round(float(d) * 3))) if float(d) > 0 else 2,
                   None, adj) for d in inj_days})

    sim = Sim(ep, p=params, injections=inj,
              platform=payload.get("platform", "mrna_lnp"),
              adjuvant=adj,
              dose_ug=(float(payload["dose_ug"])
                       if payload.get("dose_ug") else None),
              heterologous=payload.get("heterologous") or None)
    df = sim.run(n_steps=n_steps, out_hours=max(0.5, min(res_h, 48.0)))

    # --- shape the answer -------------------------------------------------
    keep = [c for c in ("antigen", "antigen_depot", "template", "IgM", "IgG",
                        "affinity", "TH_active", "TH_memory", "TC_active",
                        "TC_memory", "B_active", "B_memory", "plasma",
                        "plasma_longlived", "danger")
            if c in df.columns]
    step = max(1, len(df) // 900)
    d = df.iloc[::step]
    inj_day_list = [round(s * 8 / 24, 2) for s, _, _ in inj]
    boost = max(inj_day_list) if len(inj_day_list) > 1 else None
    summary = {
        "peak_IgG": float(df.IgG.max()),
        "peak_IgG_day": float(df.day[df.IgG.idxmax()]),
        "peak_CTL": float(df.TC_active.max()),
        "peak_B": float(df.B_active.max()),
        "B_memory_end": float(df.B_memory.iloc[-1]),
        "CTL_memory_end": float(df.TC_memory.iloc[-1]),
        "affinity_fold": float(df.affinity.max() / max(df.affinity.iloc[0], 1e-9)),
        "igg_durability": float(df.IgG.iloc[-1] / max(df.IgG.max(), 1e-9)),
    }
    if boost:
        pre, post = df[df.day < boost], df[df.day >= boost]
        if len(pre) and len(post):
            summary["boost_fold_over_baseline"] = float(
                post.IgG.max() / max(pre.IgG.iloc[-1], 1e-9))
            summary["boost_peak_day_after"] = float(
                df.day[post.IgG.idxmax()] - boost)

    cols = [c for c in ("allele", "mhc_class", "start", "peptide", "length",
                        "pred_IC50_nM", "percentile_rank", "call")
            if c in table.columns]
    top = table[cols].head(int(payload.get("limit", 80)))
    for c in ("pred_IC50_nM", "percentile_rank"):
        if c in top.columns:
            top[c] = top[c].round(3)

    return {
        "antigen": {"name": name, "length": len(seq) or None,
                    "source": source, "provenance": provenance},
        "epitopes": {
            "columns": cols,
            "rows": top.where(pd.notna(top), None).values.tolist(),
            "n": int(len(table)),
            "n_alleles": int(table.allele.nunique()) if "allele" in table else 0,
            "by_class": (table.groupby("mhc_class").size().to_dict()
                         if "mhc_class" in table else {}),
            "n_clones": int(len(ep)),
        },
        "day": d.day.round(3).tolist(),
        "series": {k: np.nan_to_num(d[k].values, posinf=0, neginf=0)
                    .round(6).tolist() for k in keep},
        "summary": summary,
        "info": {"params": pname, "injection_days": inj_day_list,
                 "platform": sim.pf.label, "adjuvant": sim.adj.label},
        "describe": PF.describe(sim.platform_key, sim.adj.key, sim.dose_ug),
    }
