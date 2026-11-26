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
    # A, C, G, T and U are all valid amino acid letters, so a pasted gene
    # looks like a protein. Real nucleotide sequence uses the whole alphabet,
    # though: requiring at least three distinct letters keeps this from
    # firing on a low-complexity protein such as a poly-alanine control,
    # which it did.
    nt = set(seq) <= set("ACGTU") and len(set(seq)) >= 3
    if len(seq) >= 24 and nt:
        kind = "RNA" if "U" in seq and "T" not in seq else "DNA"
        raise ValueError(
            f"This looks like a {kind} sequence, not a protein: all "
            f"{len(seq)} residues are A, C, G, T or U and it uses "
            f"{len(set(seq))} of them. Those are also valid amino acid "
            f"letters, so nothing here can tell them apart for you. "
            f"Translate it to protein first, then paste that.")
    if not name and header:
        name = header.lstrip(">").split("|")[-1].strip()[:40]
    return seq, name


def _from_engine(run_id, name):
    """Epitopes from a finished engine run, through engine_adapter."""
    from common import config
    import engine_adapter
    # NOT config.run_dir(): that calls os.makedirs, so merely asking about a
    # run id creates a directory for it. A read path must not write, and a
    # browser sending a stale or mistyped id was littering the workspace with
    # empty run folders.
    rid = str(run_id).split("_")[0]
    if not rid or rid in (".", "..") or "/" in rid or "\\" in rid:
        raise ValueError("That is not a run id.")
    d = os.path.join(config.workspace(), "runs", rid)
    if not os.path.isdir(d):
        raise ValueError(f"No engine run {rid!r} in the workspace.")
    try:
        ep = engine_adapter.load(d)
    except Exception as e:                                   # noqa: BLE001
        raise ValueError(
            f"Could not read epitopes from that run: {e}. A run has to have "
            f"finished and produced a pair table before it can be simulated.")
    if ep is None or not len(ep):
        raise ValueError("That engine run produced no epitopes to simulate.")
    return ep, (name or f"engine run {run_id[:8]}")


def _from_sequence(seq, alleles, name, predictor="per_allele"):
    """Epitopes from immune_sim's own models.

    `predictor` picks which one. The pan-allele core-scan model is the
    better of the two on a pathogen it has never seen - median ROC-AUC
    0.873 against 0.766 for class I - but it scores every length at every
    register, so it is several times slower on a long protein. The choice
    is offered rather than imposed, and whichever ran is reported back.
    """
    import pandas as _pd
    import predict as PR

    def per_allele(classes=None):
        try:
            models = PR.load_models(alleles or None)
        except ValueError:
            models = {}
        if classes is not None:
            models = {a: m for a, m in models.items()
                      if m["mhc_class"] in classes}
        if not models:
            return None
        return PR.scan(seq, models)

    pred, used = None, "per-allele"
    if predictor == "pan":
        try:
            pan = PR.scan_pan(seq, alleles or None)
        except Exception as e:                               # noqa: BLE001
            pan, used = None, f"per-allele (the pan model failed: {e})"
        if pan is not None and len(pan):
            # Fill in any class the pan model does not cover from the
            # per-allele models. Without this, a half-trained pan model
            # silently returns class I only - and a class I-only epitope set
            # gives NO CD4 help and therefore zero antibody, which looks
            # like a result rather than a missing model. It cost 86 seconds
            # and an answer of 0 to find that.
            have = set(pan.mhc_class.unique())
            missing = {"I", "II"} - have
            filled = per_allele(missing) if missing else None
            if filled is not None and len(filled):
                pan = _pd.concat([pan, filled], ignore_index=True)
                used = ("pan-allele core-scan for class "
                        + "/".join(sorted(have))
                        + ", per-allele for class "
                        + "/".join(sorted(missing))
                        + " (the pan model is not trained for it)")
            else:
                used = "pan-allele core-scan"
            pred = pan
    if pred is None:
        pred = per_allele()
        if pred is None:
            raise ValueError(
                "The simulation models are not installed yet. Settings -> "
                "Immune simulation -> Install trains them from the data "
                "committed in this repository; it needs no network.")
    binders = pred[pred.call != "non-binder"].sort_values("pred_IC50_nM")
    if binders.empty:
        raise ValueError(
            "No binders were predicted for this sequence and these alleles, "
            "so there is nothing to simulate. Try more alleles or a longer "
            "sequence.")
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in name)
    binders.to_csv(os.path.join(_results_dir(), f"{safe}_epitopes.csv"),
                   index=False)
    return binders, name, used


def _add_host_self(table):
    """Mark each epitope as identical / near / distinct to the host proteome.

    This is the engine's own host-proteome check, and it is what lets the
    simulation tell a self protein from a foreign one. Without it a human
    housekeeping protein simulates as a strong immunogen, which is the defect
    reported for C-ImmSim and which this model shared.

    Returns (table, note). The note is shown to the user: if the proteome has
    not been downloaded the simulation still runs, but it runs WITHOUT
    tolerance, and that has to be said rather than quietly assumed away.
    """
    try:
        from common import hostself
        if not hostself.available():
            return table, ("the host proteome is not downloaded, so central "
                           "tolerance was NOT applied: a self protein will "
                           "simulate as if it were foreign. Settings \u2192 "
                           "Host proteome downloads it, about 8 MB.")
        ix = hostself.index()
    except Exception as e:                                   # noqa: BLE001
        return table, (f"the host-proteome check could not run ({e}), so "
                       f"central tolerance was NOT applied.")
    table = table.copy()
    table["host_self"] = [hostself.check(str(p), ix)["Host_Self"]
                          for p in table.peptide]
    n_self = int((table.host_self == "identical").sum())
    frac = n_self / max(len(table), 1)
    note = (f"{n_self:,} of {len(table):,} epitopes ({frac:.0%}) are identical "
            f"to a host protein and are deleted by central tolerance")
    return table, note


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
    elif source == "uniprot":
        acc = "".join(c for c in str(payload.get("accession", ""))
                      if c.isalnum() or c in "-_")[:20]
        if not acc:
            raise ValueError("Give a UniProt accession, e.g. P0DTC2.")
        import predict as PR
        try:
            fetched, seq = PR.fetch_uniprot(acc)
        except Exception as e:                               # noqa: BLE001
            raise ValueError(f"Could not fetch {acc} from UniProt: {e}")
        name = name or (fetched.split("|")[-1].strip()[:40] or acc)
        table, name, used = _from_sequence(seq, payload.get("alleles"), name,
                                           payload.get("predictor", "per_allele"))
        provenance = (f"{acc} fetched from UniProt, then scanned with "
                      f"immune_sim's own models \u2014 {used}")
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
        table, name, used = _from_sequence(seq, payload.get("alleles"), name,
                                           payload.get("predictor", "per_allele"))
        provenance = (f"epitopes from immune_sim's own models \u2014 {used}. "
                      "An engine run would be a better input: on an unseen "
                      "pathogen these reach 0.873 ROC-AUC for class I and "
                      "0.748 for class II")

    table, tol_note = _add_host_self(table)
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

    # Show the numbers the source actually produced. An engine run's native
    # scale is percentile rank with agreement and coverage behind it; showing
    # only the nM that was derived from it would hide the real measurement and
    # make a fused eleven-predictor call look like a single affinity estimate.
    cols = [c for c in ("allele", "mhc_class", "start", "peptide", "length",
                        "bind_rank", "agreement", "coverage", "priority",
                        "pred_IC50_nM", "percentile_rank", "weight", "call")
            if c in table.columns]
    top = table[cols].head(int(payload.get("limit", 80)))
    for c in ("pred_IC50_nM", "percentile_rank", "bind_rank", "agreement",
              "coverage", "priority", "weight"):
        if c in top.columns:
            top[c] = top[c].round(3)

    return {
        "antigen": {"name": name, "length": len(seq) or None,
                    "source": source, "provenance": provenance,
                    "tolerance": tol_note,
                    "tolerance_applied": bool(ep.attrs.get("tolerance_applied")),
                    "deleted_by_tolerance":
                        int(ep.attrs.get("clones_deleted_by_tolerance", 0))},
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
