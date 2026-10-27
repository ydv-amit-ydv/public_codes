"""One form in, one table out: NetMHCpan (BA/EL) -> BigMHC (EL/IM) -> filters -> summary.

Chaining: with `chain` on, BigMHC scores only the peptide-allele pairs that passed the
NetMHCpan filters. BigMHC is the slow model (it reloads its 7 networks on every call),
so this typically cuts the run time by an order of magnitude.

Import after engine_server has put the tool folders on sys.path.
"""
import os, time

import pandas as pd

from common import alleles as allele_lib, config
import bigmhc_runner as bm
import netmhcpan_app as nm

MAX_PAIRS = 3_000_000
FILTER_KEYS = ("rank_el_max", "rank_ba_max", "aff_max", "bm_el_min", "bm_im_min")


def _num(x):
    if x in (None, ""):
        return None
    return float(x)


def normalise(params):
    """Fill defaults and validate; raises ValueError with a message meant for the user."""
    p = dict(params or {})
    models = {"netmhcpan": False, "bigmhc_el": False, "bigmhc_im": False}
    models.update({k: bool(v) for k, v in (p.get("models") or {}).items() if k in models})
    if not any(models.values()):
        raise ValueError("Select at least one model.")
    alleles = list(dict.fromkeys(allele_lib.canonical(a) for a in p.get("alleles", []) if str(a).strip()))
    if not alleles:
        raise ValueError("Select at least one allele.")
    lengths = sorted({int(x) for x in p.get("lengths", [9]) if 8 <= int(x) <= 14})
    if p.get("mode", "scan") == "scan" and not lengths:
        raise ValueError("Select at least one peptide length (8-14).")
    f = {k: _num((p.get("filters") or {}).get(k)) for k in FILTER_KEYS}
    nmp = {"ba": bool((p.get("nmp") or {}).get("ba", True)),
           "rankS": float((p.get("nmp") or {}).get("rankS", 0.5)),
           "rankW": float((p.get("nmp") or {}).get("rankW", 2.0))}
    return {"text": p.get("text", ""), "mode": p.get("mode", "scan"), "lengths": lengths,
            "alleles": alleles, "models": models, "chain": bool(p.get("chain", True)),
            "filters": f, "nmp": nmp}


def _supported(alleles):
    """Which alleles each model can score, plus human-readable notes about what is skipped."""
    notes = []
    known = set(nm.list_alleles())
    nm_ok = [a for a in alleles if not known or allele_lib.netmhcpan_name(a) in known]
    bm_ok = [a for a in alleles if a.startswith("HLA-")]
    if len(nm_ok) < len(alleles):
        notes.append("NetMHCpan does not know: " + ", ".join(a for a in alleles if a not in nm_ok))
    if len(bm_ok) < len(alleles):
        notes.append("BigMHC is human HLA-I only, skipping: " + ", ".join(a for a in alleles if a not in bm_ok))
    return nm_ok, bm_ok, notes


def preview(params):
    p = normalise(params)
    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"])
    nm_ok, bm_ok, notes = _supported(p["alleles"])
    n_pep = int(peps["Peptide"].nunique()) if not peps.empty else 0
    m = p["models"]
    out = {"n_peptides": n_pep, "n_alleles": len(p["alleles"]), "warnings": warns + notes,
           "n_pairs_netmhcpan": n_pep * len(nm_ok) if m["netmhcpan"] else 0,
           "n_pairs_bigmhc": n_pep * len(bm_ok) if (m["bigmhc_el"] or m["bigmhc_im"]) else 0,
           "chained": p["chain"] and m["netmhcpan"] and (m["bigmhc_el"] or m["bigmhc_im"])}
    out["too_big"] = max(out["n_pairs_netmhcpan"], out["n_pairs_bigmhc"]) > MAX_PAIRS
    return out


def _nmp_pass(df, f):
    ok = pd.Series(True, index=df.index)
    if f["rank_el_max"] is not None:
        ok &= df["NMP_%Rank_EL"] <= f["rank_el_max"]
    if f["rank_ba_max"] is not None and "NMP_%Rank_BA" in df:
        ok &= df["NMP_%Rank_BA"] <= f["rank_ba_max"]
    if f["aff_max"] is not None and "NMP_Aff_nM" in df:
        ok &= df["NMP_Aff_nM"] <= f["aff_max"]
    return ok


def run(params, job):
    t0 = time.time()
    p = normalise(params)
    m, f = p["models"], p["filters"]
    use_bm = m["bigmhc_el"] or m["bigmhc_im"]
    bm_models = [x for x, on in (("el", m["bigmhc_el"]), ("im", m["bigmhc_im"])) if on]

    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"])
    if peps.empty:
        raise ValueError("No valid peptides in the input.")
    peps = peps.drop_duplicates()
    nm_ok, bm_ok, notes = _supported(p["alleles"])
    warns += notes
    n_pep = int(peps["Peptide"].nunique())
    if n_pep * len(p["alleles"]) > MAX_PAIRS:
        raise ValueError(f"{n_pep:,} peptides x {len(p['alleles'])} alleles is over the {MAX_PAIRS:,} limit; "
                         "use fewer alleles or lengths, or split the sequence.")
    job.log(f"{n_pep:,} unique peptides, {len(p['alleles'])} alleles")

    if m["netmhcpan"]:
        job.stage("netmhcpan", "NetMHCpan 4.2 (EL / BA)", weight=1.0)
    for x in bm_models:
        job.stage(f"bigmhc_{x}", f"BigMHC {x.upper()}", weight=2.0)
    job.stage("merge", "Merging and filtering", weight=0.3)

    # ---- 1. NetMHCpan ---------------------------------------------------------------
    nmp = None
    if m["netmhcpan"]:
        if not nm_ok:
            raise ValueError("None of the selected alleles is known to NetMHCpan.")
        res, logs = nm.run_netmhcpan(peps["Peptide"].tolist(), [allele_lib.netmhcpan_name(a) for a in nm_ok],
                                     {**p["nmp"], "t": None, "raw": ""}, job=job)
        for l in logs[-3:]:
            job.log(l[:300])
        if res.empty:
            raise RuntimeError("NetMHCpan returned no results (see the log).")
        res["Allele"] = [allele_lib.canonical(a) for a in res["Allele"]]
        nmp = res
        job.finish_stage("netmhcpan", f"{len(res):,} predictions")

    # ---- 2. BigMHC (on NetMHCpan survivors when chained) ---------------------------------
    bmres = None
    if use_bm:
        if nmp is not None and p["chain"]:
            keep = nmp[_nmp_pass(nmp, f)]
            pairs = keep[["Peptide", "Allele"]].drop_duplicates()
            pairs = pairs[pairs["Allele"].isin(bm_ok)]
            job.log(f"chain: {len(pairs):,} of {len(nmp):,} NetMHCpan rows pass the filters -> BigMHC")
        else:
            pairs = pd.DataFrame({"Peptide": peps["Peptide"].drop_duplicates()}).merge(
                pd.DataFrame({"Allele": bm_ok}), how="cross")
        if pairs.empty:
            for x in bm_models:
                job.finish_stage(f"bigmhc_{x}", "nothing to score")
            warns.append("BigMHC had nothing to score (no pair passed the earlier filters, or no human alleles).")
        else:
            bmres = bm.predict_pairs(pairs.rename(columns={"Peptide": "peptide", "Allele": "allele"}),
                                     bm_models, job=job).rename(columns={"peptide": "Peptide", "allele": "Allele"})

    # ---- 3. merge, flag, summarise -----------------------------------------------------------
    job.update("merge", msg="merging", state="running")
    if nmp is not None:
        base = peps.merge(nmp, on="Peptide")
    else:
        base = peps.merge(pd.DataFrame({"Allele": bm_ok}), how="cross")
    if bmres is not None:
        base = base.merge(bmres, on=["Peptide", "Allele"], how="left")
    base.insert(3, "Length", base["Peptide"].str.len())

    flags = []
    if nmp is not None:
        base["Pass_NMP"] = _nmp_pass(base, f)
        flags.append("Pass_NMP")
    elif any(f[k] is not None for k in ("rank_el_max", "rank_ba_max", "aff_max")):
        warns.append("NetMHCpan filters ignored: NetMHCpan was not run.")
    if bmres is not None:
        ok = pd.Series(True, index=base.index)
        if f["bm_el_min"] is not None and "BigMHC_EL" in base:
            ok &= base["BigMHC_EL"] >= f["bm_el_min"]
        if f["bm_im_min"] is not None and "BigMHC_IM" in base:
            ok &= base["BigMHC_IM"] >= f["bm_im_min"]
        base["Pass_BigMHC"] = ok & base[[c for c in base if c.startswith("BigMHC_")]].notna().all(axis=1)
        flags.append("Pass_BigMHC")
    base["Pass_All"] = base[flags].all(axis=1) if flags else pd.Series(True, index=base.index)

    sort_cols, asc = ["Pass_All"], [False]
    for col, a in (("NMP_%Rank_EL", True), ("BigMHC_IM", False)):
        if col in base:
            sort_cols.append(col); asc.append(a)
    base = base.sort_values(sort_cols, ascending=asc, na_position="last", kind="stable").reset_index(drop=True)

    keys = ["Protein", "Start", "Peptide"]
    g = base.groupby(keys, sort=False)
    summ = g.agg(Length=("Length", "first"), Alleles_tested=("Allele", "nunique"),
                 Alleles_passed=("Pass_All", "sum")).reset_index()
    passed = base[base["Pass_All"]].groupby(keys, sort=False)["Allele"].agg(lambda a: ", ".join(sorted(a)))
    summ = summ.merge(passed.rename("Alleles_passed_list").reset_index(), on=keys, how="left")
    summ["Alleles_passed_list"] = summ["Alleles_passed_list"].fillna("")
    for col, fn, name in (("NMP_%Rank_EL", "min", "Best_NMP_%Rank_EL"), ("NMP_Aff_nM", "min", "Best_NMP_Aff_nM"),
                          ("BigMHC_EL", "max", "Max_BigMHC_EL"), ("BigMHC_IM", "max", "Max_BigMHC_IM")):
        if col in base:
            summ = summ.merge(getattr(g[col], fn)().rename(name).reset_index(), on=keys, how="left")
    order, asc = ["Alleles_passed"], [False]
    for col, a_ in (("Best_NMP_%Rank_EL", True), ("Max_BigMHC_IM", False)):
        if col in summ:
            order.append(col); asc.append(a_)
    summ = summ.sort_values(order, ascending=asc, na_position="last", kind="stable").reset_index(drop=True)

    jobdir = os.path.join(config.tmp_dir(), "jobs"); os.makedirs(jobdir, exist_ok=True)
    pairs_csv, summ_csv = (os.path.join(jobdir, f"{job.id}_{n}.csv") for n in ("pairs", "summary"))
    base.round(6).to_csv(pairs_csv, index=False)
    summ.round(6).to_csv(summ_csv, index=False)
    job.finish_stage("merge", "done")

    job.result = {"pairs_csv": pairs_csv, "summary_csv": summ_csv, "summary": {
        "peptides": n_pep, "pairs": int(len(base)), "pairs_passed": int(base["Pass_All"].sum()),
        "peptides_passed": int((summ["Alleles_passed"] > 0).sum()),
        "nmp_pass": int(base["Pass_NMP"].sum()) if "Pass_NMP" in base else None,
        "seconds": round(time.time() - t0, 1), "warnings": warns,
        "chained": bool(nmp is not None and use_bm and p["chain"])}}
    job.log(f"done in {time.time() - t0:.1f}s: {int(base['Pass_All'].sum()):,} pairs / "
            f"{int((summ['Alleles_passed'] > 0).sum()):,} peptides pass")
