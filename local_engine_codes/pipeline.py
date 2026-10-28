"""One form in, one table out.

    NetMHCpan (BA/EL) -> MixMHCpred (presentation) -> PRIME (immunogenicity) -> BigMHC (EL/IM)

NetMHCpan, MixMHCpred and PRIME are quick and run on every peptide-allele pair. BigMHC is the
slow one (it reloads its 7 networks on every call), so with `chain` on it scores only the pairs
that pass the *binding* filters (NetMHCpan and/or MixMHCpred). That typically cuts the run time
by an order of magnitude.

Import after engine_server has put the tool folders on sys.path.
"""
import os, time

import pandas as pd

from common import alleles as allele_lib, config
import bigmhc_runner as bm
import calis_immunogenicity as ci
import iedb_processing as ip
import mhcflurry_runner as mf
import mixmhcpred_runner as mx
import netmhcpan_app as nm
import prime_runner as pr
import consensus
import funnel
import report

MAX_PAIRS = 3_000_000
MODEL_KEYS = ("netmhcpan", "mixmhcpred", "mhcflurry", "prime", "iedb_imm", "netctlpan", "netchop", "bigmhc_el", "bigmhc_im")
FILTER_KEYS = ("rank_el_max", "rank_ba_max", "aff_max", "mix_rank_max", "prime_rank_max",
               "iedb_imm_min", "mhcf_pres_rank_max", "mhcf_aff_max", "ctlpan_rank_max", "netchop_min",
               "bm_el_min", "bm_im_min")


def _num(x):
    return None if x in (None, "") else float(x)


def normalise(params):
    """Fill defaults and validate; raises ValueError with a message meant for the user."""
    p = dict(params or {})
    fun = funnel.normalise(p.get("funnel"), 1, p.get("mode", "scan"))
    models = {k: False for k in MODEL_KEYS}
    if fun["mode"] == "manual":
        models.update({m: True for st in fun["stages"] for m in st["models"]})
    else:
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
    if p.get("mode", "scan") == "peptides" and (models["netctlpan"] or models["netchop"]):
        raise ValueError("NetCTLpan and NetChop need the protein the peptides come from (cleavage depends on the "
                         "surrounding residues): choose 'Cut a protein / FASTA into fragments'.")
    iedb = {"mask_choice": (p.get("iedb") or {}).get("mask_choice", "default"),
            "positions": str((p.get("iedb") or {}).get("positions", ""))}
    if iedb["mask_choice"] not in ci.MASK_CHOICES:
        raise ValueError("Positions to mask must be default, by_allele or custom.")
    if models["iedb_imm"] and iedb["mask_choice"] == "custom":
        ci.parse_positions(iedb["positions"])                  # ValueError with a readable message
    return {"text": p.get("text", ""), "mode": p.get("mode", "scan"), "lengths": lengths,
            "alleles": alleles, "models": models, "chain": bool(p.get("chain", True)) and fun["mode"] == "auto", "funnel": fun,
            "filters": f, "nmp": nmp, "iedb": iedb,
            "preset": p.get("preset") or "balanced", "rules": p.get("rules") or {}}


def _safe_supported(fn, home_key, alleles, tool):
    """Alleles a tool knows; (all, []) if the tool is not installed (the run itself reports that)."""
    try:
        home = config.model(home_key).get("home")
        if not home:
            return list(alleles), []
        known = fn(home)
    except OSError:
        return list(alleles), []
    ok = [a for a in alleles if (allele_lib.gfeller_name(a) or "") in known]
    return ok, ([f"{tool} does not support: " + ", ".join(a for a in alleles if a not in ok)] if len(ok) < len(alleles) else [])


def _supported(alleles, m):
    """Per-model supported alleles, plus notes about what each model will skip."""
    notes = []
    known = set(nm.list_alleles())
    nm_ok = [a for a in alleles if not known or allele_lib.netmhcpan_name(a) in known]
    bm_ok = [a for a in alleles if a.startswith("HLA-")]
    if m["netmhcpan"] and len(nm_ok) < len(alleles):
        notes.append("NetMHCpan does not know: " + ", ".join(a for a in alleles if a not in nm_ok))
    if (m["bigmhc_el"] or m["bigmhc_im"]) and len(bm_ok) < len(alleles):
        notes.append("BigMHC is human HLA-I only, skipping: " + ", ".join(a for a in alleles if a not in bm_ok))
    if m["mixmhcpred"]:
        notes += _safe_supported(mx.supported, "mixmhcpred", alleles, "MixMHCpred")[1]
    if m["prime"]:
        notes += _safe_supported(pr.supported, "prime", alleles, "PRIME")[1]
    return nm_ok, bm_ok, notes


def preview(params):
    p = normalise(params)
    m = p["models"]
    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"])
    nm_ok, bm_ok, notes = _supported(p["alleles"], m)
    n_pep = int(peps["Peptide"].nunique()) if not peps.empty else 0
    n_win = int(len(peps))
    use_bm = m["bigmhc_el"] or m["bigmhc_im"]
    fast = (m["mixmhcpred"] or m["prime"] or m["mhcflurry"] or m["netctlpan"])
    if m["netctlpan"] and any(not (8 <= l <= 11) for l in p["lengths"]):
        notes = notes + ["NetCTLpan scores 8-11-mers only; other lengths get no NetCTLpan score"]
    out = {"n_peptides": n_pep, "n_windows": n_win, "n_alleles": len(p["alleles"]), "warnings": warns + notes,
           "n_pairs_netmhcpan": n_pep * len(nm_ok) if m["netmhcpan"] else 0,
           "n_pairs_fast": n_pep * len(p["alleles"]) if fast else 0,
           "n_pairs_bigmhc": n_pep * len(bm_ok) if use_bm else 0,
           "chained": bool(p["chain"] and use_bm and (m["netmhcpan"] or m["mixmhcpred"] or m["mhcflurry"]))}
    out["too_big"] = max(out["n_pairs_netmhcpan"], out["n_pairs_fast"], out["n_pairs_bigmhc"]) > MAX_PAIRS
    return out


def _binding_pass(df, f):
    """NetMHCpan / MixMHCpred / MHCflurry binding filters (the ones that gate BigMHC when chained)."""
    ok = pd.Series(True, index=df.index)
    if "NMP_%Rank_EL" in df:
        if f["rank_el_max"] is not None: ok &= df["NMP_%Rank_EL"] <= f["rank_el_max"]
        if f["rank_ba_max"] is not None and "NMP_%Rank_BA" in df: ok &= df["NMP_%Rank_BA"] <= f["rank_ba_max"]
        if f["aff_max"] is not None and "NMP_Aff_nM" in df: ok &= df["NMP_Aff_nM"] <= f["aff_max"]
    if "MixMHC_%Rank" in df and f["mix_rank_max"] is not None:
        ok &= df["MixMHC_%Rank"] <= f["mix_rank_max"]
    if "MHCF_Presentation_%Rank" in df:
        if f["mhcf_pres_rank_max"] is not None: ok &= df["MHCF_Presentation_%Rank"] <= f["mhcf_pres_rank_max"]
        if f["mhcf_aff_max"] is not None: ok &= df["MHCF_Affinity_nM"] <= f["mhcf_aff_max"]
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
    nm_ok, bm_ok, notes = _supported(p["alleles"], m)
    warns += notes
    uniq = peps["Peptide"].drop_duplicates().tolist()
    n_pep = len(uniq)
    if n_pep * len(p["alleles"]) > MAX_PAIRS:
        raise ValueError(f"{n_pep:,} peptides x {len(p['alleles'])} alleles is over the {MAX_PAIRS:,} limit; "
                         "use fewer alleles or lengths, or split the sequence.")
    job.log(f"{n_pep:,} unique peptides, {len(p['alleles'])} alleles")

    if m["netmhcpan"]:
        job.stage("netmhcpan", "NetMHCpan 4.2 (EL / BA)", weight=1.0)
    if m["mixmhcpred"]:
        job.stage("mixmhcpred", "MixMHCpred 3.0 (presentation)", weight=0.5)
    if m["prime"]:
        job.stage("prime", "PRIME 2.1 (immunogenicity)", weight=0.5)
    if m["iedb_imm"]:
        job.stage("iedb_imm", "IEDB immunogenicity (Calis 2013)", weight=0.2)
    if m["mhcflurry"]:
        job.stage("mhcflurry", "MHCflurry 2 (affinity, processing, presentation)", weight=1.0)
    if m["netctlpan"]:
        job.stage("netctlpan", "NetCTLpan 1.1 (processing)", weight=0.5)
    if m["netchop"]:
        job.stage("netchop", "NetChop 3.1 (proteasome cleavage)", weight=0.2)
    for x in bm_models:
        job.stage(f"bigmhc_{x}", f"BigMHC {x.upper()}", weight=2.0)
    job.stage("merge", "Merging and filtering", weight=0.3)

    manual = p["funnel"]["mode"] == "manual"
    seqs = nm.read_sequences(p["text"]) if p["mode"] == "scan" else None
    fstats = []
    nmp = mix = pri = iedb = mhcf = ctl = chop = None
    if manual:
        base, ran, fstats = funnel.run_manual(1, p, peps, seqs, job, warns)
        nmp, mix, pri, iedb, mhcf, ctl, chop = (base if k in ran else None for k in ("netmhcpan", "mixmhcpred", "prime", "iedb_imm", "mhcflurry", "netctlpan", "netchop"))
    else:
        # ---- 1. NetMHCpan ------------------------------------------------------------------------
        nmp = None
        if m["netmhcpan"]:
            if not nm_ok:
                raise ValueError("None of the selected alleles is known to NetMHCpan.")
            res, logs = nm.run_netmhcpan(uniq, [allele_lib.netmhcpan_name(a) for a in nm_ok],
                                         {**p["nmp"], "t": None, "raw": ""}, job=job)
            for l in logs[-3:]:
                job.log(l[:300])
            if res.empty:
                raise RuntimeError("NetMHCpan returned no results (see the log).")
            res["Allele"] = [allele_lib.canonical(a) for a in res["Allele"]]
            nmp = res
            job.finish_stage("netmhcpan", f"{len(res):,} predictions")

        # ---- 2. MixMHCpred and PRIME (fast, every pair) ----------------------------------------------
        mix = pri = None
        if m["mixmhcpred"]:
            mix, n = mx.predict(uniq, p["alleles"], job=job)
            warns += n
            job.finish_stage("mixmhcpred", f"{len(mix):,} predictions")
        if m["prime"]:
            pri, n = pr.predict(uniq, p["alleles"], job=job)
            warns += [x for x in n if x not in warns]
            job.finish_stage("prime", f"{len(pri):,} predictions")

        iedb = None
        if m["iedb_imm"]:
            iedb, n = ci.predict(uniq, p["alleles"], p["iedb"]["mask_choice"], p["iedb"]["positions"])
            warns += n
            job.finish_stage("iedb_imm", f"{len(iedb):,} scores")

        seqs = nm.read_sequences(p["text"]) if p["mode"] == "scan" else None
        mhcf = ctl = chop = None
        if m["mhcflurry"]:
            mhcf, n = mf.predict(peps, p["alleles"], sequences=seqs, job=job)
            warns += n
            job.finish_stage("mhcflurry", f"{len(mhcf):,} predictions")
        if m["netctlpan"]:
            ctl, n = ip.netctlpan(seqs, p["lengths"], p["alleles"], job=job)
            warns += n
            job.finish_stage("netctlpan", f"{len(ctl):,} scores")
        if m["netchop"]:
            chop = ip.netchop_cterm(seqs, peps, job=job)
            job.finish_stage("netchop", f"{len(chop):,} scores")
        if seqs and any(set(sq) - set("ACDEFGHIKLMNPQRSTVWY") for _, sq in seqs) and (m["netctlpan"] or m["netchop"]):
            warns.append("Characters that are not amino acids were replaced by X before NetCTLpan / NetChop; "
                         "scores near them are less reliable")

        # ---- 3. assemble the pair table so far (the chain gate needs it) ---------------------------------
        if nmp is not None:
            base = peps.merge(nmp, on="Peptide")
        else:
            base = peps.merge(pd.DataFrame({"Allele": p["alleles"]}), how="cross")
        for extra in (mix, pri, iedb):
            if extra is not None and not extra.empty:
                base = base.merge(extra, on=["Peptide", "Allele"], how="left")
        for extra in (mhcf, ctl):                                     # context-dependent: keyed by position
            if extra is not None and not extra.empty:
                base = base.merge(extra, on=["Protein", "Start", "Peptide", "Allele"], how="left")
        if chop is not None and not chop.empty:
            base = base.merge(chop, on=["Protein", "Start", "Peptide"], how="left")

        # ---- 4. BigMHC. Chained: EL scores the pairs that could still become binders, IM only the pairs that did ----------
        bmres = None
        if use_bm:
            gated = p["chain"] and (nmp is not None or mix is not None or mhcf is not None)
            rules = consensus.resolve_rules(1, p["preset"], p["rules"])
            ok_allele = base["Allele"].isin(bm_ok)
            if gated:
                first = _binding_pass(base, f) & consensus.gate_mask(base, 1, rules, pending=int(m["bigmhc_el"]))
            else:
                first = pd.Series(True, index=base.index)
            steps = [("el", first & ok_allele)] if m["bigmhc_el"] else []
            if m["bigmhc_im"]:
                steps.append(("im", None if m["bigmhc_el"] and gated else first & ok_allele))      # None: decided after EL
            for model, keep in steps:
                if keep is None:
                    keep = consensus.binder_mask(base, 1, rules) & _binding_pass(base, f) & first & ok_allele
                pairs = base.loc[keep, ["Peptide", "Allele"]].drop_duplicates()
                if gated:
                    job.log(f"chain: {len(pairs):,} pairs -> BigMHC {model.upper()}")
                if pairs.empty:
                    job.finish_stage(f"bigmhc_{model}", "nothing to score")
                    warns.append(f"BigMHC {model.upper()} had nothing to score (no pair passed the earlier tools, or no human alleles).")
                    continue
                sc = bm.predict_pairs(pairs.rename(columns={"Peptide": "peptide", "Allele": "allele"}),
                                      (model,), job=job).rename(columns={"peptide": "Peptide", "allele": "Allele"})
                base = base.merge(sc, on=["Peptide", "Allele"], how="left")

    # ---- 5. flags, ordering, per-peptide summary -------------------------------------------------------
    job.update("merge", msg="merging", state="running")
    base.insert(3, "Length", base["Peptide"].str.len())
    flags = []

    def flag(name, ok, scored_cols):
        # a pair with no score from this model cannot pass it (also when the model scored nothing at all)
        scored = base[scored_cols].notna().all(axis=1) if scored_cols else pd.Series(False, index=base.index)
        base[name] = ok & scored
        flags.append(name)

    yes = pd.Series(True, index=base.index)
    if nmp is not None:
        flag("Pass_NMP", _binding_pass(base[["NMP_%Rank_EL"] + [c for c in ("NMP_%Rank_BA", "NMP_Aff_nM") if c in base]], f),
             ["NMP_%Rank_EL"])
    if mix is not None:
        have = "MixMHC_%Rank" in base
        flag("Pass_Mix", (base["MixMHC_%Rank"] <= f["mix_rank_max"]) if (have and f["mix_rank_max"] is not None) else yes,
             ["MixMHC_%Rank"] if have else [])
    if pri is not None:
        have = "PRIME_%Rank" in base
        flag("Pass_PRIME", (base["PRIME_%Rank"] <= f["prime_rank_max"]) if (have and f["prime_rank_max"] is not None) else yes,
             ["PRIME_%Rank"] if have else [])
    if mhcf is not None:
        have = "MHCF_Presentation_%Rank" in base
        ok = _binding_pass(base[[c for c in ("MHCF_Presentation_%Rank", "MHCF_Affinity_nM") if c in base]], f) if have else yes
        flag("Pass_MHCF", ok, ["MHCF_Presentation_%Rank"] if have else [])
    if ctl is not None:
        have = "CTLpan_%Rank" in base
        flag("Pass_CTLpan", (base["CTLpan_%Rank"] <= f["ctlpan_rank_max"]) if (have and f["ctlpan_rank_max"] is not None) else yes,
             ["CTLpan_%Rank"] if have else [])
    if chop is not None:
        have = "NetChop_Cterm" in base
        flag("Pass_NetChop", (base["NetChop_Cterm"] >= f["netchop_min"]) if (have and f["netchop_min"] is not None) else yes,
             ["NetChop_Cterm"] if have else [])
    if iedb is not None:
        flag("Pass_IEDB", (base["IEDB_Immunogenicity"] >= f["iedb_imm_min"]) if f["iedb_imm_min"] is not None else yes,
             ["IEDB_Immunogenicity"])
    if use_bm:
        ok = yes.copy()
        if f["bm_el_min"] is not None and "BigMHC_EL" in base: ok &= base["BigMHC_EL"] >= f["bm_el_min"]
        if f["bm_im_min"] is not None and "BigMHC_IM" in base: ok &= base["BigMHC_IM"] >= f["bm_im_min"]
        flag("Pass_BigMHC", ok, [c for c in base if c.startswith("BigMHC_")])
    skipped = [k for k, models in (("rank_el_max", "netmhcpan"), ("rank_ba_max", "netmhcpan"), ("aff_max", "netmhcpan"),
                                   ("mix_rank_max", "mixmhcpred"), ("prime_rank_max", "prime"),
                                   ("iedb_imm_min", "iedb_imm"), ("mhcf_pres_rank_max", "mhcflurry"), ("mhcf_aff_max", "mhcflurry"),
                                   ("ctlpan_rank_max", "netctlpan"), ("netchop_min", "netchop")) if f[k] is not None and not m[models]]
    if skipped:
        warns.append("Filters ignored because their model was not run: " + ", ".join(skipped))
    base["Pass_All"] = base[flags].all(axis=1) if flags else pd.Series(True, index=base.index)
    if manual:
        base["Pass_All"] &= base["Funnel_Alive"]

    sort_cols, asc = ["Pass_All"], [False]
    for col, a in (("NMP_%Rank_EL", True), ("MixMHC_%Rank", True), ("MHCF_Presentation_%Rank", True),
                   ("PRIME_%Rank", True), ("CTLpan_%Rank", True), ("BigMHC_IM", False)):
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
                          ("MixMHC_%Rank", "min", "Best_MixMHC_%Rank"), ("PRIME_%Rank", "min", "Best_PRIME_%Rank"),
                          ("IEDB_Immunogenicity", "max", "Max_IEDB_Immunogenicity"),
                          ("MHCF_Presentation_%Rank", "min", "Best_MHCF_Presentation_%Rank"), ("CTLpan_%Rank", "min", "Best_CTLpan_%Rank"),
                          ("NetChop_Cterm", "max", "Max_NetChop_Cterm"), ("BigMHC_EL", "max", "Max_BigMHC_EL"), ("BigMHC_IM", "max", "Max_BigMHC_IM")):
        if col in base:
            summ = summ.merge(getattr(g[col], fn)().rename(name).reset_index(), on=keys, how="left")
    order, asc = ["Alleles_passed"], [False]
    for col, a_ in (("Best_NMP_%Rank_EL", True), ("Best_MixMHC_%Rank", True), ("Best_MHCF_Presentation_%Rank", True), ("Best_PRIME_%Rank", True), ("Max_BigMHC_IM", False)):
        if col in summ:
            order.append(col); asc.append(a_)
    summ = summ.sort_values(order, ascending=asc, na_position="last", kind="stable").reset_index(drop=True)

    job.update("merge", msg="consensus, population coverage, ranking", state="running")
    info = {"mode": p["mode"], "lengths": p["lengths"], "alleles": p["alleles"], "peptides": n_pep,
            "models": [k for k, on in m.items() if on], "filters": {k: v for k, v in f.items() if v is not None},
            "chain": p["chain"], "funnel": {**p["funnel"], "stats": fstats}, "input_chars": len(p["text"]),
            "input_peptides": peps["Peptide"].drop_duplicates().tolist() if p["mode"] == "peptides" else None}
    cons = report.finish(base, 1, p, job, info, warns, base["Pass_All"])

    jobdir = config.run_dir(job.id)
    pairs_csv, summ_csv = (os.path.join(jobdir, f"{job.id}_{n}.csv") for n in ("pairs", "summary"))
    base.round(6).to_csv(pairs_csv, index=False)
    summ.round(6).to_csv(summ_csv, index=False)
    job.finish_stage("merge", "done")

    job.result = {**{k: v for k, v in cons.items()}, "pairs_csv": pairs_csv, "summary_csv": summ_csv, "summary": {
        "peptides": n_pep, "pairs": int(len(base)), "pairs_passed": int(base["Pass_All"].sum()),
        "peptides_passed": int((summ["Alleles_passed"] > 0).sum()),
        "seconds": round(time.time() - t0, 1), "warnings": list(dict.fromkeys(warns)),
        "chained": bool(use_bm and p["chain"] and (nmp is not None or mix is not None or mhcf is not None))}}
    job.log(f"done in {time.time() - t0:.1f}s: {int(base['Pass_All'].sum()):,} pairs / "
            f"{int((summ['Alleles_passed'] > 0).sum()):,} peptides pass")
