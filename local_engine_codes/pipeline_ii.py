"""MHC class II: one form in, one table out.

    NetMHCIIpan (EL + BA)  ->  MixMHC2pred  ->  CD4episcore  ->  TLimmuno2  ->  ImmuScope

NetMHCIIpan and MixMHC2pred score every peptide-allele pair. CD4episcore scores the peptide on its own (no HLA
typing; lower percentile = more immunogenic), so with `chain` on it only runs for peptides that pass the binding
filters for at least one allele. TLimmuno2 and ImmuScope score peptide-allele pairs (higher = more immunogenic); with
`chain` on they only score the pairs already called binders.

Import after engine_server has put the tool folders on sys.path.
"""
import os, time

import pandas as pd

import consensus
import funnel
import report
from guard import try_model

from common import alleles as allele_lib, config
import cd4episcore_runner as cd
import immuscope_runner as isc
import mixmhc2pred_runner as m2
import netmhcpan_app as nm
import netmhciipan_runner as nn
import tlimmuno2_runner as tl2

MAX_PAIRS = 3_000_000
SCAN_LENGTHS = tuple(range(11, 26))
PEPTIDE_BOUNDS = (9, 30)
MODEL_KEYS = ("netmhciipan_el", "netmhciipan_ba", "mixmhc2pred", "cd4episcore", "tlimmuno2", "immuscope")
PAIR_IMM = (("tlimmuno2", "TLimmuno2", tl2, "TLimmuno2"), ("immuscope", "ImmuScope-IM", isc, "ImmuScope_IM"))
FILTER_KEYS = ("nii_el_rank_max", "nii_ba_rank_max", "nii_ic50_max", "mix2_rank_max", "cd4_max")


def _num(x):
    return None if x in (None, "") else float(x)


def normalise(params):
    p = dict(params or {})
    fun = funnel.normalise(p.get("funnel"), 2, p.get("mode", "scan"))
    models = {k: False for k in MODEL_KEYS}
    if fun["mode"] == "manual":
        models.update({m: True for st in fun["stages"] for m in st["models"]})
    else:
        models.update({k: bool(v) for k, v in (p.get("models") or {}).items() if k in models})
    if not any(models.values()):
        raise ValueError("Select at least one model.")
    raw = [str(a) for a in p.get("alleles", []) if str(a).strip()]
    alleles = list(dict.fromkeys(allele_lib.canonical_ii(a) for a in raw))
    if not alleles:
        raise ValueError("Select at least one allele.")
    wrong = [a for a in alleles if not allele_lib.is_class_ii(a)]
    if wrong:
        raise ValueError("Not MHC class II alleles: " + ", ".join(wrong) + " (DR, DQ, DP or mouse I-A / I-E expected).")
    lengths = sorted({int(x) for x in p.get("lengths", [15]) if int(x) in SCAN_LENGTHS})
    if p.get("mode", "scan") == "scan" and not lengths:
        raise ValueError("Select at least one peptide length (11-25).")
    version = (p.get("nii") or {}).get("version", "4.3")
    if version not in nn.VERSIONS:
        raise ValueError(f"NetMHCIIpan version must be one of {nn.VERSIONS}.")
    f = {k: _num((p.get("filters") or {}).get(k)) for k in FILTER_KEYS}
    return {"text": p.get("text", ""), "mode": p.get("mode", "scan"), "lengths": lengths, "alleles": alleles,
            "models": models, "chain": bool(p.get("chain", True)) and fun["mode"] == "auto", "funnel": fun, "filters": f, "version": version,
            "preset": p.get("preset") or "balanced", "rules": p.get("rules") or {}}


def preview(params):
    p = normalise(params)
    m = p["models"]
    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"], PEPTIDE_BOUNDS)
    n_pep = int(peps["Peptide"].nunique()) if not peps.empty else 0
    n_win = int(len(peps))
    notes = []
    if m["mixmhc2pred"]:
        try:
            _, notes = m2.resolve(p["alleles"], config.model("mixmhc2pred")["home"]) if config.model("mixmhc2pred").get("home") else ({}, [])
        except OSError:
            pass
    nii = m["netmhciipan_el"] or m["netmhciipan_ba"]
    out = {"n_peptides": n_pep, "n_windows": n_win, "n_alleles": len(p["alleles"]), "warnings": warns + notes,
           "n_pairs_netmhcpan": n_pep * len(p["alleles"]) if nii else 0,
           "n_pairs_fast": n_pep * len(p["alleles"]) if m["mixmhc2pred"] else 0, "n_pairs_bigmhc": 0,
           "chained": bool(p["chain"] and (m["cd4episcore"] or m["tlimmuno2"] or m["immuscope"]) and (nii or m["mixmhc2pred"]))}
    out["too_big"] = max(out["n_pairs_netmhcpan"], out["n_pairs_fast"]) > MAX_PAIRS
    return out


def _binding_pass(df, f):
    ok = pd.Series(True, index=df.index)
    for col, key in (("NII_%Rank_EL", "nii_el_rank_max"), ("NII_%Rank_BA", "nii_ba_rank_max"),
                     ("NII_IC50_nM", "nii_ic50_max"), ("MIX2_%Rank", "mix2_rank_max")):
        if col in df and f[key] is not None:
            ok &= df[col] <= f[key]
    return ok


def n_unique(p):
    peps, _ = nm.build_peptides(p["text"], p["mode"], p["lengths"], PEPTIDE_BOUNDS)
    return int(peps["Peptide"].nunique()) if not peps.empty else 0


def stage_plan(p, n_pep):
    from common import timing
    m, pairs = p["models"], n_pep * len(p["alleles"])
    kinds = [k for k, on in (("EL", m["netmhciipan_el"]), ("BA", m["netmhciipan_ba"])) if on]
    gate = timing.GATE * 3 if p["chain"] else 1.0
    plan = []
    if kinds:
        plan.append({"name": "netmhciipan", "label": f"NetMHCIIpan {p['version']} ({' + '.join(kinds)})", "work": pairs * len(kinds),
                     "weight": timing.expected_seconds("netmhciipan", pairs * len(kinds))})
    if m["mixmhc2pred"]:
        plan.append({"name": "mixmhc2pred", "label": "MixMHC2pred 2.1", "work": pairs, "weight": timing.expected_seconds("mixmhc2pred", pairs)})
    if m["cd4episcore"]:
        plan.append({"name": "cd4episcore", "label": "CD4episcore (CD4 immunogenicity)", "work": int(n_pep * gate), "weight": timing.expected_seconds("cd4episcore", n_pep * gate)})
    for key, label, _, _ in PAIR_IMM:
        if m[key]:
            plan.append({"name": key, "label": f"{label} (CD4 immunogenicity)", "work": int(pairs * gate), "weight": timing.expected_seconds(key, pairs * gate)})
    plan.append({"name": "merge", "label": "Merging and filtering", "work": pairs, "weight": timing.expected_seconds("merge", pairs)})
    return plan


def run(params, job):
    t0 = time.time()
    p = normalise(params)
    m, f = p["models"], p["filters"]
    kinds = tuple(k for k, on in (("el", m["netmhciipan_el"]), ("ba", m["netmhciipan_ba"])) if on)
    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"], PEPTIDE_BOUNDS)
    if peps.empty:
        raise ValueError("No valid peptides in the input.")
    peps = peps.drop_duplicates()
    seqs = nm.read_sequences(p["text"]) if p["mode"] == "scan" else None
    uniq = peps["Peptide"].drop_duplicates().tolist()
    n_pep = len(uniq)
    if n_pep * len(p["alleles"]) > MAX_PAIRS:
        raise ValueError(f"{n_pep:,} peptides x {len(p['alleles'])} alleles is over the {MAX_PAIRS:,} limit; "
                         "use fewer alleles or lengths, or split the sequence.")
    job.log(f"{n_pep:,} unique peptides, {len(p['alleles'])} alleles (MHC class II)")
    for st in stage_plan(p, n_pep):
        job.stage(st["name"], st["label"], st["weight"], work=st["work"])
    scan_keys = ["Protein", "Start", "Peptide", "Allele"]

    manual = p["funnel"]["mode"] == "manual"
    failed = []                                  # models that broke and were skipped (see guard.py)
    n_selected = (int(m["netmhciipan_el"] or m["netmhciipan_ba"]) + int(m["mixmhc2pred"]) + int(m["cd4episcore"]) + int(m["tlimmuno2"]) + int(m["immuscope"]))
    fstats = []
    nii = mix = chop = None
    pair_ran = []                                # pair-wise immunogenicity tools that produced scores
    if manual:
        base, ran, fstats = funnel.run_manual(2, p, peps, seqs, job, warns)
        nii = base if ran & {"netmhciipan_el", "netmhciipan_ba"} else None
        mix = base if "mixmhc2pred" in ran else None
        chop = base if "cd4episcore" in ran else None
        pair_ran = [(k, c) for k, _, _, c in PAIR_IMM if k in ran]
    else:
        nii = None
        if kinds:
            def _nii():
                df, n = nn.predict(peps, p["alleles"], p["lengths"] or sorted({len(x) for x in uniq}), seqs, kinds, p["version"], job=job)
                if df.empty:
                    raise RuntimeError("NetMHCIIpan returned no results (see the log).")
                return df, n
            got = try_model("NetMHCIIpan", "netmhciipan", _nii, job, warns, failed)
            if got is None:
                m["netmhciipan_el"] = m["netmhciipan_ba"] = False
                kinds = ()
            else:
                nii, n = got
                warns += n
                job.finish_stage("netmhciipan", f"{len(nii):,} predictions")
        mix = None
        if m["mixmhc2pred"]:
            got = try_model("MixMHC2pred", "mixmhc2pred", lambda: m2.predict(peps, p["alleles"], seqs, job=job), job, warns, failed)
            if got is None:
                m["mixmhc2pred"] = False
            else:
                mix, n = got
                warns += n
                job.finish_stage("mixmhc2pred", f"{len(mix):,} predictions")

        if nii is not None:
            base = peps.merge(nii, on=["Protein", "Start", "Peptide"] if seqs else ["Peptide"])
        else:
            base = peps.merge(pd.DataFrame({"Allele": p["alleles"]}), how="cross")
        if mix is not None and not mix.empty:
            base = base.merge(mix, on=scan_keys, how="left")


        chop = None
        if m["cd4episcore"]:
            gated = p["chain"] and (nii is not None or mix is not None)
            if gated:
                rules = consensus.resolve_rules(2, p["preset"], p["rules"])
                use = base[_binding_pass(base, f) & consensus.gate_mask(base, 2, rules)]["Peptide"].drop_duplicates().tolist()
            else:
                use = uniq
            if gated:
                job.log(f"chain: {len(use):,} of {n_pep:,} peptides pass the binding filters -> CD4episcore")
            if use:
                chop = try_model("CD4episcore", "cd4episcore", lambda: cd.predict(use, job=job), job, warns, failed)
                if chop is None:
                    m["cd4episcore"] = False
                else:
                    base = base.merge(chop, on="Peptide", how="left")
                    job.finish_stage("cd4episcore", f"{len(chop):,} peptides scored")
            else:
                job.finish_stage("cd4episcore", "nothing to score")
                warns.append("CD4episcore had nothing to score (no peptide passed the binding filters).")

        for key, label, runner, col in PAIR_IMM:
            if not m[key]:
                continue
            gated = p["chain"] and (nii is not None or mix is not None)
            if gated:
                rules = consensus.resolve_rules(2, p["preset"], p["rules"])
                keep = consensus.binder_mask(base, 2, rules) & _binding_pass(base, f) & consensus.gate_mask(base, 2, rules)
                pairs = base.loc[keep, ["Peptide", "Allele"]].drop_duplicates()
                job.log(f"chain: {len(pairs):,} pairs -> {label}")
            else:
                pairs = base[["Peptide", "Allele"]].drop_duplicates()
            if pairs.empty:
                job.finish_stage(key, "nothing to score")
                warns.append(f"{label} had nothing to score (no pair passed the binding filters).")
                m[key] = False
                continue
            got = try_model(label, key, lambda: runner.predict(pairs, job=job), job, warns, failed)
            if got is None:
                m[key] = False
                continue
            sc, n = got
            warns += n
            base = base.merge(sc, on=["Peptide", "Allele"], how="left")
            pair_ran.append((key, col))
            job.finish_stage(key, f"{int(sc[col].notna().sum()):,} pairs scored")

    if failed and len(failed) >= n_selected:
        raise RuntimeError("Every selected model failed: " + ", ".join(failed) + ". See the log for each reason.")
    job.update("merge", msg="merging", state="running")
    base.insert(3, "Length", base["Peptide"].str.len())
    flags = []
    yes = pd.Series(True, index=base.index)

    def flag(name, ok, scored_cols):
        scored = base[scored_cols].notna().all(axis=1) if scored_cols else pd.Series(False, index=base.index)
        base[name] = ok & scored
        flags.append(name)

    if nii is not None:
        cols = [c for c in (("NII_%Rank_EL",) if "el" in kinds else ()) + (("NII_%Rank_BA",) if "ba" in kinds else ()) if c in base]
        flag("Pass_NII", _binding_pass(base[[c for c in ("NII_%Rank_EL", "NII_%Rank_BA", "NII_IC50_nM") if c in base]], f), cols)
    if mix is not None:
        have = "MIX2_%Rank" in base
        flag("Pass_Mix2", (base["MIX2_%Rank"] <= f["mix2_rank_max"]) if (have and f["mix2_rank_max"] is not None) else yes,
             ["MIX2_%Rank"] if have else [])
    if m["cd4episcore"]:
        have = "CD4_Episcore" in base
        flag("Pass_CD4", (base["CD4_Episcore"] <= f["cd4_max"]) if (have and f["cd4_max"] is not None) else yes,
             ["CD4_Episcore"] if have else [])
    skipped = [k for k, model in (("nii_el_rank_max", "netmhciipan_el"), ("nii_ba_rank_max", "netmhciipan_ba"),
                                  ("nii_ic50_max", "netmhciipan_ba"), ("mix2_rank_max", "mixmhc2pred"),
                                  ("cd4_max", "cd4episcore")) if f[k] is not None and not m[model]]
    if skipped:
        warns.append("Filters ignored because their model was not run: " + ", ".join(skipped))
    base["Pass_All"] = base[flags].all(axis=1)
    if manual:
        base["Pass_All"] &= base["Funnel_Alive"]

    sort_cols, asc = ["Pass_All"], [False]
    for col, a in (("NII_%Rank_EL", True), ("NII_%Rank_BA", True), ("MIX2_%Rank", True), ("CD4_Episcore", True)):
        if col in base:
            sort_cols.append(col); asc.append(a)
    base = base.sort_values(sort_cols, ascending=asc, na_position="last", kind="stable").reset_index(drop=True)

    keys = ["Protein", "Start", "Peptide"]
    g = base.groupby(keys, sort=False)
    summ = g.agg(Length=("Length", "first"), Alleles_tested=("Allele", "nunique"), Alleles_passed=("Pass_All", "sum")).reset_index()
    passed = base[base["Pass_All"]].groupby(keys, sort=False)["Allele"].agg(lambda a: ", ".join(sorted(a)))
    summ = summ.merge(passed.rename("Alleles_passed_list").reset_index(), on=keys, how="left")
    summ["Alleles_passed_list"] = summ["Alleles_passed_list"].fillna("")
    for col, fn, name in (("NII_%Rank_EL", "min", "Best_NII_%Rank_EL"), ("NII_%Rank_BA", "min", "Best_NII_%Rank_BA"),
                          ("NII_IC50_nM", "min", "Best_NII_IC50_nM"), ("MIX2_%Rank", "min", "Best_MIX2_%Rank"),
                          ("CD4_Episcore", "min", "CD4_Episcore"), ("ImmuScope_IM", "max", "Max_ImmuScope_IM"), ("TLimmuno2", "max", "Max_TLimmuno2")):
        if col in base:
            summ = summ.merge(getattr(g[col], fn)().rename(name).reset_index(), on=keys, how="left")
    order, asc = ["Alleles_passed"], [False]
    for col in ("Best_NII_%Rank_EL", "Best_NII_%Rank_BA", "Best_MIX2_%Rank", "CD4_Episcore"):
        if col in summ:
            order.append(col); asc.append(True)
    summ = summ.sort_values(order, ascending=asc, na_position="last", kind="stable").reset_index(drop=True)

    job.update("merge", msg="consensus, population coverage, ranking", state="running")
    info = {"mode": p["mode"], "lengths": p["lengths"], "alleles": p["alleles"], "peptides": n_pep,
            "models": [k for k, on in m.items() if on], "filters": {k: v for k, v in f.items() if v is not None},
            "chain": p["chain"], "funnel": {**p["funnel"], "stats": fstats}, "netmhciipan_version": p["version"], "input_chars": len(p["text"])}
    cons = report.finish(base, 2, p, job, info, warns, base["Pass_All"])

    jobdir = config.run_dir(job.id)
    pairs_csv, summ_csv = (os.path.join(jobdir, f"{job.id}_{n}.csv") for n in ("pairs", "summary"))
    base.round(6).to_csv(pairs_csv, index=False)
    summ.round(6).to_csv(summ_csv, index=False)
    job.finish_stage("merge", "done")
    job.result = {**cons, "pairs_csv": pairs_csv, "summary_csv": summ_csv, "summary": {
        "peptides": n_pep, "pairs": int(len(base)), "pairs_passed": int(base["Pass_All"].sum()),
        "peptides_passed": int((summ["Alleles_passed"] > 0).sum()), "seconds": round(time.time() - t0, 1),
        "warnings": list(dict.fromkeys(warns)), "failed_models": failed, "chained": bool((m["cd4episcore"] or m["tlimmuno2"] or m["immuscope"]) and p["chain"] and (nii is not None or mix is not None))}}
    job.log(f"done in {time.time() - t0:.1f}s: {int(base['Pass_All'].sum()):,} pairs / "
            f"{int((summ['Alleles_passed'] > 0).sum()):,} peptides pass")
