"""The funnel: in which order the models run, and which peptides go on to the next one.

  mode "auto"    (default) quick tools score everything; BigMHC EL scores only pairs that can still become binders,
                 BigMHC IM only binders, CD4episcore only peptides that bind (exact for the agreement rule, see consensus.gate_mask)
  mode "all"     no funnel: every model scores every peptide-allele pair, filtering happens afterwards
  mode "manual"  you define stages. Each stage runs its models on the pairs still alive, then a pass rule decides which pairs go on:
        all        every pair goes on
        agree      pairs that the stage's tools call a hit (any / majority / all of them, optional cut-off)
        top_n      the best N peptides (or pairs) by the stage's first model
        top_pct    the best P % of peptides (or pairs)
  A stage's models see the same input, so they can run side by side (performance.parallel_models).

The output table gets, for every pair: S<k>_In (was it scored in stage k), S<k>_Pass (did it go on), Funnel_Reached (stages that
scored it), Stopped_At (the stage that stopped it, 0 = it came through) and Funnel_Alive.
"""
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

import consensus
import models
from common import config
from common.jobs import THREADS

RULES = ("all", "agree", "top_n", "top_pct")
KEYS = ["Protein", "Start", "Peptide"]
MODEL_TOOLS = {"netmhcpan": ("NMP_EL", "NMP_BA"), "mixmhcpred": ("MIXMHC",), "mhcflurry": ("MHCF",), "prime": ("PRIME",),
               "iedb_imm": ("IEDB_IMM",), "netctlpan": ("CTLPAN",), "netchop": ("NETCHOP",), "bigmhc_el": ("BIGMHC_EL",),
               "bigmhc_im": ("BIGMHC_IM",), "netmhciipan_el": ("NII_EL",), "netmhciipan_ba": ("NII_BA",),
               "mixmhc2pred": ("MIX2",), "mhcnuggets": ("MHCNUG",), "cd4episcore": ("CD4",)}


def templates(cls):
    """Starting points for the builder: cheap tools first, expensive ones on the survivors."""
    if cls == 1:
        return {
            "recommended": ("Cheap first, slow tools on the survivors",
                            [{"models": ["mixmhcpred", "netmhcpan", "mhcflurry"], "pass": {"rule": "agree", "need": "any", "cutoff": 5}},
                             {"models": ["prime", "iedb_imm", "netctlpan", "netchop"], "pass": {"rule": "all"}},
                             {"models": ["bigmhc_el"], "pass": {"rule": "agree", "need": "any"}},
                             {"models": ["bigmhc_im"], "pass": {"rule": "all"}}]),
            "top": ("Best 2,000 peptides by MixMHCpred, then everything else",
                    [{"models": ["mixmhcpred"], "pass": {"rule": "top_n", "n": 2000, "unit": "peptide"}},
                     {"models": ["netmhcpan", "mhcflurry", "prime", "iedb_imm", "bigmhc_el", "bigmhc_im"], "pass": {"rule": "all"}}]),
            "two": ("One fast screen, then all the others",
                    [{"models": ["mixmhcpred"], "pass": {"rule": "agree", "need": "any", "cutoff": 3}},
                     {"models": ["netmhcpan", "mhcflurry", "prime", "iedb_imm", "netctlpan", "netchop", "bigmhc_el", "bigmhc_im"], "pass": {"rule": "all"}}])}
    return {
        "recommended": ("Cheap first, slow tools on the survivors",
                        [{"models": ["mixmhc2pred"], "pass": {"rule": "agree", "need": "any", "cutoff": 10}},
                         {"models": ["netmhciipan_el", "netmhciipan_ba"], "pass": {"rule": "agree", "need": "any"}},
                         {"models": ["cd4episcore"], "pass": {"rule": "all"}}]),
        "two": ("NetMHCIIpan screen, then the rest",
                [{"models": ["netmhciipan_el"], "pass": {"rule": "agree", "need": "any", "cutoff": 10}},
                 {"models": ["netmhciipan_ba", "mixmhc2pred", "cd4episcore"], "pass": {"rule": "all"}}])}


def normalise(spec, cls, mode):
    """-> {"mode", "stages"} validated; raises ValueError with a message for the user. `mode` = scan | peptides."""
    spec = spec or {}
    fm = spec.get("mode", "auto")
    if fm not in ("auto", "all", "manual"):
        raise ValueError("Funnel mode must be auto, all or manual.")
    if fm != "manual":
        return {"mode": fm, "stages": []}
    stages, seen = [], set()
    for i, st in enumerate(spec.get("stages") or [], 1):
        ms = [m for m in (st.get("models") or [])]
        if not ms:
            raise ValueError(f"Funnel stage {i} has no model.")
        for m in ms:
            if m not in models.SPECS or models.SPECS[m].cls != cls:
                raise ValueError(f"Funnel stage {i}: '{m}' is not a class {'I' if cls == 1 else 'II'} model.")
            if m in seen:
                raise ValueError(f"'{models.SPECS[m].label}' appears in more than one stage; a model runs once.")
            if models.SPECS[m].protein and mode == "peptides":
                raise ValueError(f"{models.SPECS[m].label} needs the protein: use 'Cut a protein' or remove it from the funnel.")
            seen.add(m)
        pr = dict(st.get("pass") or {"rule": "all"})
        rule = pr.get("rule", "all")
        if rule not in RULES:
            raise ValueError(f"Funnel stage {i}: pass rule must be one of {', '.join(RULES)}.")
        if rule == "top_n" and not (pr.get("n") and int(pr["n"]) >= 1):
            raise ValueError(f"Funnel stage {i}: top_n needs n >= 1.")
        if rule == "top_pct" and not (pr.get("pct") and 0 < float(pr["pct"]) <= 100):
            raise ValueError(f"Funnel stage {i}: top_pct needs 0 < pct <= 100.")
        if rule == "agree" and pr.get("need", "any") not in ("any", "majority", "all") and not str(pr.get("need")).isdigit():
            raise ValueError(f"Funnel stage {i}: need must be any, majority, all or a number.")
        if rule in ("top_n", "top_pct") and not models.SPECS[ms[0]].main:
            raise ValueError(f"Funnel stage {i}: cannot rank by {ms[0]}.")
        stages.append({"models": ms, "pass": pr})
    if not stages:
        raise ValueError("A manual funnel needs at least one stage with a model.")
    return {"mode": "manual", "stages": stages}


def _votes(df, alive, keys, cls, rules, pr):
    """(votes, scored) per row over the stage's tools, among the rows alive."""
    tools = [t for t in consensus.TOOLS if t.cls == cls and t.col in df.columns and any(t.key in MODEL_TOOLS[k] for k in keys)]
    votes = pd.Series(0, index=df.index)
    scored = pd.Series(0, index=df.index)
    cutoff = pr.get("cutoff")
    for t in tools:
        c = df[t.col].where(alive)
        if t.kind == "abs":
            thr = cutoff if cutoff is not None else {"binding": rules["bind_rank_max"], "immunogenicity": rules["imm_abs_rank_max"], "processing": rules["ctlpan_rank_max"]}[t.group]
            yes = c <= float(thr)
        elif t.key == "BIGMHC_EL":
            yes = c >= float(cutoff if cutoff is not None else rules["bigmhc_el_min"])
        elif t.key == "NETCHOP":
            yes = c >= float(cutoff if cutoff is not None else (rules["netchop_min"] if rules["netchop_min"] is not None else 0.5))
        else:                                                   # relative tools: top q % of the rows scored in this stage
            q = float(cutoff) if cutoff is not None else rules["imm_top_q"] * 100
            asc = t.direction == "low"
            rk = c.groupby(df["Allele"]).rank(pct=True, ascending=asc) if t.allele_dependent and "Allele" in df else c.rank(pct=True, ascending=asc)
            yes = rk * 100 <= q
        votes += yes.fillna(False).astype(int)
        scored += c.notna().astype(int)
    return votes, scored


def _passes(df, alive, stage, cls, rules):
    pr = stage["pass"]
    rule = pr["rule"]
    if rule == "all":
        return alive.copy()
    if rule == "agree":
        votes, scored = _votes(df, alive, stage["models"], cls, rules, pr)
        need = consensus._need(pr.get("need", "any") if not str(pr.get("need", "")).isdigit() else int(pr["need"]), scored)
        return alive & (scored > 0) & (votes >= need)
    spec = models.SPECS[stage["models"][0]]
    col = spec.main
    s = df[col].where(alive)
    if pr.get("unit", "peptide") == "peptide":
        key = df["Peptide"]
        best = (s.groupby(key).min() if spec.asc else s.groupby(key).max()).dropna().sort_values(ascending=spec.asc, kind="stable")
        n = int(pr["n"]) if rule == "top_n" else max(1, int(np.ceil(len(best) * float(pr["pct"]) / 100)))
        keep = set(best.index[:n])
        return alive & df["Peptide"].isin(keep) & s.notna()
    order = s.dropna().sort_values(ascending=spec.asc, kind="stable")
    n = int(pr["n"]) if rule == "top_n" else max(1, int(np.ceil(len(order) * float(pr["pct"]) / 100)))
    keep = pd.Series(False, index=df.index)
    keep.loc[order.index[:n]] = True
    return alive & keep


def run_manual(cls, p, peps, seqs, job, warns):
    """Runs the stages. -> (base, ran: set of model keys, stats list). `base` has every pair, scores where a model ran, and the flags."""
    stages = p["funnel"]["stages"]
    rules = consensus.resolve_rules(cls, p["preset"], p["rules"])
    base = peps.merge(pd.DataFrame({"Allele": p["alleles"]}), how="cross")
    alive = pd.Series(True, index=base.index)
    reached = pd.Series(0, index=base.index)
    stopped = pd.Series(0, index=base.index)
    ctx = {"seqs": seqs, "p": p}
    perf = config.performance()
    ran, stats = set(), []
    for k, st in enumerate(stages, 1):
        t0 = time.time()
        rows = base.loc[alive, KEYS + ["Allele"]]
        n_pairs, n_peps = len(rows), rows["Peptide"].nunique()
        job.log(f"funnel stage {k}: {', '.join(st['models'])} on {n_peps:,} peptides / {n_pairs:,} pairs")
        base[f"S{k}_In"] = alive
        reached += alive.astype(int)
        if n_pairs == 0:
            for m in st["models"]:
                job.finish_stage(m if m != "netmhciipan_el" and m != "netmhciipan_ba" else "netmhciipan", "nothing left to score")
            stats.append({"stage": k, "models": st["models"], "pairs_in": 0, "peptides_in": 0, "pairs_out": 0, "peptides_out": 0, "seconds": 0})
            base[f"S{k}_Pass"] = False
            continue
        groups = [[m] for m in st["models"] if not m.startswith("netmhciipan")]
        nn_keys = [m for m in st["models"] if m.startswith("netmhciipan")]
        if nn_keys:
            groups.append(nn_keys)
        workers = min(perf["parallel_models"], len(groups))
        share = max(1, (perf["threads"] or perf["cpus"]) // workers) if workers > 1 else None

        def work(g):
            if share:
                THREADS.set(share)
            return models.run(g, rows, ctx, job)
        if workers > 1:
            with ThreadPoolExecutor(max_workers=workers) as ex:
                outs = list(ex.map(work, groups))
        else:
            outs = [work(g) for g in groups]
        for g, (frames, notes) in zip(groups, outs):
            warns += [n for n in notes if n not in warns]
            for df in frames:
                if df is None or df.empty:
                    continue
                mk = models.merge_keys(df)
                df = df.drop_duplicates(mk).merge(rows[mk].drop_duplicates(), on=mk, how="inner")        # keep the pairs this stage asked for
                base = base.merge(df, on=mk, how="left")
            for m in g:
                job.finish_stage("netmhciipan" if m.startswith("netmhciipan") else m, f"{n_pairs:,} pairs")
            ran.update(g)
        new_alive = _passes(base, alive, st, cls, rules)
        stopped[alive & ~new_alive] = k
        base[f"S{k}_Pass"] = new_alive
        alive = new_alive
        stats.append({"stage": k, "models": st["models"], "rule": st["pass"], "pairs_in": n_pairs, "peptides_in": int(n_peps),
                      "pairs_out": int(alive.sum()), "peptides_out": int(base.loc[alive, "Peptide"].nunique()), "seconds": round(time.time() - t0, 1)})
    base["Funnel_Reached"] = reached
    base["Stopped_At"] = stopped
    base["Funnel_Alive"] = alive
    return base, ran, stats
