"""Run ONE model and keep its raw output, to check a tool on its own or compare it with the combined pipeline.

The table has every column the tool reports, with the position of each peptide in the protein, sorted by the tool's own
main score. Nothing is filtered or combined; the same code paths and settings as the pipeline run the tool, so a number
here is the number the pipeline used.

Import after engine_server has put the tool folders on sys.path.
"""
import time

import pandas as pd

import pipeline
import pipeline_ii
from common import alleles as allele_lib, config
import bigmhc_runner as bm
import calis_immunogenicity as ci
import cd4episcore_runner as cd
import iedb_processing as ip
import mhcflurry_runner as mf
import mhcnuggets_runner as mn
import mixmhc2pred_runner as m2
import mixmhcpred_runner as mx
import netmhcpan_app as nm
import netmhciipan_runner as nn
import prime_runner as pr

# key -> (class, label, needs alleles, needs the protein, main score column, ascending?)
MODELS = {
    "netmhcpan": (1, "NetMHCpan 4.2", True, False, "NMP_%Rank_EL", True),
    "mixmhcpred": (1, "MixMHCpred 3.0", True, False, "MixMHC_%Rank", True),
    "mhcflurry": (1, "MHCflurry 2", True, False, "MHCF_Presentation_%Rank", True),
    "prime": (1, "PRIME 2.1", True, False, "PRIME_%Rank", True),
    "iedb_imm": (1, "IEDB immunogenicity (Calis)", True, False, "IEDB_Immunogenicity", False),
    "netctlpan": (1, "NetCTLpan 1.1", True, True, "CTLpan_%Rank", True),
    "netchop": (1, "NetChop 3.1", False, True, "NetChop_Cterm", False),
    "bigmhc_el": (1, "BigMHC EL", True, False, "BigMHC_EL", False),
    "bigmhc_im": (1, "BigMHC IM", True, False, "BigMHC_IM", False),
    "netmhciipan": (2, "NetMHCIIpan (EL + BA)", True, False, "NII_%Rank_EL", True),
    "mixmhc2pred": (2, "MixMHC2pred 2.1", True, False, "MIX2_%Rank", True),
    "mhcnuggets": (2, "MHCnuggets", True, False, "MHCnuggets_%Rank", True),
    "cd4episcore": (2, "CD4episcore", False, False, "CD4_Episcore", True),
}


def normalise(params):
    p = dict(params or {})
    model = p.get("model")
    if model not in MODELS:
        raise ValueError("Choose a model: " + ", ".join(MODELS))
    cls, label, needs_alleles, needs_prot, _, _ = MODELS[model]
    mode = p.get("mode", "scan")
    if needs_prot and mode == "peptides":
        raise ValueError(f"{label} needs the protein the peptides come from: choose 'Cut a protein / FASTA into fragments'.")
    lengths = sorted({int(x) for x in p.get("lengths", [9 if cls == 1 else 15])})
    lo, hi = (8, 14) if cls == 1 else (11, 25)
    lengths = [x for x in lengths if lo <= x <= hi]
    if mode == "scan" and not lengths:
        raise ValueError(f"Select at least one peptide length ({lo}-{hi}).")
    canon = allele_lib.canonical if cls == 1 else allele_lib.canonical_ii
    alleles = list(dict.fromkeys(canon(a) for a in p.get("alleles", []) if str(a).strip()))
    if needs_alleles and not alleles:
        raise ValueError("Select at least one allele.")
    if cls == 2 and needs_alleles:
        wrong = [a for a in alleles if not allele_lib.is_class_ii(a)]
        if wrong:
            raise ValueError("Not MHC class II alleles: " + ", ".join(wrong))
    iedb = {"mask_choice": (p.get("iedb") or {}).get("mask_choice", "default"), "positions": str((p.get("iedb") or {}).get("positions", ""))}
    if model == "iedb_imm" and iedb["mask_choice"] not in ci.MASK_CHOICES:
        raise ValueError("Positions to mask must be default, by_allele or custom.")
    if model == "iedb_imm" and iedb["mask_choice"] == "custom":
        ci.parse_positions(iedb["positions"])
    version = (p.get("nii") or {}).get("version", "4.3")
    if model == "netmhciipan" and version not in nn.VERSIONS:
        raise ValueError(f"NetMHCIIpan version must be one of {nn.VERSIONS}.")
    nmp = {"ba": True, "rankS": 0.5, "rankW": 2.0}
    return {"model": model, "cls": cls, "text": p.get("text", ""), "mode": mode, "lengths": lengths, "alleles": alleles,
            "iedb": iedb, "version": version, "nmp": nmp}


def needed(p):
    """Health-check names the model depends on."""
    m = p["model"]
    return {"netmhcpan": ["netmhcpan"], "mixmhcpred": ["mixmhcpred"], "mhcflurry": ["mhcflurry"], "prime": ["mixmhcpred", "prime"],
            "iedb_imm": [], "netctlpan": ["iedb_bundle"], "netchop": ["iedb_bundle"], "bigmhc_el": ["bigmhc"], "bigmhc_im": ["bigmhc"],
            "netmhciipan": ["iedb_bundle_ii"], "mixmhc2pred": ["mixmhc2pred"], "mhcnuggets": ["mhcnuggets"], "cd4episcore": ["iedb_bundle_ii"]}[m]


def preview(params):
    p = normalise(params)
    bounds = (8, 14) if p["cls"] == 1 else pipeline_ii.PEPTIDE_BOUNDS
    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"], bounds)
    n = int(peps["Peptide"].nunique()) if not peps.empty else 0
    na = len(p["alleles"]) if MODELS[p["model"]][2] else 1
    return {"n_peptides": n, "n_windows": int(len(peps)), "n_alleles": na, "n_predictions": n * na, "warnings": warns,
            "too_big": n * na > pipeline.MAX_PAIRS}


def run(params, job):
    t0 = time.time()
    p = normalise(params)
    model, cls = p["model"], p["cls"]
    label = MODELS[model][1]
    bounds = (8, 14) if cls == 1 else pipeline_ii.PEPTIDE_BOUNDS
    peps, warns = nm.build_peptides(p["text"], p["mode"], p["lengths"], bounds)
    if peps.empty:
        raise ValueError("No valid peptides in the input.")
    peps = peps.drop_duplicates()
    seqs = nm.read_sequences(p["text"]) if p["mode"] == "scan" else None
    uniq = peps["Peptide"].drop_duplicates().tolist()
    al = p["alleles"]
    if len(uniq) * max(len(al), 1) > pipeline.MAX_PAIRS:
        raise ValueError(f"{len(uniq):,} peptides x {len(al)} alleles is over the {pipeline.MAX_PAIRS:,} limit.")
    stage = "netmhciipan" if model == "netmhciipan" else model
    job.stage(stage, label, weight=1.0)
    job.log(f"{label}: {len(uniq):,} unique peptides, {len(al)} alleles")

    if model == "netmhcpan":
        known = set(nm.list_alleles())
        ok = [a for a in al if not known or allele_lib.netmhcpan_name(a) in known]
        if len(ok) < len(al):
            warns.append("NetMHCpan does not know: " + ", ".join(a for a in al if a not in ok))
        res, logs = nm.run_netmhcpan(uniq, [allele_lib.netmhcpan_name(a) for a in ok], {**p["nmp"], "t": None, "raw": ""}, job=job)
        res["Allele"] = [allele_lib.canonical(a) for a in res["Allele"]]
    elif model == "mixmhcpred":
        res, n = mx.predict(uniq, al, job=job); warns += n
    elif model == "prime":
        res, n = pr.predict(uniq, al, job=job); warns += n
    elif model == "iedb_imm":
        res, n = ci.predict(uniq, al, p["iedb"]["mask_choice"], p["iedb"]["positions"]); warns += n
    elif model == "mhcflurry":
        res, n = mf.predict(peps, al, sequences=seqs, job=job); warns += n
    elif model == "netctlpan":
        res, n = ip.netctlpan(seqs, p["lengths"], al, job=job); warns += n
    elif model == "netchop":
        res = ip.netchop_cterm(seqs, peps, job=job)
    elif model in ("bigmhc_el", "bigmhc_im"):
        human = [a for a in al if a.startswith("HLA-")]
        if len(human) < len(al):
            warns.append("BigMHC is human HLA-I only, skipping: " + ", ".join(a for a in al if a not in human))
        grid = pd.DataFrame([(x, a) for x in uniq for a in human], columns=["peptide", "allele"])
        res = bm.predict_pairs(grid, (model.split("_")[1],), job=job).rename(columns={"peptide": "Peptide", "allele": "Allele"})
    elif model == "netmhciipan":
        res, n = nn.predict(peps, al, p["lengths"] or sorted({len(x) for x in uniq}), seqs, ("el", "ba"), p["version"], job=job); warns += n
    elif model == "mixmhc2pred":
        res, n = m2.predict(peps, al, seqs, job=job); warns += n
    elif model == "mhcnuggets":
        res, n = mn.predict(uniq, al, job=job); warns += n
    else:
        res = cd.predict(uniq, job=job)
    if res is None or res.empty:
        raise RuntimeError(f"{label} returned no results (see the log).")
    job.finish_stage(stage, f"{len(res):,} rows")

    on = [c for c in ("Protein", "Start", "Peptide") if c in res.columns]
    table = peps.merge(res, on=on, how="inner")
    table.insert(3, "Length", table["Peptide"].str.len())
    main, asc = MODELS[model][4], MODELS[model][5]
    if main in table:
        table = table.sort_values(main, ascending=asc, na_position="last", kind="stable")
    table = table.reset_index(drop=True)
    path = f"{config.run_dir(job.id)}/{job.id}_{model}.csv"
    table.round(6).to_csv(path, index=False)
    job.result = {"single": True, "model": model, "table_csv": path, "main_column": main,
                  "summary": {"model": model, "label": label, "peptides": len(uniq), "pairs": int(len(table)), "pairs_passed": 0, "peptides_passed": 0,
                              "seconds": round(time.time() - t0, 1), "warnings": list(dict.fromkeys(warns)), "chained": False}}
    job.log(f"done in {time.time() - t0:.1f}s: {len(table):,} rows")
