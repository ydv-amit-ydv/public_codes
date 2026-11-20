"""One place that knows how to run every model on a set of peptide-allele pairs.

The single-model tab and the manual funnel both call run(); the automatic pipelines keep their own (equivalent) calls.
`rows` is the set of pairs to score: a DataFrame[Protein, Start, Peptide, Allele]. Tools that take a peptide list and an allele
list (NetMHCpan, MixMHCpred, PRIME, MHCflurry, NetMHCIIpan, MixMHC2pred) are given the peptides and alleles present in `rows`;
BigMHC scores exactly the pairs. The caller keeps only the results of the pairs it asked for.

Import after engine_server has put the tool folders on sys.path.
"""
from dataclasses import dataclass

import pandas as pd

from common import alleles as allele_lib
import bigmhc_runner as bm
import calis_immunogenicity as ci
import cd4episcore_runner as cd
import immuscope_runner as isc
import iedb_processing as ip
import mhcflurry_runner as mf
import mixmhc2pred_runner as m2
import mixmhcpred_runner as mx
import netmhcpan_app as nm
import netmhciipan_runner as nn
import prime_runner as pr
import tlimmuno2_runner as tl2


@dataclass(frozen=True)
class Spec:
    key: str
    cls: int
    label: str
    alleles: bool          # scores peptide-allele pairs (False: the peptide alone)
    protein: bool          # needs the protein around the peptide
    cost: int              # 1 (seconds) .. 5 (slow): used to propose an order
    main: str              # column that ranks its results
    asc: bool              # lower is better
    health: tuple          # health-check names it depends on
    weight: float          # share of a run's progress bar


SPECS = {s.key: s for s in [
    Spec("mixmhcpred", 1, "MixMHCpred 3.0", True, False, 1, "MixMHC_%Rank", True, ("mixmhcpred",), 0.5),
    Spec("prime", 1, "PRIME 2.1", True, False, 1, "PRIME_%Rank", True, ("mixmhcpred", "prime"), 0.5),
    Spec("iedb_imm", 1, "IEDB immunogenicity (Calis)", True, False, 1, "IEDB_Immunogenicity", False, (), 0.2),
    Spec("netchop", 1, "NetChop 3.1", False, True, 2, "NetChop_Cterm", False, ("iedb_bundle",), 0.2),
    Spec("netctlpan", 1, "NetCTLpan 1.1", True, True, 3, "CTLpan_%Rank", True, ("iedb_bundle",), 0.5),
    Spec("netmhcpan", 1, "NetMHCpan 4.2", True, False, 3, "NMP_%Rank_EL", True, ("netmhcpan",), 1.0),
    Spec("mhcflurry", 1, "MHCflurry 2", True, False, 3, "MHCF_Presentation_%Rank", True, ("mhcflurry",), 1.0),
    Spec("bigmhc_el", 1, "BigMHC EL", True, False, 5, "BigMHC_EL", False, ("bigmhc",), 2.0),
    Spec("bigmhc_im", 1, "BigMHC IM", True, False, 5, "BigMHC_IM", False, ("bigmhc",), 2.0),
    Spec("mixmhc2pred", 2, "MixMHC2pred 2.1", True, False, 1, "MIX2_%Rank", True, ("mixmhc2pred",), 0.7),
    Spec("netmhciipan_el", 2, "NetMHCIIpan (EL)", True, False, 3, "NII_%Rank_EL", True, ("iedb_bundle_ii",), 1.0),
    Spec("netmhciipan_ba", 2, "NetMHCIIpan (BA)", True, False, 3, "NII_%Rank_BA", True, ("iedb_bundle_ii",), 1.0),
    Spec("cd4episcore", 2, "CD4episcore", False, False, 3, "CD4_Episcore", True, ("iedb_bundle_ii",), 0.5),
    Spec("tlimmuno2", 2, "TLimmuno2", True, False, 4, "TLimmuno2", False, ("tlimmuno2",), 1.0),
    Spec("immuscope", 2, "ImmuScope-IM", True, False, 4, "ImmuScope_IM", False, ("immuscope",), 1.0),
]}
ORDER = {1: [k for k, s in SPECS.items() if s.cls == 1], 2: [k for k, s in SPECS.items() if s.cls == 2]}


def by_cost(keys):
    return sorted(keys, key=lambda k: (SPECS[k].cost, list(SPECS).index(k)))


def run(keys, rows, ctx, job=None):
    """Score `rows` with the models in `keys` (a stage). ctx: {"seqs": [(name, seq)] or None, "p": normalised params}.
    -> (DataFrame keyed by Protein/Start/Peptide/Allele where the tool knows them, notes). Models of one call that share a
    tool (NetMHCIIpan EL + BA) run once."""
    p, seqs = ctx["p"], ctx.get("seqs")
    peps = rows[["Protein", "Start", "Peptide"]].drop_duplicates()
    uniq = peps["Peptide"].drop_duplicates().tolist()
    als = sorted(rows["Allele"].unique()) if "Allele" in rows else []
    lengths = p.get("lengths") or sorted({len(x) for x in uniq})
    notes, frames = [], []
    for key in keys:
        if key == "netmhciipan_ba" and "netmhciipan_el" in keys:
            continue                                     # one call returns both
        if key == "netmhciipan_el":
            kinds = tuple(k for k, m in (("el", "netmhciipan_el"), ("ba", "netmhciipan_ba")) if m in keys)
            df, n = nn.predict(peps, als, lengths, seqs, kinds, p.get("version", "4.3"), job=job)
        elif key == "netmhciipan_ba":
            df, n = nn.predict(peps, als, lengths, seqs, ("ba",), p.get("version", "4.3"), job=job)
        elif key == "mixmhc2pred":
            df, n = m2.predict(peps, als, seqs, job=job)
        elif key == "cd4episcore":
            df, n = cd.predict(uniq, job=job), []
        elif key in ("tlimmuno2", "immuscope"):
            runner = tl2 if key == "tlimmuno2" else isc
            df, n = runner.predict(rows[["Peptide", "Allele"]].drop_duplicates(), job=job)
        elif key == "netmhcpan":
            known = set(nm.list_alleles())
            ok = [a for a in als if not known or allele_lib.netmhcpan_name(a) in known]
            n = ["NetMHCpan does not know: " + ", ".join(a for a in als if a not in ok)] if len(ok) < len(als) else []
            df, _ = nm.run_netmhcpan(uniq, [allele_lib.netmhcpan_name(a) for a in ok], {**p["nmp"], "t": None, "raw": ""}, job=job)
            if not df.empty:
                df["Allele"] = [allele_lib.canonical(a) for a in df["Allele"]]
        elif key == "mixmhcpred":
            df, n = mx.predict(uniq, als, job=job)
        elif key == "prime":
            df, n = pr.predict(uniq, als, job=job)
        elif key == "iedb_imm":
            df, n = ci.predict(uniq, als, p["iedb"]["mask_choice"], p["iedb"]["positions"])
        elif key == "mhcflurry":
            df, n = mf.predict(peps, als, sequences=seqs, job=job)
        elif key == "netctlpan":
            df, n = ip.netctlpan(seqs, lengths, als, job=job)
        elif key == "netchop":
            df, n = ip.netchop_cterm(seqs, peps, job=job), []
        elif key in ("bigmhc_el", "bigmhc_im"):
            human = rows[rows["Allele"].str.startswith("HLA-")][["Peptide", "Allele"]].drop_duplicates()
            n = ["BigMHC is human HLA-I only, skipping: " + ", ".join(a for a in als if not a.startswith("HLA-"))] if len(human) < len(rows[["Peptide", "Allele"]].drop_duplicates()) and any(not a.startswith("HLA-") for a in als) else []
            if human.empty:
                df = pd.DataFrame(columns=["Peptide", "Allele"])
            else:
                df = bm.predict_pairs(human.rename(columns={"Peptide": "peptide", "Allele": "allele"}), (key.split("_")[1],), job=job
                                      ).rename(columns={"peptide": "Peptide", "allele": "Allele"})
        else:
            raise ValueError(f"unknown model {key}")
        notes += n or []
        frames.append(df)
    return frames, notes


def merge_keys(df):
    return [c for c in ("Protein", "Start", "Peptide", "Allele") if c in df.columns]
