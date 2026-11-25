#!/usr/bin/env python3
"""
Local C-ImmSim-like dynamic immune model, driven by our own predicted epitopes.

No external server. The epitope predictions from predict.py set the
immunodominance weights, and this module integrates the cellular response.

Compartments (per C-ImmSim's structure, as ODEs rather than agents):
    Ag            antigen
    D             danger signal (adjuvant / LPS)
    DC            dendritic cells: immature -> mature/presenting
    MA            macrophages
    TH            T helper   : naive -> activated -> memory   (class II driven)
    TC            cytotoxic T: naive -> activated -> memory   (class I driven)
    B             B cells    : naive -> activated -> plasma   (+ memory)
    IgM, IgG      antibody, with affinity maturation
    IL2           interleukin-2

Each predicted epitope becomes a clone with an immunodominance weight derived
from its predicted IC50 (and, where available, its predicted immunogenicity).

Time base matches C-ImmSim: 1 step = 8 h.
"""
import argparse
import json
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

from model import encode, inv_log50k
from platforms import resolve

HERE = Path(__file__).parent
MODELS = HERE / "models"
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

STEP_HOURS = 8.0


# ----------------------------------------------------------------- parameters
P = dict(
    # antigen
    ag_decay       = 0.05,    # /day natural clearance
    ag_clear_ab    = 2.0e-4,  # antibody-mediated clearance
    ag_uptake      = 1.0e-3,  # uptake by APC
    ag_growth      = 0.0,     # >0 for a replicating pathogen

    # follicular-DC antigen depot: immune complexes retained in the
    # germinal centre long after free antigen is cleared. Without this the
    # secondary response is starved, because pre-existing antibody clears
    # the boost within days.
    ag_deposit     = 0.08,    # free antigen -> FDC-bound depot
    agd_decay      = 0.045,   # slow release/decay of the depot

    # Endogenously expressed antigen (RNA/vector platforms). Antigen made
    # inside the cell enters the MHC class I pathway directly; exogenous
    # antigen reaches class I only by inefficient cross-presentation. This
    # split is what separates mRNA-type CD8 responses from protein + alum.
    agi_decay      = 1.20,    # turnover of intracellular antigen
    # Immune control of a replicating template: CTL kill expressing cells and
    # antibody neutralises released virions. Without this an amplifying
    # replicon (saRNA, live attenuated) saturates its cap and expresses
    # forever, which is both wrong and flattering to those platforms.
    r_clear_ctl    = 1.6e-5,  # per activated CTL per day
    r_clear_ab     = 4.0e-6,  # per unit antibody per day
    crosspres      = 0.08,    # exogenous -> class I cross-presentation
    k_agi_ctl      = 300.0,   # endogenous antigen for CTL priming

    # danger signal / adjuvant
    d_decay        = 0.35,

    # antigen presenting cells
    dc_total       = 1.0e3,
    dc_mature      = 3.0e-2,  # maturation rate, scaled by danger signal
    dc_decay       = 0.12,
    ma_total       = 1.0e3,
    ma_activate    = 1.5e-2,

    # T helper (class II)
    th_naive       = 1.0e3,
    th_prime       = 1.5e-3,  # priming by mature DC, scaled by epitope weight
    th_prolif      = 1.70,
    th_death       = 0.48,
    th_memory      = 0.085,
    th_mem_decay   = 0.012,
    th_mem_recall  = 25.0,    # memory responds faster than naive

    # cytotoxic T (class I)
    tc_naive       = 1.0e3,
    tc_prime       = 1.2e-3,
    tc_prolif      = 1.55,
    tc_death       = 0.48,
    tc_memory      = 0.075,
    tc_mem_decay   = 0.012,
    tc_mem_recall  = 25.0,

    # B cells
    b_naive        = 1.0e3,
    b_activate     = 2.0e-2,  # needs antigen AND T help
    b_prolif       = 1.75,
    b_death        = 0.17,
    b_memory       = 0.11,
    b_mem_decay    = 0.010,
    b_mem_recall   = 30.0,

    # --- saturation (half-max) constants -------------------------------
    # These must sit near the cell/molecule numbers the model actually
    # reaches, or a term silently pins to zero and the cascade starves.
    k_ag_apc       = 500.0,   # antigen for APC uptake/maturation
    k_ag_gc        = 40.0,    # antigen for germinal-centre B proliferation;
                              # low, because follicular DCs hold antigen
                              # long after free antigen has cleared
    k_danger       = 50.0,    # danger signal for APC maturation
    k_dc_prolif    = 15.0,    # mature DC for T cell proliferation
    k_il2          = 15.0,    # IL-2 for proliferation
    k_th_ctl       = 150.0,   # T help for CTL proliferation
    k_b_act        = 400.0,   # T help for B activation
    k_b_prolif     = 400.0,   # T help for B proliferation
    k_switch       = 600.0,   # T help for class switching
    k_aff_help     = 600.0,   # T help for affinity maturation
    k_aff_ag       = 60.0,    # depot antigen for affinity maturation

    # Density-dependent limit on duplication (C-ImmSim MaxLocalFact).
    # Without this the IL-2 / T-help feedback loop grows without bound.
    k_carry        = 1.0e7,   # total lymphocyte carrying capacity

    # plasma cells and antibody
    plasma_diff    = 0.30,    # activated B -> plasma
    plasma_death   = 0.35,
    igm_rate       = 6.0,     # IgM secretion per plasma cell
    igg_rate       = 22.0,    # IgG secretion, after class switch
    switch_rate    = 0.16,    # IgM -> IgG class switch, needs T help
    igm_decay      = 0.35,    # IgM half-life ~2 d
    igg_decay      = 0.030,   # IgG half-life ~23 d (C-ImmSim value)

    # affinity maturation
    aff_rate       = 0.09,    # germinal-centre affinity gain
    aff_max        = 10.0,    # max fold improvement over naive affinity

    # cytokine
    il2_rate       = 0.02,
    il2_decay      = 1.2,
)


def hill(x, k):
    return x / (k + x)


# --------------------------------------------------------------- epitope load
def load_epitopes(csv, top_i=4, top_ii=3, use_immunogenicity=True):
    """
    Take predicted epitopes and turn them into clones with weights.

    weight = binding strength (log50k) x predicted immunogenicity (if a model
    for that allele exists). This is the immunodominance hierarchy.
    """
    d = pd.read_csv(csv)
    d = d[d.call != "non-binder"].copy()
    if d.empty:
        raise SystemExit("no binders in the prediction file")

    keep = []
    for mclass, n in (("I", top_i), ("II", top_ii)):
        sub = d[d.mhc_class == mclass]
        if sub.empty:
            continue
        keep.append(sub.sort_values("pred_IC50_nM")
                       .groupby("allele").head(n))
    d = pd.concat(keep, ignore_index=True)

    # binding strength on the log50k scale
    d["bind_w"] = np.clip(1.0 - np.log(d.pred_IC50_nM.clip(1, 5e4)) / np.log(5e4),
                          0.0, 1.0)

    d["immuno_p"] = np.nan
    if use_immunogenicity:
        for allele, sub in d.groupby("allele"):
            f = MODELS / f"immuno_{allele.replace('*','_').replace(':','_')}.pkl"
            if not f.exists():
                continue
            with open(f, "rb") as fh:
                im = pickle.load(fh)
            L = im["length"]
            ok = sub[sub.peptide.str.len() == L]
            if ok.empty:
                continue
            bm = MODELS / f"{allele.replace('*','_').replace(':','_')}.pkl"
            with open(bm, "rb") as fh:
                bmod = pickle.load(fh)
            aff = bmod["model"].predict(encode(ok.peptide.values, L))
            X = np.hstack([encode(ok.peptide.values, L), aff.reshape(-1, 1)])
            d.loc[ok.index, "immuno_p"] = im["model"].predict_proba(X)[:, 1]

    # Combined weight. Where no immunogenicity model exists, fall back to
    # binding alone rather than inventing a probability.
    d["weight"] = d.bind_w * np.where(d.immuno_p.notna(), d.immuno_p, 1.0)
    d["weight"] = d.weight / d.weight.max()
    return d.sort_values(["mhc_class", "weight"], ascending=[True, False])


# -------------------------------------------------------------------- the ODE
class Sim:
    """
    One simulated host receiving one vaccine regimen.

    platform / adjuvant come from platforms.py and set how antigen is supplied
    and how strongly innate immunity is engaged. injections are
    (step, dose_ug, adjuvant_override_or_None).
    """

    def __init__(self, epitopes, p=None, injections=((2, None, None),),
                 platform="protein", adjuvant=None, dose_ug=None,
                 heterologous=None):
        self.ep = epitopes.reset_index(drop=True)
        self.p = dict(P, **(p or {}))
        self.pf, self.adj, self.dose_ug = resolve(platform, adjuvant, dose_ug)
        self.platform_key = platform
        self.heterologous = heterologous      # platform key used for the boost

        # Template decay from the published half-life, and a cap so an
        # amplifying replicon saturates rather than growing without bound.
        self.r_decay = (np.log(2.0) / (self.pf.half_life_h / 24.0)
                        if self.pf.template else 0.0)
        self.r_cap = 4.0e3

        # Adjuvant acts on three things: how much danger signal is present,
        # how hard CD8 priming is pushed (Th1 bias), and how hard B cells are
        # pushed (Th2 bias). A Th1-biased adjuvant such as AS01 raises CTL and
        # leaves antibody roughly alone; alum does the opposite.
        bias = self.adj.th1_bias
        self.ctl_gain = self.adj.ctl_boost * (0.55 + 0.90 * bias)
        self.ab_gain = self.adj.ab_boost * (0.55 + 0.90 * (1.0 - bias))

        # Danger signal decays over whichever lasts longer: the adjuvant or
        # the platform's own innate activation.
        innate_days = max(self.adj.duration_days, self.pf.innate_days, 0.5)
        self.d_decay = np.log(2.0) / innate_days

        # Depot adjuvants (alum) slow clearance of extracellular antigen.
        self.p["ag_decay"] = self.p["ag_decay"] * (1.0 - 0.75 * self.adj.depot)

        self.injections = self._build(injections)

        self.ci = self.ep.index[self.ep.mhc_class == "I"].tolist()
        self.cii = self.ep.index[self.ep.mhc_class == "II"].tolist()
        self.nI, self.nII = len(self.ci), len(self.cii)
        self.nB = max(1, self.nII)

    # --- dose -> compartments ------------------------------------------
    AG_PER_UG = 40.0          # antigen units per ug of protein
    TEMPLATE_PER_UG = 9.0     # template units per ug of RNA / vector

    def _build(self, injections):
        """Turn (step, dose_ug, adjuvant) into concrete boluses."""
        out = []
        for n, item in enumerate(injections):
            step = item[0]
            dose = item[1] if len(item) > 1 and item[1] is not None else self.dose_ug
            adj_key = item[2] if len(item) > 2 else None
            adj = self.adj
            if adj_key is not None:
                _, adj, _ = resolve(self.platform_key, adj_key, dose)

            pf = self.pf
            # A homologous boost with a viral vector is blunted by anti-vector
            # immunity; a heterologous boost is not.
            vec_loss = 1.0
            if n > 0 and pf.anti_vector > 0 and not self.heterologous:
                vec_loss = 1.0 - pf.anti_vector

            danger = (adj.magnitude + pf.innate * 100.0) * vec_loss
            if pf.template:
                out.append({"step": step, "ag": 0.0,
                            "template": dose * self.TEMPLATE_PER_UG * vec_loss,
                            "danger": danger})
            else:
                out.append({"step": step, "ag": dose * self.AG_PER_UG,
                            "template": 0.0, "danger": danger})
        return out

    # state vector layout
    def _unpack(self, y):
        i = 0
        def take(n):
            nonlocal i
            v = y[i:i + n]; i += n
            return v
        Ag   = take(1)[0]
        Agd  = take(1)[0]
        R    = take(1)[0]
        AgI  = take(1)[0]
        D    = take(1)[0]
        DCm  = take(1)[0]
        MAa  = take(1)[0]
        THa  = take(self.nII); THm = take(self.nII)
        TCa  = take(self.nI);  TCm = take(self.nI)
        Ba   = take(self.nB);  Bm  = take(self.nB)
        PL   = take(self.nB)
        IgM  = take(1)[0]
        IgG  = take(1)[0]
        Aff  = take(1)[0]
        IL2  = take(1)[0]
        return (Ag, Agd, R, AgI, D, DCm, MAa, THa, THm, TCa, TCm,
                Ba, Bm, PL, IgM, IgG, Aff, IL2)

    def y0(self):
        return np.concatenate([
            [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0],
            np.zeros(self.nII), np.zeros(self.nII),
            np.zeros(self.nI), np.zeros(self.nI),
            np.zeros(self.nB), np.zeros(self.nB), np.zeros(self.nB),
            [0.0], [0.0], [1.0], [0.0],
        ])

    def rhs(self, t, y):
        p = self.p
        (Ag, Agd, R, AgI, D, DCm, MAa, THa, THm, TCa, TCm,
         Ba, Bm, PL, IgM, IgG, Aff, IL2) = self._unpack(y)

        Ag = max(Ag, 0.0); Agd = max(Agd, 0.0)
        R = max(R, 0.0); AgI = max(AgI, 0.0); D = max(D, 0.0)
        total_lym = float(THa.sum() + TCa.sum() + Ba.sum() + PL.sum()
                          + THm.sum() + TCm.sum() + Bm.sum())
        crowd = max(0.0, 1.0 - total_lym / p["k_carry"])

        wI = self.ep.weight.values[self.ci] if self.nI else np.array([])
        wII = self.ep.weight.values[self.cii] if self.nII else np.array([])

        # --- antigen ---
        ab_eff = IgM + IgG * Aff
        dAg = (p["ag_growth"] * Ag
               - p["ag_decay"] * Ag
               - p["ag_clear_ab"] * Ag * ab_eff
               - p["ag_uptake"] * Ag * (DCm + MAa + p["dc_total"] * 0.1))
        # Template (mRNA / circRNA / vector genome): decays at its own rate,
        # optionally amplifying first (saRNA replicase, live attenuated).
        ctl_total = float(TCa.sum())
        dR = (self.pf.amplifies * R * max(0.0, 1.0 - R / self.r_cap)
              - self.r_decay * R
              - (p["r_clear_ctl"] * ctl_total
                 + p["r_clear_ab"] * (IgM + IgG * Aff)) * R)
        made = self.pf.translation * R                    # antigen translated
        # Endogenous pool feeds MHC class I; a share is secreted and joins the
        # extracellular pool where B cells and class II can see it.
        dAgI = (made * self.pf.mhc1_routing
                + p["crosspres"] * p["ag_uptake"] * Ag * (DCm + MAa)
                - p["agi_decay"] * AgI)
        dAg = dAg + made * self.pf.secreted
        dAg = dAg if Ag > 0 else max(dAg, 0.0)
        dAgd = p["ag_deposit"] * Ag - p["agd_decay"] * Agd
        dD = -self.d_decay * D

        # --- APC: maturation needs antigen AND danger signal ---
        free_dc = max(p["dc_total"] - DCm, 0.0)
        ag_ii = Ag * self.pf.mhc2_routing
        dDCm = (p["dc_mature"] * free_dc * hill(ag_ii, p["k_ag_apc"]) * (0.15 + hill(D, p["k_danger"]))
                - p["dc_decay"] * DCm)
        free_ma = max(p["ma_total"] - MAa, 0.0)
        dMAa = (p["ma_activate"] * free_ma * hill(ag_ii, p["k_ag_apc"]) * (0.2 + hill(D, p["k_danger"]))
                - p["dc_decay"] * MAa)

        # --- T helper, class II restricted ---
        help_total = float(THa.sum())
        dTHa = np.zeros(self.nII); dTHm = np.zeros(self.nII)
        for k in range(self.nII):
            prime = (p["th_prime"] * wII[k] * DCm
                     * (p["th_naive"] + p["th_mem_recall"] * THm[k]))
            prolif = (p["th_prolif"] * THa[k] * hill(IL2, p["k_il2"])
                      * hill(DCm, p["k_dc_prolif"]) * crowd)
            dTHa[k] = prime + prolif - p["th_death"] * THa[k]
            dTHm[k] = p["th_memory"] * THa[k] - p["th_mem_decay"] * THm[k]

        # --- cytotoxic T, class I restricted, needs help ---
        dTCa = np.zeros(self.nI); dTCm = np.zeros(self.nI)
        hfac = hill(help_total, p["k_th_ctl"])
        for k in range(self.nI):
            prime = (p["tc_prime"] * wI[k] * DCm * self.ctl_gain
                     * hill(AgI, p["k_agi_ctl"])
                     * (p["tc_naive"] + p["tc_mem_recall"] * TCm[k]))
            prolif = (p["tc_prolif"] * TCa[k] * hill(IL2, p["k_il2"])
                      * hill(DCm, p["k_dc_prolif"])
                      * (0.25 + 0.75 * hfac) * crowd)
            dTCa[k] = prime + prolif - p["tc_death"] * TCa[k]
            dTCm[k] = p["tc_memory"] * TCa[k] - p["tc_mem_decay"] * TCm[k]

        # --- B cells: antigen + T help ---
        dBa = np.zeros(self.nB); dBm = np.zeros(self.nB); dPL = np.zeros(self.nB)
        for k in range(self.nB):
            w = wII[k] if self.nII else 1.0
            th_k = THa[k] if self.nII else help_total
            act = (p["b_activate"] * w * self.ab_gain * hill(ag_ii, p["k_ag_apc"])
                   * (p["b_naive"] + p["b_mem_recall"] * Bm[k])
                   * hill(th_k, p["k_b_act"]))
            prolif = (p["b_prolif"] * Ba[k] * hill(th_k, p["k_b_prolif"])
                      * hill(IL2, p["k_il2"]) * hill(Agd, p["k_ag_gc"]) * crowd)
            dBa[k] = act + prolif - p["b_death"] * Ba[k] - p["plasma_diff"] * Ba[k]
            dBm[k] = p["b_memory"] * Ba[k] - p["b_mem_decay"] * Bm[k]
            dPL[k] = p["plasma_diff"] * Ba[k] - p["plasma_death"] * PL[k]

        PLs = float(PL.sum())
        # --- antibody: IgM first, IgG after class switch (T-help dependent) ---
        switch = p["switch_rate"] * hill(help_total, p["k_switch"])
        dIgM = p["igm_rate"] * PLs * (1.0 - switch) - p["igm_decay"] * IgM
        dIgG = p["igg_rate"] * PLs * switch - p["igg_decay"] * IgG

        # --- affinity maturation: driven by ongoing GC reaction ---
        dAff = (p["aff_rate"] * hill(help_total, p["k_aff_help"]) * hill(Agd, p["k_aff_ag"])
                * (p["aff_max"] - Aff) / p["aff_max"])

        # --- IL-2 from activated T helper ---
        dIL2 = p["il2_rate"] * help_total - p["il2_decay"] * IL2

        return np.concatenate([[dAg], [dAgd], [dR], [dAgI], [dD], [dDCm], [dMAa],
                               dTHa, dTHm, dTCa, dTCm,
                               dBa, dBm, dPL, [dIgM], [dIgG], [dAff], [dIL2]])

    def run(self, n_steps=540):
        """Integrate, applying injections as discontinuous boluses."""
        days_total = n_steps * STEP_HOURS / 24.0
        events = sorted(self.injections, key=lambda e: e["step"])
        y = self.y0()
        ts, ys = [], []
        t0 = 0.0

        segs = [(e["step"] * STEP_HOURS / 24.0, e) for e in events]
        segs.append((days_total, None))

        for day, ev in segs:
            if day > t0:
                sol = solve_ivp(self.rhs, (t0, day), y, method="LSODA",
                                dense_output=True, rtol=1e-6, atol=1e-8,
                                max_step=0.5)
                grid = np.arange(t0, day, 8.0 / 24.0)
                if len(grid):
                    ts.append(grid); ys.append(sol.sol(grid))
                y = sol.y[:, -1]
                t0 = day
            if ev is not None:
                y[0] += ev["ag"]        # extracellular antigen bolus
                y[2] += ev["template"]  # mRNA / circRNA / vector genome
                y[4] += ev["danger"]    # adjuvant + platform innate signal

        T = np.concatenate(ts); Y = np.hstack(ys)
        return self.to_frame(T, Y)

    def to_frame(self, T, Y):
        cols = {"day": T, "step": T * 24.0 / STEP_HOURS}
        (Ag, Agd, R, AgI, D, DCm, MAa, THa, THm, TCa, TCm,
         Ba, Bm, PL, IgM, IgG, Aff, IL2) = self._unpack(Y)
        cols.update({
            "antigen": Ag, "antigen_depot": Agd, "template": R,
            "antigen_endogenous": AgI, "danger": D, "DC_mature": DCm, "MA_active": MAa,
            "TH_active": THa.sum(axis=0), "TH_memory": THm.sum(axis=0),
            "TC_active": TCa.sum(axis=0), "TC_memory": TCm.sum(axis=0),
            "B_active": Ba.sum(axis=0), "B_memory": Bm.sum(axis=0),
            "plasma": PL.sum(axis=0),
            "IgM": IgM, "IgG": IgG, "affinity": Aff, "IL2": IL2,
            "Ab_total": IgM + IgG,
        })
        for n, k in enumerate(self.cii):
            cols[f"TH_clone_{self.ep.peptide[k]}"] = THa[n]
        for n, k in enumerate(self.ci):
            cols[f"TC_clone_{self.ep.peptide[k]}"] = TCa[n]
        return pd.DataFrame(cols)
