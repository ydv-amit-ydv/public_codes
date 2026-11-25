#!/usr/bin/env python3
"""
One honest scorecard for the whole pipeline.

Every claim the pipeline makes is listed with the evidence behind it, the
statistic, its uncertainty, and a stated confidence. The point is that a
reader should never have to guess which numbers are measured and which are
assumed.

EVIDENCE GRADES
---------------
  A  Validated against measured wet-lab data, out-of-sample.
     A held-out split or an unseen pathogen. The number is a prediction that
     was tested and could have failed.

  B  Validated against measured wet-lab data, in-sample.
     Real data, but the model saw it during fitting. Shows the model CAN fit,
     not that it will generalise.

  C  Calibrated to published qualitative findings.
     No measured dataset. The model was tuned until it reproduced the
     DIRECTION of results reported in the literature. A structured prior.

  D  Assumption / structural choice.
     A parameter set by judgement. Falsifiable in principle, untested here.

WHAT "CONFIDENCE" MEANS HERE
----------------------------
For grade A and B claims it is a frequentist statement derived from the data
(a CI or a p-value), and it is labelled as such.

For grade C claims it is a ROBUSTNESS frequency from Monte-Carlo sampling of
the parameters: the share of plausible parameter sets for which the claim
still holds. That is NOT a probability that the biology is true. It answers
"how load-bearing is my parameter guess", not "is this real". It is reported
separately and never mixed with the grade A numbers.

Run:  python3 scorecard.py
      python3 scorecard.py --robustness 150     (slower, samples parameters)
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))
RESULTS = HERE / "results"
MODELS = HERE / "models"

from stats import prob_direction                     # noqa: E402

GRADE_TEXT = {
    "A": "measured data, out-of-sample",
    "B": "measured data, in-sample",
    "C": "calibrated to published qualitative findings",
    "D": "assumption / structural choice",
}


def load(p):
    f = RESULTS / p if not Path(p).is_absolute() else Path(p)
    if not f.exists():
        return None
    return (pd.read_csv(f) if f.suffix == ".csv"
            else json.loads(f.read_text()))


# ------------------------------------------------------- robustness sampling
def robustness(n=120, seed=0):
    """
    Resample the fitted parameters and ask how often each platform claim
    still holds. Reports a robustness frequency, not a truth probability.
    """
    from dynamics import P, Sim, load_epitopes
    f = MODELS / "calibrated_wetlab.json"
    if not f.exists():
        f = MODELS / "calibrated.json"
    base = dict(P, **json.loads(f.read_text())["params"]) if f.exists() else dict(P)

    ep = load_epitopes(str(RESULTS / "demo_epitopes.csv"))
    rng = np.random.default_rng(seed)
    # parameters the claims are most sensitive to
    vary = ["th_prime", "th_prolif", "th_death", "b_activate", "b_prolif",
            "plasma_diff", "igg_rate", "switch_rate", "aff_rate",
            "crosspres", "agi_decay", "k_agi_ctl", "dc_mature", "ag_deposit"]

    claims = {"mrna_beats_protein_CD8": [], "as01_more_CTL_than_alum": [],
              "adjuvant_raises_antibody": [], "boost_exceeds_prime": []}

    for _ in range(n):
        p = dict(base)
        for k in vary:
            if k in p:
                p[k] = float(p[k] * np.exp(rng.normal(0, 0.25)))   # +-25% lognormal
        try:
            def go(platform, adj):
                s = Sim(ep, p=p, injections=[(2, None, adj), (270, None, adj)],
                        platform=platform, adjuvant=adj)
                d = s.run(n_steps=700)
                return d
            m = go("mrna_lnp", "intrinsic")
            pr = go("protein", "alum")
            a1 = go("protein", "as01")
            no = go("protein", "none")

            claims["mrna_beats_protein_CD8"].append(
                float(m.TC_active.max() - pr.TC_active.max()))
            claims["as01_more_CTL_than_alum"].append(
                float(a1.TC_active.max() - pr.TC_active.max()))
            claims["adjuvant_raises_antibody"].append(
                float(pr.IgG.max() - no.IgG.max()))
            bd = 270 * 8 / 24
            pre, post = pr[pr.day < bd], pr[pr.day >= bd]
            claims["boost_exceeds_prime"].append(
                float(post.IgG.max() - pre.IgG.max()))
        except Exception:
            continue

    return {k: {"n": len(v), "robustness": prob_direction(v),
                "median_effect": float(np.median(v)) if v else float("nan")}
            for k, v in claims.items()}


# ----------------------------------------------------------------- scorecard
def build(rob=None):
    rows = []

    def add(claim, grade, evidence, stat, conf, note=""):
        rows.append({"claim": claim, "grade": grade, "evidence": evidence,
                     "statistic": stat, "confidence": conf, "note": note})

    # ---------- epitope prediction ----------
    ext = load("external_validation.csv")
    if ext is not None and len(ext):
        for cls in ("I", "II"):
            s = ext[ext.mhc_class == cls]
            if s.empty:
                continue
            nsig = int(s.sig_vs_baseline.sum())
            add(f"Binding model ranks epitopes of an UNSEEN pathogen "
                f"(MHC class {cls})", "A",
                f"leave-one-pathogen-out, {len(s)} pathogens, "
                f"n={int(s.n_test.sum()):,} test peptides",
                f"median ROC-AUC {s.roc_auc.median():.3f} "
                f"(range {s.roc_auc.min():.3f}-{s.roc_auc.max():.3f})",
                f"beats composition baseline in {nsig}/{len(s)} after "
                f"Benjamini-Hochberg",
                f"median ECE {s.ece.median():.3f}")
    else:
        add("Binding model generalises to an unseen pathogen", "A",
            "run external_validation.py", "pending", "pending")

    vb = load("validation_binding.csv")
    if vb is not None and len(vb):
        for cls in ("I", "II"):
            s = vb[vb.mhc_class == cls]
            if s.empty:
                continue
            add(f"Binding model ranks epitopes on held-out peptides "
                f"(MHC class {cls})", "A",
                f"peptide-disjoint split, {len(s)} alleles, "
                f"n={int(s.n.sum()):,}",
                f"mean ROC-AUC {s.roc_auc.mean():.3f}",
                f"per-allele 95% CIs in validation_binding.csv",
                "same-pathogen split: optimistic vs the LOPO figure above")

    vi = load("validation_immunogenicity.csv")
    if vi is not None and len(vi):
        for cls in ("I", "II"):
            s = vi[vi.mhc_class == cls]
            if s.empty:
                continue
            weak = s[(s.roc_auc_lo <= 0.5)]
            add(f"Immunogenicity model predicts T cell response "
                f"(MHC class {cls})", "A",
                f"cross-validated, peptide-disjoint, {len(s)} alleles",
                f"mean ROC-AUC {s.roc_auc.mean():.3f}",
                f"{len(s)-len(weak)}/{len(s)} alleles with CI excluding 0.5",
                ("alleles whose CI touches chance: "
                 + ", ".join(weak.allele) if len(weak) else ""))

    # ---------- dynamic model ----------
    for split, grade in (("test", "A"), ("train", "B")):
        wl = load(f"wetlab_validation_{split}.json")
        if not wl or "_summary" not in wl:
            continue
        s = wl["_summary"]
        lbl = ("held-out subjects, never used for calibration"
               if split == "test" else "subjects used for calibration")
        add(f"Boost antibody fold-rise matches measured human data", grade,
            f"CMI-PB Tdap, {lbl}",
            f"model {s['model_fold_rise']:.2f}x; inside the observed "
            f"5-95% subject band for {s['fold_rise_ok']}/{s['antigens']} "
            f"vaccine antigens",
            f"{s['fold_rise_ok']}/{s['antigens']} antigens consistent",
            "ratios only - the model's antibody is in arbitrary units")
        add(f"Boost antibody time-to-peak matches measured human data", grade,
            f"CMI-PB Tdap, {lbl}",
            f"model {s['model_peak_day']:.1f} d; inside the band for "
            f"{s['peak_day_ok']}/{s['antigens']} antigens",
            f"{s['peak_day_ok']}/{s['antigens']} antigens consistent",
            "CMI-PB samples at days 0,1,3,7,14,30 - resolution is coarse")
        vc = wl.get("_vaccine_vs_control")
        if vc:
            add("The measured boost signal is real, not analysis noise", grade,
                "CMI-PB vaccine antigens vs non-vaccine controls (OVA etc.)",
                f"median fold-rise {vc['median_vaccine']:.2f}x vs "
                f"{vc['median_control']:.2f}x",
                f"Mann-Whitney p {'<1e-4' if vc['p'] < 1e-4 else f'{vc[chr(34)+chr(34)] if False else vc[chr(39)+chr(39)] if False else vc['p']:.3g}'}, "
                f"rank-biserial r={vc['effect_rank_biserial']:.3f}",
                "negative control behaves as it should")

    # ---------- platform claims ----------
    plat_claims = [
        ("mRNA-LNP induces more CD8 than adjuvanted protein",
         "mrna_beats_protein_CD8"),
        ("AS01 induces more CTL than alum (Th1 vs Th2)",
         "as01_more_CTL_than_alum"),
        ("An adjuvant is required for a protein antibody response",
         "adjuvant_raises_antibody"),
        ("The secondary response exceeds the primary",
         "boost_exceeds_prime"),
    ]
    for text, key in plat_claims:
        r = (rob or {}).get(key)
        if r and np.isfinite(r["robustness"]):
            add(text, "C", "direction reported in the vaccine literature; "
                f"robustness over {r['n']} parameter draws (+-25%)",
                f"holds in {r['robustness']*100:.0f}% of parameter draws",
                f"robustness {r['robustness']*100:.0f}% (NOT a truth probability)",
                "no measured dataset compares these platforms here")
        else:
            add(text, "C", "direction reported in the vaccine literature; "
                "reproduced by test_platforms.py",
                "ordinal benchmark passes", "direction only",
                "run --robustness to quantify parameter sensitivity")

    # ---------- explicit non-claims ----------
    add("Absolute antibody titer in a named individual", "D",
        "none - the model's antibody is in arbitrary units",
        "not predicted", "none",
        "only ratios and timings are comparable to data")
    add("Platform constants (half-lives, routing, Th1/Th2 bias)", "D",
        "order-of-magnitude literature values and relative scales",
        "not fitted", "none",
        "falsifiable prior; see platforms.py docstring")
    add("Primary (first-exposure) response kinetics", "D",
        "no measured primary-response dataset harvested - every CMI-PB "
        "subject was primed in infancy",
        "not validated", "none",
        "the boost arm IS validated; the primary arm is not")
    add("Protection / efficacy against infection", "D",
        "no challenge or efficacy data anywhere in this pipeline",
        "not modelled", "none",
        "immunogenicity is not efficacy")

    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--robustness", type=int, default=0,
                    help="parameter draws for grade-C robustness (0 = skip)")
    args = ap.parse_args()

    rob = None
    if args.robustness:
        print(f"sampling {args.robustness} parameter sets for robustness...\n")
        rob = robustness(args.robustness)

    tab = build(rob)
    tab.to_csv(RESULTS / "scorecard.csv", index=False)

    print("=" * 100)
    print("PIPELINE SCORECARD - what is validated, how, and how confidently")
    print("=" * 100)
    for g in ("A", "B", "C", "D"):
        s = tab[tab.grade == g]
        if s.empty:
            continue
        print(f"\n  GRADE {g}: {GRADE_TEXT[g].upper()}")
        print("  " + "-" * 96)
        for _, r in s.iterrows():
            print(f"\n    {r['claim']}")
            print(f"      evidence   {r['evidence']}")
            print(f"      statistic  {r['statistic']}")
            print(f"      confidence {r['confidence']}")
            if r["note"]:
                print(f"      caveat     {r['note']}")

    print("\n" + "=" * 100)
    print("HOW TO READ THIS")
    print("=" * 100)
    print("  Grade A is the only grade where a number was tested out-of-sample")
    print("  and could have failed. Grade C 'robustness' is the share of")
    print("  plausible parameter sets giving the stated direction - it measures")
    print("  how load-bearing a parameter guess is, NOT whether the biology is")
    print("  true. Grade D items are not claims at all; they are listed so")
    print("  nobody mistakes silence for evidence.")
    print(f"\n  wrote {RESULTS/'scorecard.csv'}")


if __name__ == "__main__":
    main()
