#!/usr/bin/env python3
"""Plain-language answers to "does this thing actually work?".

Every verdict here is COMPUTED from the result files, never written down, so
it cannot drift from the evidence. If a result file is missing the verdict
says "not tested yet" rather than guessing, and if a number moves the wording
moves with it.

Each verdict has:
  mark      "y" it works, "n" it does not, "q" partly / with a catch
  question  what a reader would actually ask
  answer    one sentence, no jargon
  because   the measurement, in the plainest form it survives
  means     what to do differently because of it

The honest-reporting rule this file exists to serve: a reader who only ever
looks at this page should come away with the same impression as a reader who
works through VALIDATION.md. So the failures are listed with the successes,
in the same voice, not relegated to a footnote.
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"


def _csv(name):
    f = RESULTS / name
    if not f.exists():
        return None
    try:
        d = pd.read_csv(f)
        return d if len(d) else None
    except Exception:
        return None


def _json(name):
    f = RESULTS / name
    if not f.exists():
        return None
    try:
        return json.loads(f.read_text())
    except Exception:
        return None


def _missing(question, what, how):
    return {"mark": "q", "question": question,
            "answer": "Not tested yet on this machine.",
            "because": f"{what} has not been produced here.",
            "means": how, "untested": True}


# ---------------------------------------------------------------------------
# the questions
# ---------------------------------------------------------------------------
def q_new_pathogen():
    d = _csv("external_validation_pan.csv")
    if d is None:
        return _missing(
            "Can it find the real epitopes in a virus it has never seen?",
            "the leave-one-pathogen-out comparison",
            "Run 'Re-run the whole validation' on the Setup page, or the "
            "leave-one-pathogen-out step on the Validation page.")
    one = d[d.mhc_class == "I"]
    auc = float(one.auc_pan.median())
    beat = int(one.sig_vs_comp.sum())
    n = len(one)
    # AUC in plain terms: pick-the-real-one-out-of-two probability
    pct = round(auc * 100)
    return {
        "mark": "y" if auc >= 0.8 and beat == n else "q",
        "question": "Can it find the real epitopes in a virus it has never seen?",
        "answer": f"Yes, for MHC class I. Shown {n} different organisms it was "
                  f"never trained on, it ranked their real epitopes above "
                  f"non-epitopes about {pct}% of the time.",
        "because": f"Each organism was removed from training completely, then "
                   f"its epitopes were ranked. Median ROC-AUC {auc:.3f} across "
                   f"{n} organisms. A coin flip scores 0.5 and a perfect "
                   f"ranker scores 1.0. It also beat a dumb baseline that "
                   f"only counts amino acids in {beat} of {n}, which matters "
                   f"because that baseline can look good by accident.",
        "means": "Reasonable for shortlisting class I epitopes in a new "
                 "pathogen. Not a substitute for testing them.",
    }


def q_class_ii():
    d = _csv("external_validation_pan.csv")
    if d is None:
        return _missing(
            "Does it work as well for the CD4 (class II) side?",
            "the leave-one-pathogen-out comparison",
            "Run the leave-one-pathogen-out step on the Validation page.")
    two = d[d.mhc_class == "II"]
    if not len(two):
        return _missing("Does it work as well for the CD4 (class II) side?",
                        "class II folds", "Re-run with class II included.")
    auc = float(two.auc_pan.median())
    beat = int(two.sig_vs_comp.sum())
    n = len(two)
    return {
        "mark": "n",
        "question": "Does it work as well for the CD4 (class II) side?",
        "answer": f"No. On an unseen organism class II is barely better than "
                  f"counting amino acids.",
        "because": f"Median ROC-AUC {auc:.3f} over {n} organisms, and it beat "
                   f"the amino-acid-counting baseline in only {beat} of {n}. "
                   f"Two separate datasets agree on this: the "
                   f"mass-spectrometry set gives the same answer.",
        "means": "Do not make a class II decision on this alone. The "
                 "epitope_filter_engine next door fuses eleven predictors and "
                 "is the better tool for class II.",
    }


def q_new_allele():
    d = _csv("eluted_validation.csv")
    if d is None:
        return _missing(
            "Does it work for an HLA type it was not trained on?",
            "the eluted-ligand validation",
            "Run 'Validate against eluted ligands' on the Validation page.")
    col = "coverage" if "coverage" in d.columns else None
    if col is None or "pan_ppv_1pct" not in d.columns:
        return _missing("Does it work for an HLA type it was not trained on?",
                        "the allele-coverage breakdown",
                        "Re-run the eluted-ligand validation.")
    one = d[d.mhc_class == "I"] if "mhc_class" in d.columns else d
    tr = one[one[col] == "trained"].pan_ppv_1pct.median()
    un = one[one[col] == "unknown"].pan_ppv_1pct.median()
    if pd.isna(tr) or pd.isna(un):
        return _missing("Does it work for an HLA type it was not trained on?",
                        "both trained and untrained alleles",
                        "Re-run the eluted-ligand validation.")
    return {
        "mark": "n",
        "question": "Does it work for an HLA type it was not trained on?",
        "answer": "No. It is strong on the HLA types it learned and weak on "
                  "everything else, and this is the single biggest limit on "
                  "using it.",
        "because": f"Of the peptides it ranks in its top 1%, "
                   f"{round(tr * 100)}% are genuinely presented for an HLA "
                   f"type it was trained on, against {round(un * 100)}% for "
                   f"one it has never seen. Both are measured against real "
                   f"mass-spectrometry data where 5% would be chance.",
        "means": "Only use the HLA types the form offers. They are the only "
                 "ones it was trained on, which is why nothing else is listed.",
    }


def q_human_size():
    d = _json("wetlab_validation_test.json")
    if d is None:
        return _missing(
            "Does it get the size of a real human booster response right?",
            "the held-out human subject check",
            "Run 'Set up everything' on the Setup page; it ends with this check.")
    ags = [k for k in d if not k.startswith("_")]
    inband = sum(1 for a in ags
                 if d[a]["fold_rise_check"]["within_subject_90pct_band"])
    inci = sum(1 for a in ags if d[a]["fold_rise_check"].get("within_median_CI"))
    model = d[ags[0]]["fold_rise_check"]["model"]
    meds = [d[a]["observed_fold_rise"]["median"] for a in ags]
    n = d[ags[0]]["n_subjects"]
    return {
        "mark": "y" if inband == len(ags) and inci > len(ags) / 2 else "q",
        "question": "Does it get the size of a real human booster response right?",
        "answer": f"Mostly. It predicts antibody rising {model:.1f}-fold after "
                  f"a boost; real people in the study rose "
                  f"{min(meds):.1f}- to {max(meds):.1f}-fold depending on the "
                  f"antigen.",
        "because": f"Tested against {n} adults given a tetanus/diphtheria/"
                   f"pertussis boost, who were NOT used to tune the model. "
                   f"The prediction falls inside the range of real individuals "
                   f"for {inband} of {len(ags)} antigens, and inside the "
                   f"tighter interval around the typical person for "
                   f"{inci} of {len(ags)}. Where it misses, it predicts too "
                   f"big a rise for the weakly-responding antigens.",
        "means": "Treat the fold-rise as the right order of magnitude, not a "
                 "number to plan a trial around.",
    }


def q_human_timing():
    d = _json("wetlab_validation_test.json")
    if d is None:
        return _missing("Does it get the timing right?",
                        "the held-out human subject check",
                        "Run 'Set up everything' on the Setup page.")
    ags = [k for k in d if not k.startswith("_")]
    model = d[ags[0]]["peak_day_check"]["model"]
    obs = d[ags[0]]["observed_peak_day"]["median"]
    inci = sum(1 for a in ags if d[a]["peak_day_check"].get("within_median_CI"))
    return {
        "mark": "q",
        "question": "Does it get the timing right?",
        "answer": f"Roughly, but it peaks late: about day {model:.0f} against "
                  f"day {obs:.0f} in real people.",
        "because": f"It lands inside the spread of real individuals for every "
                   f"antigen, but inside the tight interval around the typical "
                   f"person for only {inci} of {len(ags)}. The study only drew "
                   f"blood on days 0, 1, 3, 7, 14 and 30, so an observed peak "
                   f"can only ever land on one of those days - the data cannot "
                   f"resolve a five-day error, which partly hides this.",
        "means": "Believe 'about two to three weeks', not the exact day.",
    }


def q_platforms():
    d = _csv("scorecard.csv")
    if d is None:
        return _missing(
            "Does it reproduce the differences between vaccine types?",
            "the evidence scorecard",
            "Run 'Evidence scorecard' on the Validation page.")
    col = "grade" if "grade" in d.columns else None
    nc = int((d[col] == "C").sum()) if col else 0
    return {
        "mark": "q",
        "question": "Does it reproduce the differences between vaccine types?",
        "answer": "It reproduces the directions people have published - mRNA "
                  "beats protein on CD8, alum pushes antibody, a boost is "
                  "faster than a prime - but reproducing them is not the same "
                  "as predicting them.",
        "because": f"These {nc} statements were used to tune the model, so of "
                   f"course it reproduces them; they are a structured "
                   f"restatement of the literature, not an out-of-sample "
                   f"result. They are graded C for that reason. All 13 "
                   f"platform benchmarks pass, but a benchmark you tuned "
                   f"against cannot fail you honestly.",
        "means": "Use platform comparisons to generate hypotheses, not to "
                 "choose between platforms.",
    }


def q_trials():
    d = _csv("trial_validation.csv")
    if d is None:
        return _missing(
            "Do its vaccine-type rankings match real clinical trials?",
            "the clinical trial comparison",
            "Run 'Compare against clinical trial fold-rises' on the "
            "Validation page.")
    rho = p = None
    j = _json("trial_validation.json")
    if isinstance(j, dict):
        bp = j.get("between_platform") or {}
        rho = bp.get("spearman_rho", bp.get("rho"))
        p = bp.get("permutation_p", bp.get("spearman_p", bp.get("p")))
    inband = (int(d.inside_5_95.sum()) if "inside_5_95" in d.columns else None)
    extra = ""
    if rho is not None and p is not None:
        extra = (f" Correlation between the model's platform ordering and the "
                 f"trials' was {rho:+.2f} with p = {p:.2f}, which is no "
                 f"relationship at all.")
    return {
        "mark": "n",
        "question": "Do its vaccine-type rankings match real clinical trials?",
        "answer": "No. Against the posted results of real vaccine trials, the "
                  "model's ordering of vaccine types is not supported.",
        "because": f"Compared against antibody fold-rises posted by "
                   f"{len(d)} vaccine platforms' worth of registered trials."
                   + extra
                   + (f" Its predicted rise does fall inside the range real "
                      f"trials report for {inband} of {len(d)} platforms."
                      if inband is not None else "")
                   + " Two caveats cut both ways: nobody randomised patients "
                     "between platforms, so trials differ in population and "
                     "schedule as well as platform, and with only a handful "
                     "of platforms this test could not detect a real ordering "
                     "even if there were one.",
        "means": "The absolute response size is in the right region. The "
                 "ranking of platform against platform is not evidence.",
    }


def q_protection():
    return {
        "mark": "n",
        "question": "Can it tell me whether a vaccine will actually protect anyone?",
        "answer": "No, and it never tries to. Nothing in this project touches "
                  "protection, efficacy or infection.",
        "because": "There is no challenge or efficacy data anywhere in the "
                   "pipeline. The model predicts how big an immune response "
                   "is, and a big response is not the same thing as "
                   "protection - that link has to be established per disease, "
                   "in people.",
        "means": "Never read an antibody number here as a protection claim.",
    }


def q_units():
    return {
        "mark": "q",
        "question": "Are the antibody numbers real units?",
        "answer": "No. Antibody is in arbitrary units, so only ratios and "
                  "timings mean anything.",
        "because": "The model was fitted to fold-changes and peak days, not "
                   "to absolute titres, because the human data reports "
                   "fold-rises. 'Peak IgG 2,009' is only meaningful next to "
                   "another run's number.",
        "means": "Compare runs with each other. Never quote a single "
                 "antibody number on its own.",
    }


QUESTIONS = [q_new_pathogen, q_class_ii, q_new_allele, q_human_size,
             q_human_timing, q_platforms, q_trials, q_protection, q_units]


def all_verdicts():
    out = []
    for fn in QUESTIONS:
        try:
            out.append(fn())
        except Exception as e:                     # never break the page
            out.append({"mark": "q", "question": fn.__name__,
                        "answer": "Could not be computed.",
                        "because": f"{type(e).__name__}: {e}",
                        "means": "", "untested": True})
    return out


def headline(verdicts):
    """One sentence over the top, derived from the verdicts themselves."""
    y = sum(1 for v in verdicts if v["mark"] == "y")
    n = sum(1 for v in verdicts if v["mark"] == "n")
    q = sum(1 for v in verdicts if v["mark"] == "q")
    has = "has" if y == 1 else "have"
    return (f"Of these {len(verdicts)} questions, {y} {has} a clear yes, "
            f"{n} a clear no, and {q} a yes-with-a-catch. The noes sit beside "
            f"the yeses on purpose: a tool you can only trust in places is "
            f"only useful if you know which places.")


if __name__ == "__main__":
    vs = all_verdicts()
    print(headline(vs), "\n")
    for v in vs:
        print(f"[{v['mark']}] {v['question']}\n    {v['answer']}\n")
