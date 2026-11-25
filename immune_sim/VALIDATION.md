# Validation

What is validated, against what data, and how confidently. Read the grades
before the numbers.

## Evidence grades

| Grade | Meaning |
|---|---|
| **A** | Measured wet-lab data, **out-of-sample**. A held-out split or an unseen pathogen. The number is a prediction that could have failed. |
| **B** | Measured wet-lab data, **in-sample**. Real data, but the model saw it while fitting. Shows the model *can* fit, not that it generalises. |
| **C** | **Calibrated to published qualitative findings.** No measured dataset. Tuned until it reproduced the *direction* of literature results. A structured prior. |
| **D** | **Assumption.** Set by judgement. Falsifiable in principle, untested here. |

```bash
python3 scorecard.py --robustness 100
```

---

## Grade A — tested out-of-sample

### Binding prediction on an unseen pathogen

`external_validation.py` retrains with one organism **entirely removed**, then
predicts its epitopes. This is the real use case: a novel pathogen arrives and
must be ranked by a model that has never seen it.

| | Class I | Class II |
|---|---|---|
| Pathogens | 7 | 6 |
| Test peptides | 10,399 | 6,114 |
| Median ROC-AUC | **0.769** (0.662–0.855) | **0.733** (0.470–0.792) |
| Beats amino-acid-composition baseline (BH-corrected) | **7/7** | **2/6** |
| Above chance (BH-corrected) | 7/7 | 4/6 |
| Median ECE | 0.111 | 0.148 |

**This is the number to quote, not 0.916.** The random-split figure in the
main README (0.916 class I) is optimistic: peptides from the same organism are
compositionally similar and often near-duplicates, so a random split leaks.
Honest generalisation to a new pathogen is **~0.77**.

**Class II barely works on a novel pathogen.** It beats a plain amino-acid
composition baseline in only 2 of 6 pathogens, and one case (Influenza A,
n=206) came out at 0.470 — below chance. Class II calls on an unseen pathogen
should be treated as weak evidence. Use the
[`epitope_filter_engine`](../epitope_filter_engine) instead where possible.

### Boost kinetics against measured human data

`wetlab_validation.py` tests the dynamic model against **CMI-PB** — an openly
published NIH dataset of 166 adults given a Tdap boost, with antibody titers
sampled at days 0, 1, 3, 7, 14 and 30.

Harvest it with `python3 harvest_cmipb.py` (open PostgREST API, no key).

Subjects are split in half by ID. Calibration sees only the training half;
`--split test` evaluates on subjects never used.

| Outcome | Literature-calibrated | Wet-lab, single-config | Wet-lab + cross-config |
|---|---|---|---|
| Fold-rise inside observed subject band | **1/6** | 6/6 | **4/6** |
| Time-to-peak inside band | 6/6 | 6/6 | **6/6** |
| Behaves correctly in other configurations | — | **no** | **yes** |

The middle column is a warning, not an achievement. Fitting the single CMI-PB
configuration scored 6/6 while producing parameters under which an
unadjuvanted vaccine was 800x *stronger* than an adjuvanted one. Adding
cross-configuration constraints costs two antigens of single-config fit and
buys a model that is not nonsense elsewhere. The right column is the one
shipped.

The first run was a genuine out-of-sample **failure**: the literature-tuned
model predicted a 71× boost against an observed median of 1.7–11.4×, an
overestimate of roughly 10×. Recalibrating against the measured training half
fixed it, and the fix holds on the held-out half.

**Negative control.** Antigens *not* in the vaccine (OVA, Bet v 1, Fel d 1)
rise 1.18× while vaccine antigens rise 2.91× — Mann-Whitney p < 1e-4,
rank-biserial r = 0.744. The signal is real, not an artefact of the analysis.

#### A measurement trap worth knowing

"Boost ratio" is ambiguous and the two readings differ by an order of
magnitude:

- peak-after-boost ÷ **peak-after-prime** ≈ 5.9×
- peak-after-boost ÷ **pre-boost baseline** ≈ 71× (before recalibration)

Wet-lab studies report the second. Calibrating against the first and comparing
against the second is how the 10× error entered. `calibrate.py` now measures
`boost_fold_over_baseline` explicitly.

---

## Grade C — direction only, with robustness

No measured dataset in this pipeline compares vaccine platforms. These claims
reproduce documented directions (`test_platforms.py`), and `scorecard.py
--robustness` resamples the fitted parameters by ±25% to ask **how
load-bearing the parameter guesses are**.

> Robustness is **not** a probability that the biology is true. It is the
> share of plausible parameter sets giving the stated direction. It answers
> "is this conclusion an artefact of my tuning?", not "is this real".

This check earned its place immediately: it flagged two claims at 19% and 3%
robustness, which turned out to be a **real bug** — calibrating on a single
configuration (protein + alum) produced parameters under which an
*unadjuvanted* vaccine was 800× stronger than an adjuvanted one, because
nothing in the target set ever simulated the alternatives. `calibrate.py` now
carries cross-configuration constraints:

```
adjuvant_effect        alum IgG / no-adjuvant IgG   in [5, 1000]
boost_over_prime_peak  secondary peak / primary     in [1.2, 20]
mrna_ctl_over_protein  mRNA CTL / protein CTL       in [5, 1e4]
```

After the fix, robustness over 100 parameter draws (+-25%):

| Claim | Before | After |
|---|---|---|
| mRNA induces more CD8 than adjuvanted protein | 100% | **100%** |
| AS01 induces more CTL than alum | 99% | **100%** |
| An adjuvant is required for a protein antibody response | 19% | **100%** |
| The secondary response exceeds the primary *in absolute peak* | 3% | **57%** |

### The one claim that remains weak

**Secondary peak > primary peak holds in only 57% of parameter draws.** For
strongly priming platforms such as mRNA the model's primary peak can exceed
its secondary (40,612 vs 11,578), which is not what mRNA vaccine trials show.

The related quantity that *is* validated is the **fold-rise over the pre-boost
baseline** - which is what wet-lab studies actually report, and what
`wetlab_validation.py` tests against CMI-PB. Timing is also right
(21 d -> 13 d). So the model gets the shape and the clinically reported
statistic right, and understates absolute secondary magnitude. Stated in the
benchmark's own docstring rather than tuned away.

---

## Grade D — explicitly not claimed

- **Absolute titer in a named individual.** The model's antibody is in
  arbitrary units. Only ratios and timings are comparable to data.
- **Primary response kinetics.** Every CMI-PB subject was primed in infancy,
  so every time course is a boost. The boost arm is validated; the primary arm
  is not.
- **Platform constants.** Half-lives and doses are order-of-magnitude
  literature values; routing and bias are relative scales. Not fitted.
- **Protection or efficacy.** No challenge or efficacy data anywhere in this
  pipeline. Immunogenicity is not efficacy.

---

## Statistical methods

`stats.py` — all in one place so every number arrives with an interval:

- **DeLong** variance for ROC AUC: analytic CIs, and paired tests between a
  model and a baseline on the same cases.
- **Percentile bootstrap** for medians, fold-rises and metric differences.
- **Permutation tests** against chance by label shuffling.
- **Benjamini–Hochberg** FDR correction across alleles and pathogens — with 13
  alleles tested, uncorrected p-values would be misleading.
- **Brier score** with Murphy's reliability/resolution decomposition, and
  **expected calibration error**, because a ranking that is well ordered can
  still be badly calibrated as a probability.
- **Mann-Whitney U** with rank-biserial effect size for vaccine vs control.

## Reproducing

```bash
python3 harvest_iedb.py && python3 harvest_cmipb.py
python3 model.py && python3 immunogenicity.py
python3 external_validation.py
python3 calibrate.py --wetlab --iters 1200
python3 wetlab_validation.py --split test --params models/calibrated_wetlab.json
python3 scorecard.py --robustness 100
```

## The honest summary

The **epitope prediction layer** is validated out-of-sample on a novel
pathogen: solid for class I (~0.77), weak for class II.

The **boost kinetics** of the dynamic model are validated out-of-sample
against 166 real human subjects, for one formulation (Tdap, protein + alum),
in ratio and timing but not absolute magnitude.

Everything about **platform comparison** — mRNA vs protein, circRNA vs mRNA,
adjuvant choice — is a structured prior reproducing literature directions, not
a validated prediction. Treat it as hypothesis generation.
