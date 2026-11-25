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

## Grade A — two further datasets, neither in training

### Mass-spectrometry eluted ligands: what a cell actually presents

`harvest_eluted.py` / `validate_eluted.py`. **42,641 peptides** recovered from
the surface of live cells by immunoprecipitation and mass spectrometry, across
**42 alleles and 131 publications**, with 810,179 decoys sampled from the same
parent proteins at the same lengths.

This is independent in *assay technology*, not merely in split. The models are
trained on purified-MHC IC50 values; natural presentation additionally
requires proteasomal cleavage, TAP transport and peptide loading, none of
which an affinity assay measures.

| Class I alleles | n | ROC-AUC | PPV@1% |
|---|---|---|---|
| **In the training set** | 9 | **0.948** | **0.833** |
| First-field group seen only | 5 | 0.701 | 0.181 |
| Locus seen only | 9 | 0.709 | 0.135 |
| Nothing seen | 9 | 0.705 | 0.083 |

**On a trained allele the model is strong**: 83% of its top 1% of peptides are
real eluted ligands, at a prevalence of 5%. **On an allele it has not been
trained on it is not usable** for ranking, whatever the ROC-AUC suggests.

> Decoys are peptides mass spectrometry did not detect, not peptides measured
> as unpresented. Absolute ROC-AUC therefore depends on the decoy ratio.
> **PPV@1% is the figure to quote** — it is far less sensitive to it, and it is
> what a designer actually does: take the top few peptides and test them.

### The allele one-hot bug this dataset found

The allele had been encoded as a flat one-hot over training alleles, with an
unseen allele falling back to index 0 — so **HLA-B\*44:03 was scored as if it
were HLA-A\*01:01**. "Pan-allele" meant pooled training, not transfer.
`encoding.AlleleEncoder` replaces it with locus / first-field group / exact
allele, so B\*44:03 inherits the B\*44 group from B\*44:02 and degrades to a
group average rather than to an arbitrary allele.

Paired over the same alleles:

| | one-hot | hierarchical | improved |
|---|---|---|---|
| Trained, ROC-AUC | 0.931 | **0.942** | 9/11 |
| Trained, PPV@1% | 0.772 | **0.809** | 10/11 |
| Untrained, ROC-AUC | 0.608 | **0.679** | 20/32 |
| Untrained, PPV@1% | 0.087 | 0.121 | 16/32 |

The first three are real improvements. The fourth is not: Wilcoxon
**p = 0.45**, and 16 of 32 is a coin flip. The pathological fallback is gone,
but **zero-shot transfer to a new allele is unsolved** and needs real MHC
pseudo-sequences, which this pipeline does not have offline.

### Does the extra model complexity pay off?

On identical cases, the pan core-scan ensemble against per-allele `model.py`:

| | class I | class II |
|---|---|---|
| Random split | 0.916 → **0.923** | 0.782 → **0.794** |
| **Eluted ligands (paired)** | **+0.000**, better in 3/9 | **-0.048**, better in 0/2 |

So the core scan, pooling and all lengths buy a small gain on a random split
and **nothing measurable on the independent assay** — on class II they may
cost a little, on two alleles. Stated because it is the result, not the hoped
result: `model.py` remains a reasonable choice, and the
[`epitope_filter_engine`](../epitope_filter_engine) consensus is still the
preferred path.

### Measured vaccine fold-rise from 221 clinical trials

`harvest_trials.py` / `validate_trials.py`. **14,364 posted GMFR measurements
from 221 completed trials** on ClinicalTrials.gov. GMFR is the quantity
`calibrate.py` fits, so platform predictions face numbers for the first time.

Clustered by trial, because one trial posts a value per strain, arm and
timepoint — treating 3,500 mRNA rows from 46 trials as independent shrank the
median's CI to 3.60–3.93 and made the model look falsified. Bootstrapping over
trials gives 2.83–5.71.

| Platform | Trials | Observed median [CI] | Model | Percentile |
|---|---|---|---|---|
| protein | 58 | 8.32× [6.75–11.51] | 2.16× | 10th |
| mRNA-LNP | 46 | 4.20× [2.83–5.71] | 2.19× | 9th |
| live-attenuated | 18 | 2.42× [1.79–4.20] | 1.26× | 11th |
| inactivated | 11 | 4.75× [2.20–7.05] | 3.01× | 27th |
| saRNA | 4 | 7.75× [2.75–11.40] | 1.24× | 0th |

Two findings, one reassuring and one not.

**Reassuring**: two wholly independent datasets agree on the order of
magnitude of a vaccine boost — CMI-PB's 166 subjects give a median 2.88×, the
221 trials give 5.71×. The gap is explicable rather than contradictory: every
CMI-PB subject was primed in infancy, so 2.88× is a boost in an already-immune
adult, while the trial set is dominated by primary series in naive subjects.
The model was fitted to the CMI-PB design and sits in the trials' lower tail —
it **under-predicts a primary series**, which is the direction of the error.

**Not reassuring**: the between-platform ordering is **Spearman ρ = -0.10,
permutation p = 0.95 — no relationship at all.** The platform rankings this
model produces are *not* corroborated by the only measured data available.

They have not been tuned to fix that, for a reason that must be weighed
against the result: platform is **not randomised** across trials. A protein
row is often a conjugate vaccine in naive infants; an mRNA row is usually a
COVID-19 vaccine in previously exposed adults. Antigen, age, prior immunity,
schedule and assay all differ. With five platforms the test also has almost no
power. So this is not a clean refutation — but it is evidence against, and it
replaces the "100% robustness" framing that previously stood behind these
claims with nothing better than a confounded null.

**Platform comparison should be treated as hypothesis generation only.**

---

## Three structural defects the targets exposed

Each was found by adding a constraint, not by inspection, and each needed a
change to the model's structure rather than a parameter tweak. They are
recorded here because every one of them passed the previous benchmark suite.

### 1. One carrying capacity for all lymphocytes

`crowd = 1 - total_lymphocytes / k_carry` throttled every proliferating
population together. mRNA-LNP drives a very large CD8 response — 1.3e7 CTL
against a `k_carry` of 1e7 — which drove the term to zero and **switched off B
cell proliferation entirely**. The model produced an mRNA antibody boost of
**0.99×**: no boost at all, flatly contrary to trial data, while passing all
ten benchmarks then in the suite.

T zone and germinal centre are distinct resources; a germinal centre does not
run out of room because the paracortex is full. Capacity is now per lineage
(`f_carry_t`, `f_carry_b`) with a weak shared ceiling.

### 2. B cell proliferation could not exceed B cell death

Proliferation was the **product** of four independent saturating gates (T
help, IL-2, depot antigen, crowding), each below 1, so the effective rate
never beat death plus plasma differentiation. The B pool could not expand at
all: peak activated B was **13 cells** and peak memory B **2.5**, against
100,000 CD4 memory in the same run. The antibody curve still had the right
shape, because a handful of plasma cells with a fitted secretion rate of 22–60
units per cell per day can produce any titre you like.

A germinal centre is a burst of 6–10 divisions over a week, so the gates are
now combined as a **geometric mean**: each can still throttle the reaction
without three of them multiplying the rate to nothing.

**This is why the cell-count targets exist.** Antibody is in arbitrary units,
so a model scored only on titre and timing has a free parameter to hide in.
Four bands close it:

```
peak_B_active      1e3 - 1e6      peak_plasma_total   1e2 - 1e5
peak_B_memory      1e2 - 1e5      igg_durability      0.05 - 0.90
```

### 3. No long-lived plasma cells, so no durable antibody

With a single short-lived plasma pool, modelled IgG always decayed back to
baseline — so no durability statement could be made at all, and the model
could not represent the single most important property of a vaccine. A
bone-marrow long-lived plasma cell compartment was added, exported from mature
germinal centres in proportion to achieved affinity and T help.

### What the fixes bought

| Quantity | Before | After | Target |
|---|---|---|---|
| Peak activated B | 13 | **1,499** | 1e3–1e6 |
| Peak memory B | 2.5 | **5,044** | 1e2–1e5 |
| Peak plasma cells | — | **1,365** | 1e2–1e5 |
| IgG retained at end of follow-up | ~0 | **55%** | 5–90% |
| mRNA antibody boost (fold over baseline) | **0.99×** | **5.70×** | 2–60× |
| mRNA secondary ÷ primary peak | **0.28×** | **4.57×** | ≥1× |

The last row is the weakness documented in the previous version of this file
as holding in only 57% of parameter draws. It was not a tuning problem. It was
the shared carrying capacity.

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
