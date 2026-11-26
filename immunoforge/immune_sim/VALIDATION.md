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
[`immunoforge`](../immunoforge) instead where possible.

### Pan-allele core-scan on an unseen pathogen

`external_validation_pan.py` repeats the leave-one-pathogen-out test for
`model_pan.py` — one shared model across all alleles, a sliding 9-mer binding
register chosen by iterative realignment, and a hierarchical allele encoding
(locus / first-field group / exact allele) so an untrained allele falls back
rather than being silently scored as a different one.

16 pathogens, 8 per class, each removed from training entirely.

| | Class I | Class II |
|---|---|---|
| Pathogens | 8 | 8 |
| Usable rows | 36,792 | 18,274 |
| Median ROC-AUC, pan | **0.873** (0.849–0.961) | **0.748** (0.596–0.803) |
| Median ROC-AUC, `model.py` on the shared rows | 0.766 | 0.754 |
| Median difference | **+0.123** | **+0.013** |
| Pan better | 8/8 | 5/8 |
| Significant after BH (paired DeLong) | **7/8** | 3/8 |
| Beats composition baseline after BH | **8/8** | **2/8** |

Per pathogen:

| Pathogen | Class | n | pan | `model.py` | diff | p |
|---|---|---|---|---|---|---|
| Vaccinia virus WR | I | 7,442 | 0.849 | 0.772 | +0.112 | <1e-4 |
| SARS coronavirus Tor2 | I | 2,437 | 0.857 | 0.781 | +0.098 | <1e-4 |
| Homo sapiens | I | 1,888 | 0.861 | 0.764 | +0.112 | <1e-4 |
| Zaire ebolavirus | I | 1,573 | 0.886 | 0.751 | +0.135 | <1e-4 |
| Influenza A virus | I | 1,488 | 0.959 | 0.769 | +0.202 | <1e-4 |
| Giardia lamblia | I | 1,040 | 0.961 | 0.769 | +0.191 | <1e-4 |
| Mycobacterium tuberculosis | I | 812 | 0.952 | 0.662 | +0.285 | <1e-4 |
| Dengue virus 2 | I | 622 | 0.849 | 0.738 | +0.052 | 0.142 |
| Phleum pratense | II | 2,571 | 0.737 | 0.745 | −0.009 | 0.377 |
| Homo sapiens | II | 2,044 | 0.725 | 0.721 | +0.054 | 4.5e-4 |
| Mammarenavirus juninense | II | 868 | 0.803 | 0.792 | +0.011 | 0.516 |
| Mammarenavirus brazilense | II | 840 | 0.788 | 0.764 | +0.024 | 0.170 |
| LCM virus | II | 817 | 0.754 | 0.781 | −0.027 | 0.118 |
| Mammarenavirus guanaritoense | II | 817 | 0.747 | 0.786 | −0.040 | 0.028 |
| Mammarenavirus machupoense | II | 789 | 0.748 | 0.732 | +0.016 | 0.393 |
| SARS coronavirus Tor2 | II | 788 | 0.596 | 0.470 | +0.126 | <1e-4 |

**Class I: the added complexity is worth it.** Median AUC on a pathogen the
model has never seen rises from 0.766 to 0.873, and the gain is significant in
7 of 8 after correction. The one non-significant case (Dengue 2) still favours
the pan model.

**Class II: it is not.** The median difference is +0.013 and the pan model is
*worse* on 3 of 8 pathogens, one of them significantly so. More tellingly it
beats a plain amino-acid-composition baseline in only **2 of 8** — the same
result as the simpler model, and the same conclusion the independent
mass-spectrometry set reached. Two independent lines of evidence now say the
class II layer on an unseen pathogen is close to composition-matching. It
should not be used to make decisions on its own.

The SARS class II row is worth reading carefully rather than as a win: the
pan model's 0.596 beats `model.py`'s 0.470 (below chance) by a wide margin and
with a tiny p-value, but 0.596 is itself barely better than guessing.

### Boost kinetics against measured human data

`wetlab_validation.py` tests the dynamic model against **CMI-PB** — an openly
published NIH dataset of 166 adults given a Tdap boost, with antibody titers
sampled at days 0, 1, 3, 7, 14 and 30.

Harvest it with `python3 harvest_cmipb.py` (open PostgREST API, no key).

Subjects are split in half by ID. Calibration sees only the training half;
`--split test` evaluates on subjects never used.

| Outcome | Literature-calibrated | Wet-lab, single-config | Wet-lab + cross-config | **+ timing targets (shipped)** |
|---|---|---|---|---|
| Fold-rise inside observed subject band | **1/6** | 6/6 | 4/6 | **6/6** |
| Time-to-peak inside band | 6/6 | 6/6 | 6/6 | **6/6** |
| Behaves correctly in other configurations | — | **no** | yes | **yes** |
| Platform benchmarks passed | — | — | 12/13 | **13/13** |

The shipped column adds two cross-configuration constraints on mRNA response
*timing* (`mrna_primary_peak_day` 8–35 d, `mrna_secondary_peak_day` 3–21 d).
They were added because the long-lived plasma cell compartment, introduced to
fix antibody durability, let total IgG climb for two months and put mRNA's
primary peak at **day 61** — while protein+alum stayed correct, because it was
the only arm whose timing anything constrained. Adding the two targets raised
the training objective from 0.0002 to **1.6703**, confirming the violation was
real and not a reporting artefact; re-searching brought it to 2.5e-05 and
recovered the two antigens of fold-rise fit lost to the cross-config
constraints.

Held-out numbers under the shipped parameters: model fold-rise **3.86×**
(observed medians 1.70–10.56× across the six antigens), model peak day
**19.3 d**.

**Two caveats the "6/6 inside the band" line hides.** The subject band is the
5–95% range across individuals and it is very wide — 0.95–34× for fold-rise,
1–31 d for peak day — so sitting inside it is a weak test. Measured against
the bootstrap CI of the *median* instead:

- fold-rise: the model's 3.86× is inside the median CI for 3 of 6 antigens.
  It sits at the 87th and 94th percentile for DT and TT, so it **over-predicts
  the weakly-responding antigens**.
- peak day: the model's 19.3 d is **outside** the median CI (14–14 d) for
  **all six**. CMI-PB samples on days 0/1/3/7/14/30, so an observed peak can
  only land on one of those days and the median is pinned at 14 with no
  spread; the model is ~5 days late against a grid that cannot resolve 5 days.
  This is a real discrepancy, partially masked by the sampling schedule.

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
| **In the training set** | 9 | **0.943** | **0.835** |
| First-field group seen only | 5 | 0.735 | 0.143 |
| Locus seen only | 9 | 0.704 | 0.169 |
| Nothing seen | 9 | 0.682 | 0.125 |

| Class II alleles | n | ROC-AUC | PPV@1% |
|---|---|---|---|
| In the training set | 2 | 0.637 | 0.140 |
| First-field group seen only | 1 | 0.605 | 0.167 |
| Locus seen only | 2 | 0.606 | 0.132 |
| Nothing seen | 6 | 0.534 | 0.061 |

**On a trained class I allele the model is strong**: 84% of its top 1% of
peptides are real eluted ligands, at a prevalence of 5%. **On an allele it has
not been trained on it is not usable** for ranking, whatever the ROC-AUC
suggests. **Class II is weak even on a trained allele** — 0.637 ROC-AUC and
14% precision in the top 1%, against a 5% base rate.

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
| Random split | 0.916 → **0.921** | 0.782 → **0.765** |
| **Eluted ligands (paired)** | **+0.000**, better in 3/9 | **−0.041**, better in 0/2 |
| **Leave-one-pathogen-out** | 0.766 → **0.873**, 8/8 | 0.754 → 0.748, 5/8 |

On the rebuilt, reproducible figures the core scan does not beat `model.py`
on a class II random split at all (0.765 against 0.782). Where it does win
decisively is class I on an unseen pathogen, +0.107. On class II it buys
**nothing measurable on the independent assay** — on class II they may
cost a little, on two alleles. Stated because it is the result, not the hoped
result: `model.py` remains a reasonable choice, and the
[`immunoforge`](../immunoforge) consensus is still the
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
| Peak activated B | 13 | **1,269** | 1e3–1e6 |
| Peak memory B | 2.5 | **1,487** | 1e2–1e5 |
| Peak plasma cells | — | **474** | 1e2–1e5 |
| IgG retained at end of follow-up | ~0 | **70%** | 5–90% |
| mRNA antibody boost (fold over baseline) | **0.99×** | **2.04×** | 2–60× |
| mRNA primary peak day | **61 d** | **within 8–35 d** | 8–35 d |

The last row is the weakness documented in the previous version of this file
as holding in only 57% of parameter draws. It was not a tuning problem. It was
the shared carrying capacity.

---

## A reproducibility failure, found by rebuilding what was already committed

Rebuilding `model_pan.py` so its report and its weight files came from the
same fit did **not** reproduce the committed figures:

| Random split, ROC-AUC | Committed | Rebuilt | Change |
|---|---|---|---|
| class I, GBM member | 0.920 | 0.917 | −0.003 |
| class I, ensemble | 0.923 | **0.921** | −0.002 |
| class II, GBM member | 0.788 | 0.770 | −0.018 |
| class II, ensemble | 0.794 | **0.765** | **−0.029** |

Every class II allele moved in the same direction (HLA-DRB1\*01:01
0.788→0.757, \*04:01 0.762→0.705, \*07:01 0.787→0.778, \*15:01 0.758→0.739),
so this is not noise in one allele. Everything downstream moved with it: on
the eluted-ligand set, class II trained-allele ROC-AUC went 0.648→0.637 and
PPV@1% 0.202→0.140.

What was checked and ruled out:

- **Not nondeterminism in the current tree.** Three separate processes gave
  bit-identical results (0.770 / 0.708 / 0.765), and forcing
  `OMP_NUM_THREADS=1` gave the same again, so it is not thread-level
  floating-point reduction order.
- **Not the data.** `data/binding_affinity.csv` is unchanged since
  2026-11-25 13:56 and both runs report the same 19,152 class II rows and the
  same 4 alleles.
- **Not the committed code.** `git diff HEAD` is empty for `model_pan.py`,
  `encoding.py`, `model.py` and `stats.py`.

The probable cause: `results/eluted_validation_onehot.csv` exists, written
17:59 that day, and **no script in the tree writes it** — it was produced by
editing the allele encoding to the flat one-hot form, running it, and
reverting. The committed `pan_report.json` was most likely written during
that window, from a tree state that no longer exists, and never regenerated
afterwards. That is consistent with the direction of the change but it is not
proven, and the honest statement is that the earlier figure cannot be
reproduced and the cause was not established with certainty.

**What this means for any future claim about class II.** The reproducible
figure is now 0.765, and the gap to the unreproducible 0.794 is 0.029. A
class II improvement smaller than that cannot be distinguished from whatever
produced this discrepancy on a single run. Class II changes must therefore be
judged on repeated runs with an interval, or on the leave-one-pathogen-out
and eluted-ligand sets, not on one random-split number.

The leave-one-pathogen-out results above are **unaffected**:
`external_validation_pan.py` retrains from scratch inside each fold and never
reads `pan_report.json` or the shipped weights.

---

## A fourth defect, this one in the test rather than the model

With the recalibrated parameters the `live_attenuated_is_most_durable`
benchmark failed. The claim was that a replicating live attenuated vaccine
gives "the longest antigen exposure", measured by days the antigen stays above
**1% of that arm's own peak**. That metric normalises away the thing being
compared: an arm with a 10× higher peak is scored against a 10× higher bar, so
it measures peak shape, not persistence. Live attenuated read 15.7 d against
protein+alum's 19.3 d purely because its peak is 6.4e3 against 6.3e2.

Rescoring on a shared absolute threshold shows the original claim was simply
wrong about the tail:

| Threshold | Live attenuated | Protein + alum |
|---|---|---|
| above 1% of the larger peak (64 units) | **15.7 d** | 11.7 d |
| above 1% of the smaller peak (6.3 units) | 19.3 d | 19.3 d |
| above 1 unit | 22.0 d | **25.0 d** |

An alum depot sustains trace antigen at least as long as replication does.
What live attenuated dominates is the *amount* of antigen over time
(cumulative exposure 44,395 vs 3,489, 12.7×) and the memory it leaves (CTL
memory 450,007 vs 3,984, 113×). The benchmark now asserts those two and
states in its own text that the tail is comparable and not claimed.

This is recorded because it is the same failure mode as the single-config fit
and the unconstrained mRNA timing: a metric that cannot distinguish the
quantity of interest will pass for the wrong reason until something changes
underneath it. The benchmark had been passing on a scale artefact.

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
pathogen across 16 leave-one-pathogen-out folds. Class I is solid and the
pan-allele core-scan model is a real improvement: median AUC **0.873** against
0.766 for the simpler per-allele model, significant in 7 of 8 pathogens after
correction, beating a composition baseline in 8 of 8. Class II is weak and the
extra complexity does not help it: median AUC 0.748 against 0.754, better in
only 5 of 8, and beating a composition baseline in only 2 of 8. On the
independent mass-spectrometry set the same comparison came out +0.000 for
class I and **−0.048** for class II. Zero-shot transfer to an allele never
seen in training remains unsolved (PPV@1% 0.087 → 0.121, Wilcoxon p = 0.45).

The **boost kinetics** of the dynamic model are validated out-of-sample
against 166 real human subjects, for one formulation (Tdap, protein + alum),
in ratio and timing but not absolute magnitude. Fold-rise and time-to-peak
both sit inside the observed subject band for 6/6 antigens, but the band is
wide: against the bootstrap CI of the median the model is inside for 3/6 on
fold-rise and **0/6 on peak day**, running ~5 days late against a sampling
grid too coarse to resolve 5 days. All 13 platform benchmarks pass, one of
them only after a scale artefact in the benchmark itself was fixed. The
221-trial clinical dataset does **not** corroborate the platform ordering
(Spearman ρ = −0.10, p = 0.95, 5 platforms, not randomised).

Everything about **platform comparison** — mRNA vs protein, circRNA vs mRNA,
adjuvant choice — is a structured prior reproducing literature directions, not
a validated prediction. Treat it as hypothesis generation.
