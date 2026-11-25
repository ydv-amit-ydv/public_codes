# immune_sim

An agent-style dynamic model of the immune response, in the spirit of
C-ImmSim, driven by **predicted epitopes** and **vaccine platform physiology**.

Give it an antigen and a formulation — mRNA-LNP, circRNA, protein + AS01,
adenovirus — and it simulates the T cell, B cell and antibody response over
months: priming, class switching, affinity maturation, memory, and the faster
larger secondary response to a boost.

```
epitope source  ->  dynamic model  ->  T / B / antibody time course + figures
```

---

## Where epitopes come from

Two sources, in order of preference.

### 1. epitope_filter_engine (preferred)

The [`epitope_filter_engine`](../epitope_filter_engine) in this repository
fuses **eleven published predictors** — NetMHCpan 4.2, MixMHCpred, MHCflurry 2,
BigMHC EL/IM, PRIME, NetCTLpan for class I; NetMHCIIpan 4.3, MixMHC2pred,
CD4episcore, ImmuScope, TLimmuno2 for class II — with population-weighted HLA
coverage. That consensus is better than anything in this directory, so use it
when a run exists:

```bash
python3 simulate.py --engine ../epitope_filter_engine/engine_workspace/runs/<run_id> \
                    --platform mrna_lnp
```

`engine_adapter.py` reads the engine's `*_pairs.csv` / `*_peptides.csv`, keeps
its percentile ranks (a better cross-allele scale than nM), and uses
`Priority_Score` — which already fuses coverage, agreement, binding and
immunogenicity — as the immunodominance weight.

### 2. Built-in IEDB-trained model (fallback)

For when the engine has not been run. Harvests IEDB directly and trains its
own predictors.

```bash
python3 harvest_iedb.py      # 67,365 affinity + 60,820 T cell records
python3 model.py             # per-allele binding models
python3 immunogenicity.py    # T cell response models
python3 predict.py P0DTC2    # scan an antigen
```

| Model | Class I | Class II |
|---|---|---|
| Binding, random held-out peptides (`model.py`) | 0.916 | 0.782 |
| Binding, random held-out peptides (`model_pan.py`) | 0.923 | 0.794 |
| **Binding, unseen pathogen (leave-one-pathogen-out)** | **0.769** | **0.733** |
| **Eluted ligands, trained allele** — ROC-AUC / PPV@1% | **0.948 / 0.833** | 0.648 / 0.202 |
| **Eluted ligands, untrained allele** — ROC-AUC / PPV@1% | **0.705 / 0.083** | 0.547 / 0.044 |
| Immunogenicity (cross-validated) | 0.714 | 0.620 |

**Quote the second row, not the first.** A random peptide split leaks:
peptides from one organism are compositionally similar and often near
duplicates. Retraining with an organism entirely removed — the real situation
when a novel pathogen arrives — gives ~0.77 for class I.

**Class II barely generalises.** On an unseen pathogen it beats a plain
amino-acid-composition baseline in only 2 of 6 pathogens, and one case came
out below chance. Another reason the engine is the preferred path.

**Only use an allele the model was trained on.** Against 42,641
mass-spectrometry eluted ligands from 131 publications - peptides a live cell
really presented, an assay wholly independent of the IC50 data used for
training - a trained allele reaches ROC-AUC 0.948 with **83% precision in its
top 1%**. An untrained allele falls to PPV@1% 0.083. The hierarchical allele
encoding in `encoding.py` removes the worst failure mode (a flat one-hot
scored an unseen allele as if it were a different allele entirely) but does
not fix transfer: zero-shot prediction for a new allele needs real MHC
pseudo-sequences, which this pipeline does not have.

**The pan model's extra complexity does not pay off out of assay.** Against
`model.py` on identical eluted-ligand cases the difference is +0.000 ROC-AUC
for class I and -0.048 for class II. It is kept because the core scan and
all-length training are the right structure, and because it wins slightly on a
random split - not because it has been shown to be better where it counts.

See [VALIDATION.md](VALIDATION.md) for the full statistics.

---

## Vaccine platforms

The central idea of `platforms.py`: **how antigen is supplied changes the
response as much as which antigen it is.**

| Platform | Template t½ | MHC-I routing | Intrinsic innate | Note |
|---|---|---|---|---|
| `protein` | — (bolus) | 0.08 | very low | Needs adjuvant; strong Ab, weak CTL |
| `peptide` | — | 0.25 | very low | Depends entirely on adjuvant |
| `mrna_lnp` | 14 h | 0.95 | 0.55 | Endogenous translation → strong CD8 |
| `circrna_lnp` | 48 h | 0.95 | 0.40 | No free ends → long expression |
| `sarna_lnp` | 20 h + replicon | 0.95 | 0.75 | Dose-sparing |
| `adenovirus` | 30 h | 0.95 | 0.65 | Anti-vector immunity blunts homologous boost |
| `mva` | 24 h | 0.90 | 0.70 | Common heterologous boost |
| `dna` | 72 h | 0.90 | 0.30 | Durable, weakly expressed |
| `inactivated` | — | 0.12 | 0.25 | Native epitopes, class II route |
| `live_attenuated` | 36 h + replication | 1.00 | 0.80 | Most durable |

Exogenous protein reaches MHC class I only by cross-presentation
(`crosspres = 0.08`), whereas an RNA template is translated inside the cell and
enters class I directly. That single distinction is why protein + alum gives
good titers and almost no CTL, and why mRNA vaccines give strong CD8.

### Adjuvants

`none`, `alum`, `mf59`, `as03`, `as01`, `cpg1018`, `polyic`, `mpla`, `lps`,
`saponin`, `intrinsic`. Each sets danger-signal magnitude and duration, a
**Th1/Th2 bias**, and separate CTL and antibody multipliers — so alum (Th1 bias
0.25) drives antibody and little CTL, while AS01 (0.80) drives both.

```bash
python3 platforms.py                # describe every platform
python3 platforms.py mrna_lnp as01  # one combination
```

---

## Running it

```bash
# protein + alum, prime and boost 90 days apart
python3 simulate.py results/demo_epitopes.csv --platform protein --adjuvant alum \
                    --inject 2,270 --steps 700

# mRNA-LNP at 30 ug
python3 simulate.py results/demo_epitopes.csv --platform mrna_lnp --dose-ug 30

# adenovirus prime, MVA boost
python3 simulate.py results/demo_epitopes.csv --platform adenovirus --heterologous mva

python3 simulate.py --list-platforms
```

Writes a time course CSV, a six-panel response figure (antigen/antibody,
T cells, B cells, innate + cytokine, affinity maturation, immunodominance by
clone) and a prime/boost comparison.

---

## Web GUI

```bash
pip install -r requirements.txt
cd webapp && python3 app.py        # http://127.0.0.1:5000
```

Simulate, compare up to six formulations, predict epitopes from a UniProt
accession, and read the validation dashboard — all calling the same functions
as the CLI, so the GUI cannot drift from the validated pipeline. Details in
[`webapp/README.md`](webapp/README.md).

---

## Model structure

An ODE system with compartments matching C-ImmSim's: extracellular antigen,
**follicular-DC antigen depot**, **template** (mRNA/circRNA/vector),
**endogenous antigen**, danger signal, mature DC, macrophages, per-clone TH,
CTL and B cells with their memory pools, plasma cells, IgM, IgG, antibody
affinity, IL-2. One step is 8 h, as in C-ImmSim.

Three mechanisms that are easy to get wrong and are handled explicitly:

- **Density-dependent limitation** (C-ImmSim's `MaxLocalFact`). The IL-2 /
  T-help loop is positive feedback; without a carrying capacity it diverges.
- **The FDC antigen depot.** On a boost, pre-existing antibody clears free
  antigen within days. Germinal centres run off antigen retained on follicular
  dendritic cells, so a model gated on free antigen starves the secondary
  response.
- **Immune clearance of replicating templates.** saRNA and live attenuated
  amplify; CTL killing and antibody neutralisation must control them, or they
  saturate and express indefinitely.

### Calibration

Hand-tuning this many coupled rates does not converge — each fix displaces the
imbalance elsewhere. `calibrate.py` instead states nine quantitative targets
and searches:

| Target | Band | Achieved |
|---|---|---|
| Primary IgG peak | 10–21 d | 17.7 d |
| Secondary IgG peak | 4–10 d | 9.3 d |
| Secondary/primary ratio | 5–50× | 5.8× |
| Class switching (IgG:IgM shift) | 2–200× | 2.5× |
| T cell contraction by 60 d | 80–99.5% | 97.9% |
| Affinity maturation | 2–10× | 2.3× |

Final score **0.0** — all nine met.

```bash
python3 calibrate.py --iters 1500
```

---

## Validation

```bash
python3 validate.py        # ROC + PR curves, bootstrap CIs, for both predictors
python3 test_platforms.py  # 10 wet-lab benchmark comparisons
```

The benchmarks assert the **direction** of documented findings, two-sided and
with physiological bounds:

| Benchmark | Result |
|---|---|
| mRNA beats protein on CD8 | CTL 18,806 vs 3 |
| Alum is Th2, AS01 is Th1 | CTL/IgG ratio 0.0041 vs 0.0005 |
| Adjuvant is required | IgG 5,881 vs 47 |
| circRNA expresses longer | template 11 d vs 8.7 d |
| saRNA is dose-sparing | 2.7× at a tenth of the dose |
| Anti-vector blunts homologous boost | IgG 8,836 vs 6,323 |
| Boost beats prime | 10,192 → 11,294; peak 16 d → 8 d |
| Dose-response monotonic | 8,832 < 10,544 < 12,125 |
| Live attenuated most durable | CTL memory 52,152 vs 3 |
| All platforms physiological | 10/10 within bounds |

---

## Validation

Full detail in **[VALIDATION.md](VALIDATION.md)**, with evidence grades.
Headlines:

- **Epitope prediction** is tested out-of-sample on 7 unseen pathogens
  (10,399 class I test peptides): median ROC-AUC 0.769, beating a composition
  baseline in 7/7 after Benjamini-Hochberg correction.
- **Boost kinetics** are tested against **CMI-PB** — 166 real adults given a
  Tdap boost, titers at days 0/1/3/7/14/30. Subjects are split in half;
  calibration sees one half, evaluation the other.
  The first attempt **failed honestly**: the literature-tuned model predicted
  a 71x boost against an observed 1.7-11.4x. Recalibrating on measured data
  fixed it, and the fix holds on held-out subjects (6/6 antigens).
- **Negative control**: non-vaccine antigens rise 1.18x vs 2.91x for vaccine
  antigens, p<1e-4.
- **Platform comparisons are NOT validated against data** — they reproduce
  literature directions. `scorecard.py --robustness` quantifies how much they
  depend on parameter guesses.

```bash
python3 harvest_cmipb.py          # open NIH dataset, no key needed
python3 external_validation.py    # leave-one-pathogen-out
python3 wetlab_validation.py --split test --params models/calibrated_wetlab.json
python3 scorecard.py --robustness 100
```

## What this model can and cannot claim

**Can.** Rank formulations and schedules. Predict which epitopes dominate and
which HLA types respond. Show the direction and rough timing of a response.
Explore dose, interval, adjuvant choice and heterologous boosting.

**Cannot.** Predict an absolute antibody titer in a real subject. The ROC/AUC
figures validate the **epitope prediction layer only**. The dynamic model's
kinetics are calibrated against qualitative benchmarks and are **not fitted to
any measured time-course data** — none was available, because IEDB holds
epitopes and assays, not longitudinal titers. `validate.py` prints this caveat
on every run.

Grounding the kinetics empirically needs post-vaccination titer time courses,
which live in **ImmPort / ImmuneSpace** and require registration. That is the
single most valuable next step, and would turn the benchmark targets from
assumptions into fitted parameters with confidence intervals.

Platform constants in `platforms.py` are a structured, falsifiable prior:
half-lives and dose ranges are order-of-magnitude literature values, the
routing and bias scales are relative values set so the ranking reproduces
documented comparisons. They are not measured constants, and
`test_platforms.py` states exactly which comparisons each one is held to.

---

## Files

| File | Purpose |
|---|---|
| `platforms.py` | Vaccine platforms and adjuvants |
| `dynamics.py` | The ODE immune model |
| `simulate.py` | CLI runner and figures |
| `calibrate.py` | Parameter search against benchmarks |
| `engine_adapter.py` | Read epitope_filter_engine output |
| `harvest_iedb.py` | Pull real data from IEDB |
| `model.py` | Binding models |
| `immunogenicity.py` | T cell response models |
| `predict.py` | Antigen → epitopes → C-ImmSim parameters |
| `validate.py` | ROC / PR validation |
| `test_platforms.py` | Wet-lab benchmark suite |
| `encoding.py` | Core-scan peptide features (register + flanks + physicochemical) |
| `model_pan.py` | Pan-allele core-scanning ensemble (the better binding model) |
| `external_validation_pan.py` | Leave-one-pathogen-out, pan vs `model.py` |
| `harvest_cmipb.py` | Measured human titer time courses (CMI-PB) |
| `wetlab_validation.py` | Held-out subject comparison |
| `external_validation.py` | Leave-one-pathogen-out for `model.py` |
| `stats.py` | DeLong, bootstrap, permutation, BH-FDR, calibration |
| `scorecard.py` | Every claim graded A-D, with robustness sampling |
| `webapp/` | Flask GUI — simulate, compare, predict, validation dashboard |

Data and trained weights are gitignored; the scripts regenerate them.
Requires `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`.
