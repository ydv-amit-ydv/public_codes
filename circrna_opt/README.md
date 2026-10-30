# circopt — dual-species, structure-aware CDS optimiser for circRNA vaccines

Redesigns the **CDS of a circular RNA** using synonymous codon changes only, so that it
translates well in **both mouse and human** (one species primary, the other held above a floor)
while keeping the molecule's secondary structure clean — in particular free of the long
double-stranded stretches that PKR/MDA5 sense and that block scanning.

Everything functional you paste in (IRES, UTRs, spacers) is **protected** and never edited.
The protein sequence is never changed: that is checked, not assumed.

## Why it is built the way it is

* **Circular folding is the primary check.** Every fold in the loop runs ViennaRNA in `--circ`
  mode; the linear fold is only ever a secondary comparison. Closing the circle changes the
  energy — on the example construct, −38.2 kcal/mol circular vs −44.0 linear.
* **"Minimise dsRNA" is a parameter, not a belief.** Total base-pairing is only a soft term.
  The hard driver is `dsrna_min_bp` (default 30 bp): nucleotides sitting in *contiguous*
  duplex, which is what has an actual innate-immunity mechanism behind it. If you want to test
  the opposite assumption — that more structure means a longer-lived, better-expressed molecule
  (the LinearDesign/DERNA direction) — set `structure_goal: "minimize_mfe"` and compare.
* **Codon fitness saturates.** Once the weighted geometric-mean CAI reaches `cai_target`
  (default 0.92) the codon term stops pushing. Without this the optimiser trades a lot of
  base-pairing for the last few points of CAI, because the highest-CAI codons in both species
  are GC-rich. `cai_floor` stays a hard constraint on the *weaker* species. Set
  `cai_target: 1.0` for plain maximisation.
* **Two-tier folding is what makes long molecules affordable.** The GA scores thousands of
  candidates with a fast circular MFE fold plus an RNAplfold accessibility screen around the
  start codon; only the handful of elites per cycle get an exact circular partition function.
* **The fast screen optimises what the run is judged on.** An earlier version scored the GA on
  a cheap objective and ranked on the exact one, so the GA never saw start-codon accessibility
  and could not fix it. The screen now includes it, benchmarked against the exact PF.

## Install

See [INSTALL.md](INSTALL.md). Short version:

```bash
conda create -y -n circopt python=3.11 && conda activate circopt
conda install -y -c conda-forge -c bioconda derna       # optional, see below
pip install -r requirements.txt
python -m circopt.cli selftest
```

## Use

```bash
# CLI
python -m circopt.cli run examples/example_params.json --out runs/demo
python -m circopt.cli run my_params.json --out runs/mouse_first --set primary_species=mouse --set rounds=6
python -m circopt.cli params > my_params.json          # dump every default to edit

# Web runner
python -m circopt.app                                   # http://127.0.0.1:5000
```

`examples/make_example.py` builds a deliberately bad demo construct (rare codons, a planted
34 bp synonymous inverted repeat, a pairable start context) so you can see the optimiser work.

## What a cycle does

```
cycle 1..rounds
  ├─ GA over synonymous codons           fast circular MFE fold + plfold start accessibility
  ├─ exact circular MFE + partition fn   on the cycle's elites only
  ├─ Boltzmann-sample the circular       <- the per-cycle "ensemble scan"
  │  ensemble of the cycle's best
  └─ re-balance weights from what the elites still fail on, then seed the next cycle
final
  ├─ re-score every elite on YOUR weights (the adaptive weights only steer the search,
  │  so totals from different cycles stay comparable)
  ├─ greedy single-codon polish on the exact objective
  └─ noLP/noGU stress test, G-quadruplex check, linear comparison, ranking
```

## Outputs (in `--out`)

| file | contents |
|---|---|
| `ranked_variants.csv` | every metric for the original + top N, both species, all penalty terms |
| `comparison_original_vs_best.csv` | original → best with deltas |
| `variants.fasta` | the sequences, headers carrying the key metrics |
| `best_structure.txt` | sequence, MFE dot-bracket, centroid |
| `ensemble_by_cycle.txt` | sampled circular structures per cycle + for the final best |
| `run.json` | full record: params, per-cycle log, weights, Pareto set |
| `plots/` | bp-probability map, mountain plot, entropy/accessibility tracks, convergence, Pareto |

## Measured performance

Timings from this build (ViennaRNA 2.7.2, 4 cores for the measurement, projections scaled by core
count). The exact circular partition function dominates and is the reason for the two-tier design.

| molecule | fast candidate eval | exact circular MFE + PF | peak RSS per exact fold |
|---|---|---|---|
| 618 nt | ~0.3 s | ~8 s | — |
| 2873 nt | ~3.3 s | ~93 s (66 s of it the PF) | 243 MB |

A 2873 nt molecule at the default search size (population 40 x 12 generations x 4 cycles)
projects to **~68 min on 4 cores, ~25 min on 24 cores**. Memory is not a constraint: only
`elite_k + 1` exact folds run at once, so ~1.7 GB at 2873 nt.

For a ~3000 nt CDS with a 1-3 h budget, spend it on search rather than on bigger folds:

```json
{ "rounds": 6, "generations_per_round": 20, "population": 64,
  "elite_k": 8, "top_n": 8, "polish_codons": 10,
  "inner_span": 250, "access_window": 200, "access_span": 150 }
```

## Key parameters

| group | parameters |
|---|---|
| species | `primary_species`, `secondary_species`, `primary_weight`, `cai_floor`, `cai_target`, `codon_metric` (`cai`/`tai`), `trna_copy_file`, `min_codon_freq` |
| structure | `structure_goal`, `mfe_per_nt_target`, `dsrna_min_bp`, `stem_target`, `stem_cap`, `long_range_nt`, `min_paired_fraction`, `red_prob` |
| initiation | `start_unpaired_min_prob` (AUG), `start_context_mean_min`, `start_upstream`, `start_downstream`, `early_orf_nt` |
| composition | `gc_min`, `gc_max`, `gc_window`, `gc_window_min`, `gc_window_max`, `max_poly_u`, `max_homopolymer`, `ir_kmer` |
| forbidden | `avoid_re_sites`, `re_sites`, `avoid_slippery`, `avoid_cryptic_splice`, `avoid_g4`, `extra_forbidden` (IUPAC) |
| folding | `temperature`, `dangles`, `no_lonely_pairs`, `no_gu`, `no_gu_closure`, `gquad`, `inner_span`, `access_window`, `access_span`, `stress_test`, `linear_compare` |
| search | `rounds`, `generations_per_round`, `population`, `elite_k`, `top_n`, `polish_codons`, `mutation_rate`, `seed`, `workers` |
| DERNA | `use_derna`, `derna_bin`, `derna_lambdas`, `derna_max_aa`, `derna_timeout_s` |

All hard constraints are reported per variant in the `violations` column rather than silently
dropping candidates, so you can see *why* nothing feasible was found if that happens.

## DERNA

[DERNA](https://doi.org/10.1089/cmb.2023.0283) (Gu, Qi & El-Kebir 2024) is integrated as a
**seed and comparison arm**, not as the optimiser, for three measured reasons:

* it optimises toward *minimal* MFE (maximally structured) — the opposite direction to the
  default here, which is exactly why it is useful as a contrast;
* its `-c` custom codon table is unusable: six plausible CSV layouts all abort with
  `std::invalid_argument`, and the source is not mirrored on GitHub, so it runs with its
  built-in human table;
* it is O(n³) and single-threaded: ~27 s for a 70 aa protein, and its λ sweep takes minutes.
  `derna_max_aa` (default 400) stops it being switched on accidentally for a long CDS.

With `use_derna: true` its designs are grafted onto designable slots only, seeded into the
population, and recorded in `run.json`.

## A tension you should expect in the results

Minimising base-pairing and wanting a *confident* fold (low ensemble diversity, low positional
entropy, centroid close to the MFE structure) pull against each other, and the brief asks for
both. A weakly paired molecule has many near-degenerate structures almost by definition, so on
the example construct the optimised variant improves nearly everything while ensemble diversity
rises 125 -> 203 and mean entropy 0.58 -> 0.97. That is the expected physics, not a regression.
Decide which you care about and set `weights.entropy` / `weights.diversity` against
`weights.paired` / `weights.dsrna` accordingly; both are reported for every variant either way.

The same applies to CAI against structure: the highest-CAI codons in human and mouse are
GC-rich, so pushing CAI up costs base-pairing. `cai_target` exists to stop that trade running
further than you want, and `plots/pareto.png` shows the frontier you are choosing a point on.

## Honest limitations

* **Structure prediction is a model.** ViennaRNA's nearest-neighbour thermodynamics at 37 °C
  with no ions, no proteins and no co-transcriptional folding. A clean dot-bracket is a
  hypothesis to test, not a measured property. Check the top variants in another tool and,
  ideally, by probing (SHAPE/DMS).
* **A fully unstructured RNA is not the goal.** `min_paired_fraction` (default 0.10) is a floor
  against hydrolytic instability, and it is a crude one.
* **The IRES keeps its native fold only because it is protected.** This package does not verify
  that the IRES still folds correctly in its new context — check that separately.
* **CAI/tAI are proxies** for translation efficiency, and A/U enrichment can reduce expression
  for reasons CAI does not capture, which is why both species' CAI are tracked and floored.
* **Circ-specific immunogenicity is mostly a purification problem.** Removing linear precursors,
  5′ triphosphates and IVT dsRNA by-products dominates, and no amount of sequence design fixes
  it. Structure design only addresses the intramolecular contribution.
* **Directed repair only covers local constraints.** GC windows and forbidden motifs are fixed
  deterministically; start-codon accessibility and duplex length are left to the search, because
  they depend on the whole fold rather than on a few positions.
* The dual-species claim rests on codon-usage tables; mouse and human usage are very similar,
  so the achievable differentiation between them is genuinely small.
