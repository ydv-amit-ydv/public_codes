# Epitope Filter Engine: the complete guide

What the engine does, how a result is produced, and what every term and column header means.
The README says how to install and run it; this document explains the output.

Contents

1. [What the engine does, in one page](#1-what-the-engine-does-in-one-page)
2. [How a run flows](#2-how-a-run-flows)
3. [The models: what each one measures](#3-the-models-what-each-one-measures)
4. [How the tools are combined (votes, agreement, binder, candidate)](#4-how-the-tools-are-combined)
5. [The rules and the strategies](#5-the-rules-and-the-strategies)
6. [Population coverage](#6-population-coverage)
7. [Clusters, regions and the minimum panel](#7-clusters-regions-and-the-minimum-panel)
8. [Host-proteome identity](#8-host-proteome-identity-self-similarity)
9. [Every output file and every column](#9-every-output-file-and-every-column)
10. [Glossary](#10-glossary)
11. [What the numbers do not tell you](#11-what-the-numbers-do-not-tell-you)

---

## 1. What the engine does, in one page

You give it a protein (or peptides) and a set of HLA alleles. It cuts the protein into peptides, asks many prediction programs how
well each peptide is presented by each allele and how immunogenic it looks, combines the answers, weighs them by how common each
allele is in the chosen population, and returns a ranked list of candidate T-cell epitopes, grouped into regions of the protein
and reduced to the smallest set of peptides that covers the population.

* **MHC class I** (HLA-A, -B, -C) presents peptides of 8-14 residues to **CD8 T cells** (cytotoxic, "CTL"). 
* **MHC class II** (HLA-DR, -DQ, -DP) presents peptides of 11-25 residues to **CD4 T cells** (helper, "HTL").
* A **CTL+HTL region** is a stretch of the protein where a class I and a class II candidate overlap: one piece of sequence that both
  T-cell arms can recognise. These are usually the first choice for a vaccine.
* Everything runs on your machine. The server never imports torch or TensorFlow; each tool runs as its own program.

What it does **not** do: it predicts nothing about antibodies (B cells), and it does not decide safety on its own (allergenicity,
toxicity and antigenicity filters are separate steps, see section 11).

## 2. How a run flows

```
protein / peptides
   │  cut into windows (class I: 8-14-mers, class II: 11-25-mers), every position
   ▼
quick tools score every peptide-allele pair  (MixMHCpred, NetMHCpan, MHCflurry, PRIME, IEDB immunogenicity; class II: MixMHC2pred, NetMHCIIpan)
   │  the gate: only pairs that can still pass go on  (exact for the agreement rule, nothing a full run would find is lost)
   ▼
slow tools score the survivors  (NetCTLpan, NetChop, BigMHC; class II: CD4episcore, TLimmuno2, ImmuScope)
   ▼
consensus: votes per tool → agreement → binder / candidate calls per pair
   ▼
peptide table: one row per peptide over all its alleles; population coverage; host-proteome check; priority score
   ▼
clusters (overlapping windows become one epitope) → representatives → minimum panel
   ▼
regions (class I and class II candidates that overlap) → report, CSVs, zip
```

* **Automatic funnel** (default): slow tools only see pairs that can still pass the rule.
* **All models**: every tool scores everything (slow, complete).
* **Manual funnel**: you build the stages (which tools, which pass rule) yourself.
* **Re-rank**: the predictions of a finished run are kept; choosing another strategy re-applies the rules in seconds.
* **Reuse finished results** (prediction cache): a peptide-allele pair a tool has already scored is not scored again.
* A model that fails is **skipped with a warning**; the run goes on with the others.

## 3. The models: what each one measures

Lower or higher is better as stated. "%Rank" is a **percentile against random natural peptides** for that allele: %Rank 1 means
the peptide is in the best 1 % of random peptides.

### MHC class I

| Tool | What it predicts | Scale | Better |
|---|---|---|---|
| NetMHCpan 4.2 | binding / eluted-ligand likelihood; BA = binding affinity | `NMP_%Rank_EL`, `NMP_%Rank_BA` (%), `NMP_Aff_nM` (nM) | lower |
| MixMHCpred 3.0 | presentation (ligand motif model) | `MixMHC_%Rank` | lower |
| MHCflurry 2 | binding affinity, antigen processing and presentation | `MHCF_Presentation_%Rank`, `MHCF_Affinity_nM` | lower |
| PRIME 2.1 | likelihood of **recognition by CD8 T cells** (immunogenicity) | `PRIME_%Rank` | lower |
| IEDB immunogenicity (Calis 2013) | amino-acid positions that favour immunogenicity ("positions to mask") | `IEDB_Immunogenicity`, a score | higher |
| NetCTLpan 1.1 | MHC binding + TAP transport + proteasomal cleavage, 8-11-mers only | `CTLpan_%Rank` | lower |
| NetChop 3.1 | C-terminal cleavage by the proteasome | `NetChop_Cterm`, 0-1 | higher |
| BigMHC EL | presentation (eluted-ligand), deep learning | `BigMHC_EL`, 0-1 | higher |
| BigMHC IM | immunogenicity, deep learning | `BigMHC_IM`, 0-1 | higher |

### MHC class II

| Tool | What it predicts | Scale | Better |
|---|---|---|---|
| NetMHCIIpan 4.3 (EL) | eluted-ligand presentation | `NII_Score_EL` (0-1), `NII_%Rank_EL` | score higher / %Rank lower |
| NetMHCIIpan 4.3 (BA) | binding affinity | `NII_IC50_nM`, `NII_%Rank_BA` | lower |
| MixMHC2pred 2.1 | ligand presentation | `MIX2_%Rank` | lower |
| CD4episcore (IEDB) | CD4 T-cell immunogenicity of the peptide alone (no HLA) | `CD4_Episcore` percentile | **lower** |
| TLimmuno2 | CD4 immunogenicity of a peptide-allele pair (binding network feeding an immunogenicity network); raw score | `TLimmuno2`, 0-1 | higher |
| ImmuScope-IM | CD4 immunogenicity of a peptide-allele pair, mean of 10 networks | `ImmuScope_IM`, 0-1 | higher |

TLimmuno2 and ImmuScope read at most 21 residues; longer pairs are left unscored (and the run says so).

## 4. How the tools are combined

**Vote.** Every tool gives each pair a yes/no vote:

* percentile tools (NetMHCpan, MixMHCpred, MHCflurry, NetMHCIIpan, MixMHC2pred, PRIME, NetCTLpan) vote **yes** when %Rank ≤ the cut-off
  (`bind_rank_max`, default 2 for class I, 5 for class II);
* BigMHC-EL votes yes when its score ≥ `bigmhc_el_min` (default 0.2);
* tools with no absolute scale for immunogenicity (BigMHC-IM, IEDB/Calis, CD4episcore, TLimmuno2, ImmuScope) vote **relatively**: yes when the
  pair is in the top `imm_top_q` (default 25 %) of the **binders of that allele** in this run (CD4episcore, which has no allele, ranks the peptide
  among the binder peptides). A relative vote needs at least 8 binder pairs to rank among; with fewer the tool gives no vote and the run says so.

**Agreement** = votes / tools that scored the pair. A tool that did not score a pair is not counted against it.

**Binder** (a pair): the binding tools agree as the rule says (`bind_min_agree`: *any one*, *majority* (default), *all*, or a number). 
**Strong** binder: a binder whose median binding %Rank ≤ `strong_rank_max` (0.5 class I, 2 class II).

**Immunogenic (`Imm_OK`)**: the immunogenicity tools agree as `imm_min_agree` says. **`Proc_OK`**: a processing tool called it processed (only
when `require_processing` is on).

**Pair candidate (`Pair_Candidate`)** = Binder and Imm_OK and Proc_OK and your per-tool filters.

**Candidate (a peptide)** = a peptide with at least `min_alleles` candidate alleles and `min_coverage` population coverage, and (if you chose
a removing host rule) not matching the host proteome.

**Priority score** (0-100): the peptide's percentile standing on four components, weighted (`rank_weights`, default coverage 0.35,
agreement 0.25, binding 0.20, immunogenicity 0.20). 100 is the best of this run. It is relative to the peptides of this run, not an absolute
probability.

## 5. The rules and the strategies

Strategy = a named set of rules. "Adjust the rules" overrides single fields (blank = the strategy's value).

| Rule | Meaning |
|---|---|
| `bind_rank_max` | %Rank at or below which a percentile tool votes "binder" |
| `strong_rank_max` | median %Rank at or below which a binder is "strong" |
| `bigmhc_el_min` | BigMHC-EL score at or above which it votes yes |
| `bind_min_agree` | how many binding tools must say yes: any / majority / all / number |
| `imm_top_q` | an immunogenicity tool votes yes for the top this-share of the binders |
| `imm_abs_rank_max` | cut-off for immunogenicity tools that have a percentile (PRIME) |
| `imm_min_agree` | how many immunogenicity tools must say yes: off / any / majority / all |
| `ctlpan_rank_max`, `netchop_min`, `require_processing` | processing filters (off by default) |
| `min_alleles` | a candidate must be a candidate in at least this many alleles |
| `min_coverage` | ...and reach this population coverage |
| `collapse`, `overlap_fraction` | merge overlapping windows into one epitope when they overlap by this fraction of the shorter one (default 0.7) |
| `cluster_by` | `peptide` (class I) or `core` (class II: the predicted 9-residue cores overlap) |
| `host_mode`, `host_near_min` | host-proteome check: off / flag / exclude_identical / exclude_near; the "near" threshold (default 8/9) |
| `panel_target`, `panel_cap` | minimum panel stops at this coverage (default 95 %) or this many peptides (25) |
| `rank_by` | consensus (priority score), coverage, promiscuity (number of alleles), binding, immunogenicity |
| `rank_weights` | weights of the four priority components |

Strategies: **Balanced** (default), **Strict** (every binding tool, a tight cut-off, top 10 % immunogenic, at least 2 alleles, identical host
peptides removed), **Max coverage** (broadly presented peptides, ranked by coverage), **IEDB cut-offs** (percentile ≤ 1 class I, ≤ 10 class II, any
one tool), **Wide net** (a first look). The Strategy box on the page prints what the chosen one means.

## 6. Population coverage

* Each allele has a frequency in the chosen population (default: India, six regions, weights in a CSV). 
* **Coverage of a set of alleles** = the share of people who carry at least one of them. Within a locus the frequencies add; a person has two
  haplotypes, so P(carries none) = (1 − sum)². Loci are treated as independent (an approximation; linkage between loci is ignored).
  This is computed per region and averaged with the region weights.
* `Allele_Data` says where an allele's frequency came from: `exact` (typed allele), `dqb1_component` (DQA1/DQB1 pair: the DQB1 chain stands in),
  `ld_proxy` (DRB3/4/5: estimated from the DRB1 alleles they travel with), `none` (no data: weight 0, counted in agreement but not in coverage;
  HLA-DP and mouse alleles).
* Without a population file, coverage is the plain share of the chosen alleles.

## 7. Clusters, regions and the minimum panel

* **Cluster**: overlapping candidate windows are one epitope. Class I: windows overlapping ≥ 70 % of the shorter one. Class II: the predicted
  **binding cores** overlap ≥ 70 % (so the many 15-25-mers built around one core become one cluster). The best-ranked peptide of a cluster is its
  **representative**.
* **Minimum panel**: greedy; start with nothing, add the representative that raises population coverage most, repeat until `panel_target` or
  `panel_cap`. `gain` is the coverage that peptide adds; `cumulative_coverage` the coverage so far.
* **Region**: candidate peptides of **both classes** that overlap on the protein. Type **CTL+HTL** (both), **CTL only**, **HTL only**. Region
  score = mean priority score of the best CTL and HTL peptide, +10 for CTL+HTL. Regions are ranked CTL+HTL first, then by score.

## 8. Host-proteome identity (self-similarity)

A candidate identical to a host protein may be tolerated by the immune system or risk autoimmunity, so it is a poor vaccine candidate.
Each candidate is compared with UniProt's reference proteome of the host (Settings → Host proteome → Download; human by default, mouse optional).

* The unit is the **9-mer** (the MHC-bound core; the whole peptide for 8-mers). A longer peptide is judged by its best 9-mer window.
* **`identical`**: a 9-mer of the peptide occurs exactly in a host protein. **`near`**: a 9-mer matches with one mismatch. **`distinct`**: neither.
* The search is exact for matches with up to 2 mismatches (a 9-mer with two changes keeps one of its three 3-residue blocks, and the host windows
  are looked up by block). Anything with 3 or more mismatches is reported as distinct, with a blank `Host_Best_Identity` (below 7/9).
* Only candidate peptides are checked. Settings → Host proteome has a "Check peptides now" box for any peptides you paste.
* It checks sequence, not expression: whether the host protein is made in a relevant tissue is a separate question.

## 9. Every output file and every column

Files of a run (Projects tab, or the zip): `master.csv`, `peptides.csv`, `candidates.csv`, `regions.csv`, `panel.csv`, `summary.csv` /
`pairs.csv`, per-model CSVs, `run_details.json`, and the HTML report (a single self-contained page).

### 9.1 `pairs.csv`: one row per peptide × allele, the raw tool outputs

| Column | Meaning |
|---|---|
| `Protein`, `Start` | the protein the window was cut from and its 1-based start position (peptide input mode: `peptide`) |
| `Peptide`, `Length` | the sequence and its length |
| `Allele` | the HLA allele (class II pairs `HLA-DQA1*05:01/DQB1*02:01`; mouse `H2-IAb`) |
| `NMP_Score_EL`, `NMP_%Rank_EL` | NetMHCpan eluted-ligand score and percentile |
| `NMP_Score_BA`, `NMP_%Rank_BA`, `NMP_Aff_nM` | NetMHCpan binding-affinity score, percentile and predicted IC50 in nM |
| `NMP_BindLevel` | NetMHCpan's own label: SB strong binder, WB weak binder |
| `MixMHC_Score`, `MixMHC_%Rank` | MixMHCpred score and percentile |
| `MHCF_Affinity_nM`, `MHCF_Affinity_%Rank` | MHCflurry binding affinity (nM) and its percentile |
| `MHCF_Processing` | MHCflurry's antigen-processing score |
| `MHCF_Presentation`, `MHCF_Presentation_%Rank` | MHCflurry presentation score and percentile |
| `PRIME_Score`, `PRIME_%Rank` | PRIME immunogenicity score and percentile; `PRIME_%RankBinding` is MixMHCpred's binding rank used inside PRIME |
| `IEDB_Immunogenicity` | Calis immunogenicity score (higher = more) |
| `CTLpan_MHC`, `CTLpan_TAP`, `CTLpan_Cleavage`, `CTLpan_Combined`, `CTLpan_%Rank` | NetCTLpan: MHC-binding, TAP-transport and cleavage components, the combined score and its percentile |
| `NetChop_Cterm` | probability that the proteasome cuts at the peptide's C-terminus |
| `BigMHC_EL`, `BigMHC_IM` | BigMHC presentation and immunogenicity (0-1) |
| `NII_Core` | NetMHCIIpan's predicted 9-residue binding core inside the peptide |
| `NII_Score_EL`, `NII_%Rank_EL` | NetMHCIIpan eluted-ligand score and percentile |
| `NII_IC50_nM`, `NII_%Rank_BA` | NetMHCIIpan predicted affinity (nM) and percentile |
| `MIX2_%Rank` | MixMHC2pred percentile |
| `CD4_Episcore` | CD4episcore percentile of the peptide (no allele; lower = more immunogenic) |
| `TLimmuno2`, `ImmuScope_IM` | class II pair immunogenicity scores (0-1, higher = more) |
| `Pass_<tool>` (`Pass_NMP`, `Pass_Mix`, `Pass_NII`, `Pass_CD4`, ...) | true when that tool scored the pair and met the per-tool filter you set (no filter = true if scored) |
| `Pass_All` | true when every selected tool passed. In a manual funnel it also needs `Funnel_Alive` |
| `S<k>_In`, `S<k>_Pass`, `Funnel_Reached`, `Stopped_At`, `Funnel_Alive` | manual funnel only: stage k scored it / let it through; stages that scored it; the stage that stopped it (0 = it came through); still alive at the end |

### 9.2 `master.csv`: pairs with the consensus calls (the main table)

Everything in `pairs.csv`, plus:

| Column | Meaning |
|---|---|
| `Priority_Rank` | rank of the peptide (1 = best) by the chosen `rank_by` |
| `Allele_Freq` | population frequency of the allele (national) |
| `Allele_Data` | where that frequency came from: exact / dqb1_component / ld_proxy / none |
| `Pair_Call` | `strong binder`, `binder` or blank |
| `Pair_Candidate` | the pair passes binder + immunogenicity + processing + your filters |
| `Binder`, `Strong` | the binder call and the strong-binder call (section 4) |
| `Bind_Consensus_%Rank` | median %Rank over the binding tools that scored the pair |
| `Vote_<TOOL>` | 1 = yes, 0 = no, blank = not scored / no vote (`Vote_NII_EL`, `Vote_MIX2`, `Vote_PRIME`, ...) |
| `TopPct_<TOOL>` | for relative tools: the pair's position among the binders in percent (lower = better) |
| `Plausible` | the pair is a binder, so immunogenicity tools rank it among binders |
| `Bind_Votes`, `Bind_Tools` | votes / tools in the binding group; `Imm_Votes`, `Imm_Tools`; `Proc_Votes`, `Proc_Tools` likewise |
| `Agree_Votes`, `Agree_Tools`, `Agree_Fraction` | all groups together; fraction = votes / tools |
| `Agree_List` | the tools that said yes, semicolon-separated |
| `Imm_TopPct_Consensus` | mean of the `TopPct_*` of the immunogenicity tools (lower = better) |
| `Imm_OK`, `Proc_OK`, `Filters_OK` | immunogenicity / processing / your filters satisfied |
| `Priority_Score`, `Candidate`, `Is_Representative`, `Cluster_ID`, `Region`, `In_Panel`, `Panel_Order` | copied from the peptide table |
| `Coverage_Binder`, `Coverage_Candidate`, `N_Alleles_Binder`, `Agree_Fraction_PopW`, `Bind_Consensus_PopW` | peptide-level values, see 9.3 |

### 9.3 `peptides.csv` and `candidates.csv`: one row per peptide

`candidates.csv` is the same table restricted to the cluster representatives, in rank order.

| Column | Meaning |
|---|---|
| `Priority_Rank`, `Priority_Score` | rank and 0-100 score (section 4) |
| `Protein`, `Start`, `End`, `Peptide`, `Length` | where the peptide sits (1-based, inclusive) |
| `Alleles_Tested` | how many alleles were scored for it |
| `N_Alleles_Binder`, `N_Alleles_Strong`, `N_Alleles_Candidate` | in how many alleles it is a binder / strong binder / pair candidate |
| `Alleles_Binder`, `Alleles_Candidate` | those alleles, listed |
| `Best_Allele` | the binder allele with the lowest binding %Rank |
| `Bind_Consensus_Best` | that allele's median binding %Rank |
| `Agree_Fraction_Best` | its tool agreement; `Imm_Votes_Best` / `Imm_Tools_Best` / `Proc_Votes_Best` its votes |
| `Best_Core`, `Core_Start`, `Core_End` | class II: the predicted 9-residue core in that allele and its protein positions |
| `Imm_TopPct_Best` | best (lowest) mean immunogenicity position among its binder alleles, in % |
| `Agree_Fraction_PopW` | tool agreement averaged over all tested alleles, weighted by allele frequency |
| `Bind_Consensus_PopW` | population-weighted geometric mean of the binding %Rank |
| `Coverage_Binder`, `Coverage_Candidate` | population coverage of the binder alleles / candidate alleles (also `_<region>` columns per region) |
| `Candidate` | passes the peptide-level rule (alleles, coverage, host rule) |
| `Mass_Da`, `GRAVY`, `Net_Charge_pH7.4`, `pI`, `Aromatic_Frac`, `Hydrophobic_Frac`, `N_Cys`, `N_Met`, `N_Trp` | synthesis properties: monoisotopic-style mass, Kyte-Doolittle hydropathy (GRAVY), net charge at pH 7.4, isoelectric point, aromatic and hydrophobic fractions, counts of Cys, Met, Trp |
| `Synthesis_Notes` | warnings for making the peptide: Cys (oxidation), Met/Trp, N-terminal Gln/Glu, Asp-Pro / Asp-Gly, long hydrophobic runs, GRAVY > 1 |
| `Host_Self`, `Host_Exact_Kmers`, `Host_Best_Identity`, `Host_Best_Protein`, `Host_Best_Window`, `Host_Best_Start`, `Host_Excluded` | host-proteome check, section 8 (`Host_Excluded` = removed from the candidates by your rule) |
| `Cluster_ID`, `Cluster_Size`, `Is_Representative` | the cluster it belongs to, its size, whether this peptide stands for it |
| `Region` | the protein span of the cluster, `Protein:first-last` |
| `Cluster_Representative` | the peptide that stands for this one (itself when it is the representative) |
| `Cluster_Overlap` | overlap with the representative (cores for class II) as a share of the shorter one |
| `Cluster_Members` | on a representative: the other peptides in its cluster, `SEQ@start` |
| `In_Panel`, `Panel_Order`, `Panel_Cum_Coverage` | in the minimum panel; order added; coverage reached when it was added |
| `Region_ID`, `Region_Type` | the epitope region it falls in and its type (CTL+HTL, CTL only, HTL only) |

### 9.4 `regions.csv`

| Column | Meaning |
|---|---|
| `Region_Rank`, `Region_ID` | rank (CTL+HTL first, then score) and the label R001... |
| `Type` | CTL+HTL, CTL only, HTL only |
| `Protein`, `Start`, `End`, `Length`, `Sequence` | the stretch of protein covered by the region's overlapping candidates |
| `Region_Score` | mean priority score of the best CTL and HTL peptide, +10 for CTL+HTL |
| `CTL_Peptides`, `HTL_Peptides` | how many candidate windows of each class lie in it |
| `CTL_Best_Peptide`, `CTL_Best_Allele`, `CTL_Coverage`, `CTL_Priority_Score` | the best class I peptide of the region, its best allele, its population coverage and priority score (`HTL_*` the same for class II) |
| `In_Panel` | one of the region's peptides is in a minimum panel |
| `Synthesis_Notes` | the synthesis warnings for the whole region sequence (5-60 residues) |

### 9.5 `panel.csv`

`order` (1 = added first), `peptide`, `protein`, `start`, `length`, `gain` (coverage added), `cumulative_coverage` (coverage so far),
`new_alleles` (the alleles that peptide newly covers), `region`.

### 9.6 `summary.csv` (class II pipeline) and the result page

`Alleles_tested`, `Alleles_passed` (alleles with `Pass_All`), `Alleles_passed_list`, then the best value per tool: `Best_NII_%Rank_EL`,
`Best_NII_%Rank_BA`, `Best_MIX2_%Rank`, `CD4_Episcore`, `Max_ImmuScope_IM`, `Max_TLimmuno2` (class I: `Best_NMP_%Rank_EL`, `Max_BigMHC_IM`, ...).

The page shows: **unique peptides scored**, **epitope regions**, **CTL+HTL regions**, panel coverage per arm (with the number of peptides),
**candidates like host proteins**, run time. "Chained" = the slow tools only scored the pairs the quick tools let through.

### 9.7 Single-model tab

One tool alone, everything it reports, sorted by its main score; nothing is combined or filtered. Use it to compare a number with what the
pipeline used.

### 9.8 `run_details.json`

Everything needed to reproduce a run: inputs, alleles, tools and versions, the rules in force (every field), the population used, warnings.

## 10. Glossary

| Term | Meaning |
|---|---|
| **Allele** | one version of an HLA gene, e.g. `HLA-A*02:01` (A = gene, 02 = allele group, 01 = protein). More than ~35,000 exist. |
| **Epitope** | the short peptide a T-cell receptor actually recognises. |
| **MHC / HLA** | the molecules that display peptides on cell surfaces. HLA is the human MHC. |
| **Class I / II** | class I (A, B, C) shows peptides from inside the cell to CD8 T cells; class II (DR, DQ, DP) shows peptides taken up from outside to CD4 T cells. |
| **CTL / HTL** | cytotoxic / helper T lymphocyte (CD8 / CD4). |
| **Binding core** | the 9 residues of a class II peptide that sit in the groove; the rest hangs out of both ends. |
| **%Rank / percentile rank** | position of a score among scores of random natural peptides for the same allele. Lower = better; 2 % is a common "binder" line for class I. |
| **IC50 (nM)** | the concentration that displaces half of a reference peptide; lower = tighter binding. 50 nM strong, 500 nM weak by convention. |
| **EL / BA** | eluted ligand (trained on peptides found on cells, so it learns processing and presentation) / binding affinity (trained on test-tube binding). |
| **Presentation** | the peptide is processed, transported and displayed by the MHC. |
| **Immunogenicity** | the displayed peptide actually triggers a T-cell response. Much harder to predict than binding. |
| **Processing / TAP / proteasome** | the cell chops proteins (proteasome), moves the pieces into the ER (TAP) where they can meet MHC class I. |
| **Pseudosequence** | the ~34 MHC residues that contact the peptide, used by pan-allele networks to handle alleles they were not trained on. |
| **Pan-allele** | one network for all alleles. |
| **Consensus** | the combined call of several tools. |
| **Vote / agreement** | section 4. |
| **Gate / chain / funnel** | running slow tools only on pairs that can still pass. |
| **Representative** | the best peptide of a cluster of overlapping ones. |
| **Minimum panel** | the smallest set of representatives reaching the coverage goal. |
| **Population coverage** | share of the population that carries at least one of the alleles a peptide is presented by. |
| **Hardy-Weinberg** | the genetics rule used to turn allele frequencies into the share of people carrying them. |
| **Promiscuous epitope** | one presented by many alleles. |
| **Tolerance / autoimmunity** | the immune system ignores what resembles the body's own proteins; a vaccine peptide that matches host proteins may be ignored or may provoke an attack on the host protein. |
| **Reference proteome** | UniProt's curated set of proteins of one organism, one sequence per gene. |
| **k-mer** | a stretch of k residues; a 9-mer is nine in a row. |
| **Cache** | stored tool results; the same pair is not scored twice. |
| **Watchdog** | the engine kills a tool that shows no activity for the stall time (Settings) instead of hanging. |
| **Guarded model** | a model whose failure is reported and skipped; the run continues. |

## 11. What the numbers do not tell you

* Predictions are not measurements. Immunogenicity tools in particular have modest accuracy (AUC around 0.6-0.8 depending on the data set);
  agreement across tools raises confidence but cannot confirm a response.
* Binding is necessary, not sufficient: T-cell receptor repertoire, tolerance, dose, adjuvant and delivery all matter.
* Coverage numbers assume independent loci and the frequencies in your population file; alleles with no data count in agreement but not in coverage.
* "Distinct" in the host check means no 9-mer with ≥ 7 of 9 identical residues was found; it does not exclude cross-reactivity with a
  structurally similar self peptide.
* Not part of the engine (yet): antibody (B-cell) epitopes, class II processing (MHCII-NP needs Docker), and the safety filters:
  antigenicity (VaxiJen), allergenicity, toxicity, and conservation across variants. These are meant to be added as post-filters on the
  candidate table.
