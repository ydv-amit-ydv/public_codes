"""The engine's documentation, as data, and the page built from it.

One source feeds two things: the full documentation page (served at /guide) and the short "?" help on every table header
and result tile in the web page (served at /api/docs). Edit the data here; both follow.
"""
import html

E = html.escape

REF_NOTE = ("References are those given in each tool's own documentation or publication. Where a tool states a preferred citation, "
            "use that one in a paper.")

# ------------------------------------------------------------------------------------------------ models
# key, name, full name, class, role, version, by, inputs, outputs [(column, meaning)], method, scale, use, limits, licence, refs
MODELS = [
    dict(key="netmhcpan", name="NetMHCpan", full="Network-based MHC prediction, pan-specific", cls="I", role="Binding and presentation",
         version="4.2 (method of 4.1)", by="DTU Health Tech (Nielsen group), Denmark; run here through the DTU program you installed",
         inputs="Peptides of 8-14 residues and HLA class I alleles (e.g. HLA-A*02:01). Any allele with a known protein sequence works (pan-specific).",
         outputs=[("NMP_Score_EL", "Eluted-ligand score, 0-1; higher = more likely presented"),
                  ("NMP_%Rank_EL", "Percentile of that score among random natural peptides for the allele; lower = better"),
                  ("NMP_Score_BA", "Binding-affinity score, 0-1"),
                  ("NMP_%Rank_BA", "Percentile of the affinity score"),
                  ("NMP_Aff_nM", "Predicted IC50 in nM"),
                  ("NMP_BindLevel", "SB = strong binder, WB = weak binder (by the %Rank thresholds you set in the program)")],
         method=("An ensemble of feed-forward neural networks. The input is the peptide (encoded with the BLOSUM62 matrix) and the MHC "
                 "'pseudosequence', the 34 residues of the HLA molecule that touch the peptide. The networks were trained jointly on in-vitro binding "
                 "affinities (BA) and on peptides eluted from cells and read by mass spectrometry (EL). The EL output is a likelihood-like score from 0 to 1. "
                 "The BA output is a transformed affinity, s = 1 - log(IC50) / log(50000), so IC50 = 50000^(1 - s) nM. %Rank is the share of random natural "
                 "peptides of the same length that score higher than the query for that allele."),
         use="Votes 'binder' when %Rank_EL is at or below the binding cut-off (default 2 %). Its %Rank also enters the median binding rank used for 'strong binder'.",
         limits="The scores describe presentation on the cell surface; they do not say a T cell will respond. Rare alleles with few training peptides are less reliable.",
         licence="Free for academic use from DTU; not redistributed here.",
         refs=["Reynisson B, Alvarez B, Paul S, Peters B, Nielsen M. NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. Nucleic Acids Res 2020;48:W449-W454."]),
    dict(key="mixmhcpred", name="MixMHCpred", full="Mixture MHC predictor", cls="I", role="Presentation",
         version="3.0", by="Gfeller laboratory, University of Lausanne",
         inputs="Peptides of 8-14 residues and HLA class I alleles.",
         outputs=[("MixMHC_Score", "Raw score (log-likelihood-like)"), ("MixMHC_%Rank", "Percentile among random peptides for the allele; lower = better")],
         method=("Each allele is described by a position weight matrix (a table of how much each amino acid is favoured at each peptide position) plus a "
                 "preference for peptide length. The matrices were learned from HLA immunopeptidomics data by motif deconvolution: a computer procedure that "
                 "separates the mixture of motifs present in a sample that carries several alleles. A peptide's score combines the matrix values at its positions; "
                 "the %Rank compares it with random peptides."),
         use="Votes 'binder' at the binding cut-off. Its binding score is also used inside PRIME.",
         limits="Only alleles with a defined matrix are supported (the engine lists the unsupported ones, e.g. HLA-C*07:06).",
         licence="Free for academic use; for-profit use needs a licence (see the tool's repository).",
         refs=["Gfeller D, Schmidt J, Croce G, et al. Improved predictions of antigen presentation and TCR recognition with MixMHCpred2.2 and PRIME2.0 reveal potential natural neoantigens. Nat Biotechnol 2023;41:1790-1798.",
               "Gfeller D, Guillaume P, Michaux J, et al. The length distribution and multiple specificity of naturally presented HLA-I ligands. J Immunol 2018;201:3705-3716.",
               "Bassani-Sternberg M, Chong C, Guillaume P, et al. Deciphering HLA-I motifs across HLA peptidomes improves neo-antigen predictions and identifies allostery regulating HLA specificity. PLoS Comput Biol 2017;13:e1005725."]),
    dict(key="mhcflurry", name="MHCflurry", full="MHCflurry (neural network predictor of class I presentation)", cls="I", role="Affinity, processing and presentation",
         version="2.0", by="O'Donnell, Rubinsteyn and Laserson",
         inputs="Peptides of 8-15 residues, alleles, and (for the processing model) the residues that flank the peptide in its protein.",
         outputs=[("MHCF_Affinity_nM", "Predicted IC50 in nM"), ("MHCF_Affinity_%Rank", "Percentile of the affinity"),
                  ("MHCF_Processing", "Antigen-processing score from the flanking residues"),
                  ("MHCF_Presentation", "Presentation score, 0-1"), ("MHCF_Presentation_%Rank", "Percentile of the presentation score; lower = better")],
         method=("Three parts: (1) an ensemble of neural networks for binding affinity, trained on affinity and mass-spectrometry data; (2) a second network that "
                 "learns from the residues around a peptide whether the cell's machinery tends to produce it (antigen processing); (3) a logistic regression that "
                 "combines the affinity and processing scores into the presentation score."),
         use="Votes 'binder' when the presentation %Rank is at or below the binding cut-off.",
         limits="Class I only. The processing part needs the protein context, so loose peptides get a presentation score without it.",
         licence="Apache-2.0 (open source).",
         refs=["O'Donnell TJ, Rubinsteyn A, Laserson U. MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. Cell Syst 2020;11:42-48."]),
    dict(key="prime", name="PRIME", full="PRedictor of IMmunogenic Epitopes", cls="I", role="Immunogenicity (CD8 T-cell recognition)",
         version="2.1", by="Gfeller laboratory, University of Lausanne",
         inputs="Peptides of 8-14 residues and HLA class I alleles (needs MixMHCpred installed).",
         outputs=[("PRIME_Score", "Immunogenicity score"), ("PRIME_%Rank", "Percentile among random peptides; lower = more likely recognised by CD8 T cells"),
                  ("PRIME_%RankBinding", "MixMHCpred's binding rank used inside PRIME")],
         method=("PRIME multiplies two things: how well MixMHCpred says the peptide binds the allele, and how favourable its amino acids are at the positions that "
                 "face the T-cell receptor. The second part was learned from peptides known to be recognised or ignored by CD8 T cells. The final score is turned into "
                 "a percentile in the same way as MixMHCpred's."),
         use="Votes at the immunogenicity percentile cut-off (default 5 %).",
         limits="Known problem: PRIME 2.1 never finishes for HLA-C*07:06; the engine skips that allele for PRIME and says so.",
         licence="Free for academic use; for-profit use needs a licence.",
         refs=["Schmidt J, Smith AR, Magnin M, et al. Prediction of neo-epitope immunogenicity reveals TCR recognition determinants and provides insight into immunoediting. Cell Rep Med 2021;2:100194.",
               "Gfeller D, et al. Nat Biotechnol 2023;41:1790-1798 (PRIME 2.0, as above)."]),
    dict(key="iedb_imm", name="IEDB class I immunogenicity", full="IEDB class I immunogenicity predictor (Calis et al.)", cls="I", role="Immunogenicity baseline",
         version="2013 model", by="Immune Epitope Database (IEDB); Calis et al.",
         inputs="Peptides of 8-14 residues; a choice of which positions to mask.",
         outputs=[("IEDB_Immunogenicity", "Immunogenicity score; higher = more immunogenic")],
         method=("Each of the 20 amino acids has a fixed immunogenicity value (a published scale; for example Trp and Ile high, Lys and Met low). Each of the nine "
                 "peptide positions has a weight (0, 0, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0 from position 1 to 9; longer peptides repeat the weight 0.30 for the extra middle "
                 "positions). The score is the sum over positions of weight x value, leaving out the masked positions: the MHC anchor residues, whose effect binding "
                 "already captures (by default positions 1, 2 and the C-terminal residue; the engine's other mask choices are per-allele and custom). Higher = more immunogenic."),
         use="A relative vote: yes when the pair is in the top share (default 25 %) of the binders of that allele.",
         limits="A simple, allele-independent baseline; it is built in and always available.",
         licence="Published method; implemented in the engine.",
         refs=["Calis JJA, Maybeno M, Greenbaum JA, et al. Properties of MHC class I presented peptides that enhance immunogenicity. PLoS Comput Biol 2013;9:e1003266."]),
    dict(key="netctlpan", name="NetCTLpan", full="Network-based cytotoxic T lymphocyte epitope predictor, pan-specific", cls="I", role="Processing and presentation",
         version="1.1 (IEDB bundle)", by="DTU Health Tech, distributed in IEDB's T-cell class I tool bundle",
         inputs="A protein sequence (the windows are cut by the tool), peptide lengths 8-11, and alleles.",
         outputs=[("CTLpan_MHC", "MHC-binding component"), ("CTLpan_TAP", "TAP-transport component"), ("CTLpan_Cleavage", "Proteasomal C-terminal cleavage component"),
                  ("CTLpan_Combined", "Combined score"), ("CTLpan_%Rank", "Percentile of the combined score; lower = better")],
         method=("Combines three steps of the pathway into one score: how well the peptide binds the MHC (NetMHCpan), how efficiently the TAP transporter would carry "
                 "it into the endoplasmic reticulum (a weight-matrix model), and whether the proteasome cuts at its C-terminus (NetChop C-terminal network). "
                 "The combined score is converted to a percentile."),
         use="A processing vote (off by default; turn on 'require processing' to demand it).",
         limits="Only 8-11-mers get a score; other lengths get none, and the run says so.",
         licence="Free for academic use via DTU/IEDB.",
         refs=["Stranzl T, Larsen MV, Lundegaard C, Nielsen M. NetCTLpan: pan-specific MHC class I pathway epitope predictions. Immunogenetics 2010;62:357-368."]),
    dict(key="netchop", name="NetChop", full="Network predictor of proteasomal cleavage", cls="I", role="Processing (cleavage)",
         version="3.1 (IEDB bundle)", by="DTU Health Tech, distributed in IEDB's T-cell class I tool bundle",
         inputs="The protein sequence the peptide comes from (it needs the context).",
         outputs=[("NetChop_Cterm", "Probability, 0-1, that the proteasome cleaves after the peptide's last residue; higher = better")],
         method="A neural network that reads the residues around each position of a protein and predicts cleavage after it, trained on proteasome digestion and MHC-ligand data. The engine takes the value at the peptide's C-terminal residue.",
         use="A processing vote when 'netchop_min' is set.",
         limits="Cannot score loose peptides (no protein context); the C-terminus is the part that matters for class I.",
         licence="Free for academic use via DTU/IEDB.",
         refs=["Nielsen M, Lundegaard C, Lund O, Kesmir C. The role of the proteasome in generating cytotoxic T-cell epitopes: insights obtained from improved predictions of proteasomal cleavage. Immunogenetics 2005;57:33-41."]),
    dict(key="bigmhc", name="BigMHC", full="BigMHC (deep-learning MHC class I predictor)", cls="I", role="Presentation (EL) and immunogenicity (IM)",
         version="as installed (ensemble of 7 networks)", by="Karchin laboratory, Johns Hopkins University",
         inputs="Peptide-allele pairs (human HLA class I only).",
         outputs=[("BigMHC_EL", "Presentation likelihood, 0-1; higher = better"), ("BigMHC_IM", "Immunogenicity likelihood, 0-1; higher = better")],
         method=("An ensemble of seven deep neural networks that read the peptide and the MHC pseudosequence. BigMHC EL was trained on peptides eluted from cells; "
                 "BigMHC IM starts from the EL networks and is trained further, by transfer learning, on peptides tested for T-cell responses."),
         use="EL votes 'binder' at score >= 0.2 (default). IM is a relative vote among the binders. IM runs only on pairs already called binders, because it is slow.",
         limits="Slow on CPU; human alleles only.",
         licence="See its repository for the licence.",
         refs=["Albert BA, Yang Y, Shao XM, et al. Deep neural networks predict class I major histocompatibility complex epitope presentation and transfer learn neoepitope immunogenicity. Nat Mach Intell 2023;5:861-872."]),
    dict(key="netmhciipan", name="NetMHCIIpan", full="Network-based MHC class II prediction, pan-specific", cls="II", role="Binding and presentation",
         version="4.3 (also 4.2 and IEDB default)", by="DTU Health Tech; run through IEDB's T-cell class II standalone",
         inputs="Peptides of 9-30 residues (the engine cuts 11-25) and class II alleles: DR (DRB1/3/4/5), DQ and DP as alpha/beta pairs, mouse I-A/I-E.",
         outputs=[("NII_Core", "The predicted 9-residue binding core inside the peptide"), ("NII_Score_EL", "Eluted-ligand score, 0-1"),
                  ("NII_%Rank_EL", "Percentile of the EL score; lower = better"), ("NII_IC50_nM", "Predicted affinity, nM"),
                  ("NII_%Rank_BA", "Percentile of the affinity score")],
         method=("Like NetMHCpan, but a class II molecule holds a 9-residue core with the peptide ends hanging out, so the network first finds the best core "
                 "position and then scores it. The input is the peptide, the pseudosequence(s) of the alpha and beta chains, and (in 4.3) the flanking residues. "
                 "Trained on binding affinities and on mass-spectrometry eluted ligands across DR, DQ and DP. %Rank compares the score with random peptides."),
         use="Two binding votes (EL and BA) at the class II binding cut-off (default 5 %). The core is used to cluster overlapping peptides.",
         limits="Predictions for DQ and DP are less accurate than for DR. Version 4.3 and IEDB's default give different numbers.",
         licence="Free for academic and non-profit use (IEDB / DTU).",
         refs=["Nilsson JB, Kaabinejadian S, Yari H, et al. Accurate prediction of HLA class II antigen presentation across all loci using tailored data acquisition and refined machine learning. Sci Adv 2023;9:eadj6367.",
               "Reynisson B, et al. Nucleic Acids Res 2020;48:W449-W454 (NetMHCIIpan-4.0, as above)."]),
    dict(key="mixmhc2pred", name="MixMHC2pred", full="Mixture MHC class II predictor", cls="II", role="Presentation",
         version="2.1", by="Gfeller laboratory, University of Lausanne",
         inputs="Peptides of any length (the engine sends 11-25) and class II alleles (DR, DQ, DP pairs); uses the residues flanking a peptide when it comes from a protein.",
         outputs=[("MIX2_%Rank", "Percentile among random peptides; lower = better")],
         method=("Position weight matrices for the 9-residue core, for the residues just outside it, and for peptide length, learned from HLA class II "
                 "immunopeptidomes by motif deconvolution (including an alternative, reversed binding mode found for some alleles). The score is converted to a percentile."),
         use="One binding vote at the class II cut-off.",
         limits="Alleles without a matrix are reported as unsupported.",
         licence="Free for academic use; for-profit use needs a licence.",
         refs=["Racle J, Guillaume P, Schmidt J, et al. Machine learning predictions of MHC-II specificities reveal alternative binding mode of class II epitopes. Immunity 2023;56:1359-1375.",
               "Racle J, Michaux J, Rockinger GA, et al. Robust prediction of HLA class II epitopes by deep motif deconvolution of immunopeptidomes. Nat Biotechnol 2019;37:1283-1286."]),
    dict(key="cd4episcore", name="CD4episcore", full="CD4 T-cell epitope immunogenicity score", cls="II", role="Immunogenicity of the peptide (no HLA needed)",
         version="IEDB T-cell class II standalone (v0.2.3-beta)", by="IEDB / La Jolla Institute (Dhanda et al.)",
         inputs="A peptide (the tool takes a protein, so the engine sends each peptide as a one-window sequence).",
         outputs=[("CD4_Episcore", "Percentile; LOWER = more immunogenic")],
         method=("Combines predicted binding of the peptide to a reference panel of seven common HLA-DR molecules with a score derived from experimentally tested CD4 "
                 "epitopes, as described in the paper. The panel makes the score independent of the donor's HLA type; the result is a percentile against random peptides."),
         use="A relative immunogenicity vote (the peptide is ranked among the binder peptides). With the chain on, it only runs for peptides that bind some allele.",
         limits="Peptide-level only, so it cannot tell which allele presents the peptide. On new Intel CPUs its helper programs may crash (see the README).",
         licence="Free for academic and non-profit use (IEDB).",
         refs=["Dhanda SK, Karosiene E, Edwards L, et al. Predicting HLA CD4 immunogenicity in human populations. Front Immunol 2018;9:1369."]),
    dict(key="tlimmuno2", name="TLimmuno2", full="Transfer-learning predictor of MHC class II immunogenicity, version 2", cls="II", role="Immunogenicity of a peptide-allele pair",
         version="(models from the repository)", by="Liu laboratory (XSLiuLab), ShanghaiTech University",
         inputs="Peptide-allele pairs, peptides up to 21 residues, DR/DQ/DP alleles named as in the NetMHCIIpan pseudosequence table.",
         outputs=[("TLimmuno2", "Immunogenicity score, 0-1; higher = more immunogenic")],
         method=("Transfer learning: a network first trained to predict binding affinity (BA) is reused. Its penultimate layer, together with the BLOSUM62-encoded "
                 "peptide (21 positions) and MHC pseudosequence (34 residues), feeds a second network trained on peptides with measured CD4 T-cell responses. "
                 "The engine uses the raw output; the authors' own percentile against 90,001 background peptides is not computed (slow; the engine ranks "
                 "within the candidates instead)."),
         use="A relative vote: yes when the pair is in the top share of the binders of that allele. Runs only on pairs called binders.",
         limits="Needs TensorFlow 2.15 (Keras 2), so its own Python environment. Peptides longer than 21 residues are not scored.",
         licence="See its repository for the licence.",
         refs=["Wang G, Wu T, Ning W, et al. TLimmuno2: predicting MHC class II antigen immunogenicity through transfer learning. Brief Bioinform 2023;24:bbad116. doi:10.1093/bib/bbad116"]),
    dict(key="immuscope", name="ImmuScope", full="ImmuScope (self-iterative multiple-instance learning)", cls="II", role="Immunogenicity of a peptide-allele pair",
         version="IM (10-network ensemble)", by="Shen, Song and colleagues (code and weights on GitHub and Zenodo)",
         inputs="Peptide-allele pairs, peptides up to 21 residues, DR/DQ/DP alleles from its pseudosequence table (about 5,600 names).",
         outputs=[("ImmuScope_IM", "Instance probability averaged over 10 networks, 0-1; higher = more immunogenic")],
         method=("A weakly supervised deep-learning model. Mass-spectrometry data on class II mostly comes from samples with several alleles, so which allele presented "
                 "which peptide is unknown. Multiple-instance learning treats each sample as a 'bag' of possible pairs and iteratively decides which pair is the real one, "
                 "with a triplet loss that pulls true pairs together. The presentation model is then fine-tuned to predict CD4 epitopes and immunogenicity. "
                 "The engine averages the 10 immunogenicity networks."),
         use="A relative vote among the binders of the allele. Runs only on pairs called binders (about 4.5 ms per pair).",
         limits="Peptides longer than 21 residues are not scored. Its own EL and CD4-epitope models exist but are not used.",
         licence="Code GPL-3.0; weights CC-BY-4.0 (Zenodo).",
         refs=["Shen L-C, Zhang Y, Wang Z, Littler DR, Liu Y, Tang J, Rossjohn J, Yu D-J, Song J. Self-iterative multiple-instance learning enables the prediction of CD4+ T cell immunogenic epitopes. Nat Mach Intell 2025. doi:10.1038/s42256-025-01073-z"]),
]

# ------------------------------------------------------------------------------------------------ engine methods
METHODS = [
    ("Votes, agreement and the binder call",
     ["Each tool that scored a peptide-allele pair casts a yes/no vote.",
      "Percentile tools (NetMHCpan, MixMHCpred, MHCflurry, NetMHCIIpan, MixMHC2pred, PRIME, NetCTLpan) vote yes when %Rank <= the cut-off: 2 for class I and 5 for class II "
      "binding by default, 5 for the immunogenicity percentile of PRIME.",
      "BigMHC-EL votes yes when its score >= 0.2.",
      "Tools without an absolute scale (BigMHC-IM, IEDB/Calis, CD4episcore, TLimmuno2, ImmuScope) vote relatively: yes when the pair is in the top q % (default 25) of the "
      "binders of that allele in the run. The ranking needs at least 8 binder pairs; with fewer, the tool gives no vote and the run says so.",
      "Agreement = votes / tools that scored the pair (a tool that did not score it is not counted against it).",
      "A pair is a binder when the binding tools agree as the rule says: any one, a majority (default) or all. It is a strong binder when it is a binder and its median binding "
      "%Rank <= 0.5 (class I) or 2 (class II).",
      "A pair is a candidate pair when it is a binder, the immunogenicity tools agree as the rule says, processing is satisfied if required, and your per-tool filters pass."]),
    ("Priority score",
     ["For every peptide the engine takes four numbers: population coverage of its binder alleles, population-weighted agreement, population-weighted binding rank and best "
      "immunogenicity position among its binder alleles.",
      "Each is converted to a percentile among the peptides of this run (1 = best), the four are averaged with the weights coverage 0.35, agreement 0.25, binding 0.20, "
      "immunogenicity 0.20 (changeable), and the result is multiplied by 100.",
      "So Priority_Score is relative to this run: 100 means best of the run, not a probability of anything. Changing the peptides in the run changes the scores."]),
    ("Population coverage",
     ["Each allele has a frequency in the chosen population (default India, six regions). Coverage of a set of alleles = the share of people who carry at least one of them.",
      "Within a locus the frequencies add. A person has two haplotypes, so the chance of carrying none of the alleles of a locus is (1 - sum of frequencies)^2 (Hardy-Weinberg).",
      "Loci are treated as independent: P(covered) = 1 - the product over loci of P(none at that locus). This ignores linkage between loci and is an approximation.",
      "It is computed per region and averaged with the region weights. Alleles without frequency data count in agreement but not in coverage.",
      "Frequency source (Allele_Data): exact = a typed allele; dqb1_component = a DQA1/DQB1 pair represented by its DQB1 chain; ld_proxy = DRB3/4/5 estimated from the DRB1 alleles "
      "they travel with; none = no data (HLA-DP, mouse and untyped alleles)."]),
    ("Clusters, representatives and the minimum panel",
     ["Overlapping candidate windows are one epitope. Class I: windows that overlap by at least 70 % of the shorter one are merged. Class II: the predicted 9-residue binding cores "
      "must overlap by at least 70 %, so the many 15-25-mers built around one core become a single cluster.",
      "The best-ranked peptide of a cluster is its representative.",
      "The minimum panel is chosen greedily: start empty, add the representative that increases population coverage the most, repeat until the goal (95 % by default) or the cap "
      "(25 peptides) is reached. 'gain' is the coverage a peptide adds; 'cumulative coverage' the total so far."]),
    ("Epitope regions (CTL+HTL)",
     ["Candidate peptides of both classes that overlap on the protein (by 70 % of the shorter one) are grouped. The region is CTL+HTL when it contains a class I and a class II candidate, "
      "otherwise CTL only or HTL only.",
      "Region score = the mean of the Priority_Score of the best class I and best class II peptide, plus 10 for a CTL+HTL region. Regions are ranked CTL+HTL first, then by score."]),
    ("Host-proteome identity (self-similarity)",
     ["The reference proteome of the host (UniProt, one sequence per gene; human about 20,600 proteins and 11.4 million residues; mouse optional) is downloaded once and indexed.",
      "The unit of comparison is the 9-mer, the part of the peptide held in the MHC groove (the whole peptide for 8-mers). A longer peptide is judged by its best 9-mer window.",
      "identical: one of the peptide's 9-mers occurs exactly in a host protein. near: a 9-mer differs from a host 9-mer in one residue (identity >= 8/9, adjustable). distinct: neither.",
      "Search: a 9-mer with at most two changes keeps at least one of its three 3-residue blocks intact, so host windows are looked up by block and compared in full. "
      "This finds every match with up to two mismatches (checked against brute force in the tests); worse matches are reported as distinct and the identity is left blank (below 7/9).",
      "Rules: do not check, flag only (default), remove identical, remove identical and near. Only candidate peptides are checked.",
      "It compares sequences. It does not know whether the host protein is made in a tissue where tolerance matters."]),
    ("Peptide properties for synthesis",
     ["Mass = sum of residue masses + water. GRAVY = mean Kyte-Doolittle hydropathy (above 1 = poorly soluble). Net charge at pH 7.4 from the side-chain and terminal pKa values; pI = the pH of zero net charge.",
      "Notes warn about Cys (oxidation, disulfide formation), Met and Trp (oxidation), N-terminal Gln/Glu (pyroglutamate), Asp-Pro and Asp-Gly (acid-labile, aspartimide), hydrophobic runs of 5 or more, and GRAVY above 1."]),
    ("The funnel (running slow tools only where needed)",
     ["Quick tools score every pair first. A slow tool then scores only the pairs that can still become binders or candidates under your rules.",
      "The gate is exact for the agreement rule: a pair is kept when the votes it already has plus every vote still pending could reach the number the rule needs. A pair that fails "
      "this can never pass, so skipping it changes nothing in the result.",
      "Manual mode lets you define the stages yourself; 'all models' disables the gate."]),
]

# ------------------------------------------------------------------------------------------------ columns
# (column or prefix, meaning). A name ending in * is a prefix.
COLUMNS = [
    ("Protein", "Name of the protein (FASTA header) the window was cut from; 'peptide' for loose peptides."),
    ("Start", "1-based position of the peptide's first residue in the protein."),
    ("End", "Position of the peptide's last residue."),
    ("Peptide", "The peptide sequence (one-letter amino-acid code)."),
    ("Length", "Number of residues."),
    ("Allele", "HLA allele. Class II pairs are written HLA-DQA1*05:01/DQB1*02:01; mouse alleles as H2-IAb."),
    ("Allele_Freq", "Frequency of the allele in the chosen population (national)."),
    ("Allele_Data", "Where the allele's frequency came from: exact, dqb1_component, ld_proxy or none (no data)."),
    ("Alleles_Tested", "How many alleles were scored for the peptide."),
    ("Alleles_tested", "How many alleles were scored for the peptide."),
    ("Alleles_passed", "In how many alleles every selected tool passed."),
    ("Alleles_passed_list", "Those alleles."),
    ("NMP_Score_EL", "NetMHCpan eluted-ligand score, 0-1 (higher = more likely presented)."),
    ("NMP_%Rank_EL", "NetMHCpan eluted-ligand percentile among random peptides (lower = better)."),
    ("NMP_Score_BA", "NetMHCpan binding-affinity score, 0-1."),
    ("NMP_%Rank_BA", "NetMHCpan binding-affinity percentile."),
    ("NMP_Aff_nM", "NetMHCpan predicted IC50 in nanomolar (lower = tighter binding)."),
    ("NMP_BindLevel", "NetMHCpan label: SB strong binder, WB weak binder."),
    ("MixMHC_Score", "MixMHCpred raw score."), ("MixMHC_%Rank", "MixMHCpred percentile (lower = better)."),
    ("MHCF_Affinity_nM", "MHCflurry predicted IC50 in nM."), ("MHCF_Affinity_%Rank", "MHCflurry affinity percentile."),
    ("MHCF_Processing", "MHCflurry antigen-processing score (from the flanking residues)."),
    ("MHCF_Presentation", "MHCflurry presentation score, 0-1."), ("MHCF_Presentation_%Rank", "MHCflurry presentation percentile (lower = better)."),
    ("PRIME_Score", "PRIME immunogenicity score."), ("PRIME_%Rank", "PRIME percentile (lower = more likely recognised by CD8 T cells)."),
    ("PRIME_%RankBinding", "MixMHCpred binding rank used inside PRIME."),
    ("IEDB_Immunogenicity", "Calis class I immunogenicity score (higher = more immunogenic)."),
    ("CTLpan_MHC", "NetCTLpan MHC-binding component."), ("CTLpan_TAP", "NetCTLpan TAP-transport component."),
    ("CTLpan_Cleavage", "NetCTLpan proteasomal-cleavage component."), ("CTLpan_Combined", "NetCTLpan combined score."),
    ("CTLpan_%Rank", "NetCTLpan percentile of the combined score (lower = better)."),
    ("NetChop_Cterm", "NetChop probability (0-1) that the proteasome cuts after the peptide's last residue."),
    ("BigMHC_EL", "BigMHC presentation likelihood, 0-1."), ("BigMHC_IM", "BigMHC immunogenicity likelihood, 0-1."),
    ("NII_Core", "NetMHCIIpan predicted 9-residue binding core."), ("NII_Score_EL", "NetMHCIIpan eluted-ligand score, 0-1."),
    ("NII_%Rank_EL", "NetMHCIIpan eluted-ligand percentile (lower = better)."), ("NII_IC50_nM", "NetMHCIIpan predicted affinity in nM."),
    ("NII_%Rank_BA", "NetMHCIIpan binding-affinity percentile."), ("MIX2_%Rank", "MixMHC2pred percentile (lower = better)."),
    ("CD4_Episcore", "CD4episcore percentile of the peptide (lower = more immunogenic); needs no allele."),
    ("TLimmuno2", "TLimmuno2 immunogenicity score, 0-1 (higher = more immunogenic)."),
    ("ImmuScope_IM", "ImmuScope immunogenicity score, 0-1, mean of 10 networks (higher = more immunogenic)."),
    ("Best_NMP_%Rank_EL", "Lowest NetMHCpan EL percentile over the alleles."), ("Best_MixMHC_%Rank", "Lowest MixMHCpred percentile over the alleles."),
    ("Best_NII_%Rank_EL", "Lowest NetMHCIIpan EL percentile over the alleles."), ("Best_NII_%Rank_BA", "Lowest NetMHCIIpan BA percentile over the alleles."),
    ("Best_MIX2_%Rank", "Lowest MixMHC2pred percentile over the alleles."),
    ("Max_ImmuScope_IM", "Highest ImmuScope score over the alleles."), ("Max_TLimmuno2", "Highest TLimmuno2 score over the alleles."),
    ("Max_BigMHC_IM", "Highest BigMHC-IM over the alleles."), ("Max_BigMHC_EL", "Highest BigMHC-EL over the alleles."),
    ("Pass_All", "True when every selected tool scored the pair and met your per-tool filter (and, in a manual funnel, survived every stage)."),
    ("Pass_*", "True when that tool scored the pair and met the filter you set for it (no filter = true if it scored the pair)."),
    ("Vote_*", "That tool's vote on the pair: 1 yes, 0 no, empty = it did not score the pair or gave no vote."),
    ("TopPct_*", "Position of the pair among the binders of its allele for a relative tool, in percent (lower = better)."),
    ("Funnel_Alive", "Manual funnel: the pair passed every stage."), ("Funnel_Reached", "Manual funnel: the stages that scored it."),
    ("Stopped_At", "Manual funnel: the stage that stopped it (0 = it came through)."),
    ("S*", "Manual funnel: S<k>_In = scored in stage k; S<k>_Pass = let through stage k."),
    ("Priority_Rank", "Rank of the peptide, 1 = best, by the 'rank by' setting."),
    ("Priority_Score", "0-100 relative standing of the peptide in this run: weighted mean of four percentiles (coverage, agreement, binding, immunogenicity)."),
    ("Pair_Call", "'strong binder' or 'binder' for the pair, else empty."),
    ("Pair_Candidate", "The pair is a binder, immunogenic enough, processed if required, and passes your filters."),
    ("Binder", "The binding tools agree as the rule says (default: a majority)."),
    ("Strong", "A binder whose median binding %Rank is at or below the strong cut-off (0.5 class I, 2 class II)."),
    ("Bind_Consensus_%Rank", "Median %Rank over the binding tools that scored the pair."),
    ("Bind_Consensus_Best", "That median for the peptide's best allele."),
    ("Bind_Consensus_PopW", "Population-weighted geometric mean of the binding %Rank over the tested alleles."),
    ("Agree_Votes", "Number of yes votes over all tool groups."), ("Agree_Tools", "Number of tools that scored the pair."),
    ("Agree_Fraction", "Agree_Votes / Agree_Tools."), ("Agree_List", "The tools that voted yes."),
    ("Agree_Fraction_Best", "Agreement for the peptide's best allele."),
    ("Agree_Fraction_PopW", "Agreement averaged over the tested alleles, weighted by allele frequency."),
    ("Bind_Votes", "Yes votes of the binding tools."), ("Bind_Tools", "Binding tools that scored the pair."),
    ("Imm_Votes", "Yes votes of the immunogenicity tools."), ("Imm_Tools", "Immunogenicity tools that scored the pair."),
    ("Proc_Votes", "Yes votes of the processing tools."), ("Proc_Tools", "Processing tools that scored the pair."),
    ("Imm_Votes_Best", "Immunogenicity yes votes for the best allele."), ("Imm_Tools_Best", "Immunogenicity tools for the best allele."),
    ("Proc_Votes_Best", "Processing yes votes for the best allele."),
    ("Imm_TopPct_Consensus", "Mean of the immunogenicity tools' TopPct for the pair (lower = better)."),
    ("Imm_TopPct_Best", "Best (lowest) Imm_TopPct_Consensus over the peptide's binder alleles."),
    ("Imm_OK", "The immunogenicity tools agree as the rule says."), ("Proc_OK", "Processing satisfied (or not required)."),
    ("Filters_OK", "Your per-tool filters passed."), ("Plausible", "The pair is a binder, so immunogenicity tools rank it among binders."),
    ("Candidate", "The peptide passes the peptide-level rule: enough alleles, enough coverage, and not removed by the host rule."),
    ("N_Alleles_Binder", "In how many alleles the peptide is a binder."), ("N_Alleles_Strong", "...a strong binder."),
    ("N_Alleles_Candidate", "...a candidate pair."), ("Alleles_Binder", "Those binder alleles."), ("Alleles_Candidate", "Those candidate alleles."),
    ("Best_Allele", "The binder allele with the lowest binding %Rank."),
    ("Best_Core", "Class II: predicted 9-residue core of the best allele."), ("Core_Start", "Protein position of the core's first residue."),
    ("Core_End", "Protein position of the core's last residue."),
    ("Coverage_Binder", "Population coverage of the alleles for which the peptide is a binder (a column per region follows)."),
    ("Coverage_Candidate", "Population coverage of the alleles for which it is a candidate."),
    ("Mass_Da", "Monoisotopic-style mass of the peptide in daltons."), ("GRAVY", "Grand average of hydropathy (Kyte-Doolittle); above 1 = poorly soluble."),
    ("Net_Charge_pH7.4", "Net charge at pH 7.4."), ("pI", "Isoelectric point."), ("Aromatic_Frac", "Fraction of F, W, Y."),
    ("Hydrophobic_Frac", "Fraction of hydrophobic residues."), ("N_Cys", "Number of cysteines (oxidation, disulfides)."),
    ("N_Met", "Number of methionines (oxidation)."), ("N_Trp", "Number of tryptophans (oxidation)."),
    ("Synthesis_Notes", "Warnings for making the peptide (Cys, Met/Trp, N-terminal Gln/Glu, Asp-Pro/Asp-Gly, hydrophobic runs, poor solubility)."),
    ("Host_Self", "Comparison with the host proteome: identical (an exact 9-mer match), near (one mismatch), distinct."),
    ("Host_Exact_Kmers", "How many of the peptide's 9-mers occur exactly in the host proteome."),
    ("Host_Best_Identity", "Best share of identical residues over a 9-mer window against any host 9-mer (empty = below 7/9)."),
    ("Host_Best_Protein", "UniProt entry holding the closest host match."), ("Host_Best_Window", "The host 9-mer that is closest."),
    ("Host_Best_Start", "Its position in that host protein."), ("Host_Excluded", "Removed from the candidates by your host rule."),
    ("Cluster_ID", "Number of the cluster (overlapping peptides = one epitope)."), ("Cluster_Size", "Peptides in the cluster."),
    ("Is_Representative", "This peptide stands for its cluster."), ("Region", "Protein span of the cluster, Protein:first-last."),
    ("Cluster_Representative", "The peptide that stands for this one."), ("Cluster_Overlap", "Overlap with the representative as a share of the shorter one."),
    ("Cluster_Members", "On a representative: the other peptides of its cluster, SEQ@start."),
    ("In_Panel", "Part of the minimum panel."), ("Panel_Order", "Order in which it was added to the panel."),
    ("Panel_Cum_Coverage", "Population coverage reached when it was added."),
    ("Region_ID", "Label of the epitope region (R001...)."), ("Region_Type", "CTL+HTL, CTL only or HTL only."), ("Type", "CTL+HTL, CTL only or HTL only."),
    ("Region_Rank", "Rank of the region (CTL+HTL first, then score)."), ("Region_Score", "Mean Priority_Score of the best CTL and HTL peptide, +10 for CTL+HTL."),
    ("Sequence", "Protein sequence covered by the region."),
    ("CTL_Peptides", "Class I candidate windows inside the region."), ("HTL_Peptides", "Class II candidate windows inside the region."),
    ("CTL_Best_Peptide", "Best class I peptide of the region."), ("HTL_Best_Peptide", "Best class II peptide of the region."),
    ("CTL_Best_Allele", "Its best allele."), ("HTL_Best_Allele", "Its best allele."),
    ("CTL_Coverage", "Its population coverage."), ("HTL_Coverage", "Its population coverage."),
    ("CTL_Priority_Score", "Its priority score."), ("HTL_Priority_Score", "Its priority score."),
    ("Synthesis_Notes", "Warnings for making the peptide or region."),
    ("Best_CTLpan_%Rank", "Lowest NetCTLpan percentile over the alleles."), ("Best_MHCF_Presentation_%Rank", "Lowest MHCflurry presentation percentile over the alleles."),
    ("Best_NII_IC50_nM", "Lowest NetMHCIIpan predicted IC50 (nM) over the alleles."), ("Best_NMP_Aff_nM", "Lowest NetMHCpan predicted IC50 (nM) over the alleles."),
    ("Best_PRIME_%Rank", "Lowest PRIME percentile over the alleles."), ("Max_IEDB_Immunogenicity", "Highest Calis immunogenicity score over the alleles."),
    ("Max_NetChop_Cterm", "Highest NetChop C-terminal cleavage probability over the alleles."),
    ("Core", "Predicted 9-residue binding core inside the peptide."),
    ("protein", "Protein the panel peptide comes from."), ("start", "Its first position in the protein."), ("length", "Its length."),
    ("region", "Protein span of its cluster."),
    ("order", "Order of addition to the minimum panel."), ("peptide", "Panel peptide."), ("gain", "Population coverage this peptide adds."),
    ("cumulative_coverage", "Coverage reached so far."), ("new_alleles", "Alleles this peptide newly covers."),
]

# ------------------------------------------------------------------------------------------------ result tiles
TILES = {
    "peptides": "Distinct peptide sequences scored (windows of the protein, or the peptides you gave).",
    "regions": "Epitope regions: stretches of the protein where candidate peptides overlap.",
    "dual": "CTL+HTL regions: a class I and a class II candidate overlap, so one stretch of sequence serves both T-cell arms.",
    "panel": "Population coverage reached by the minimum panel, with the number of peptides it needs.",
    "host": "Candidates that look like the host's own proteins: identical (an exact 9-mer) plus near (one mismatch), out of the candidates checked.",
    "time": "Total run time. 'Chained' means the slow tools scored only the pairs the quick tools let through.",
}

# ------------------------------------------------------------------------------------------------ glossary
TERMS = [
    ("HLA", "Human Leukocyte Antigen", "The human version of the MHC. Genes named HLA-A, -B, -C (class I) and HLA-DR, -DQ, -DP (class II)."),
    ("MHC", "Major Histocompatibility Complex", "Molecules that display peptide fragments on the cell surface for T cells to inspect."),
    ("Allele", "", "One version of a gene. HLA-A*02:01: gene A, allele group 02, protein subtype 01. Thousands exist; people differ in which they carry."),
    ("Class I / class II", "", "Class I shows peptides made inside the cell to CD8 T cells; class II shows peptides taken up from outside to CD4 T cells."),
    ("CTL", "Cytotoxic T Lymphocyte (CD8+)", "Kills infected or abnormal cells. Recognises class I peptides of 8-14 residues."),
    ("HTL", "Helper T Lymphocyte (CD4+)", "Coordinates the immune response and helps B cells and CTLs. Recognises class II peptides of 11-25 residues."),
    ("TCR", "T-Cell Receptor", "The receptor on a T cell that recognises a peptide-MHC complex."),
    ("Epitope", "", "The part of an antigen that the immune system recognises; here the peptide a T cell sees."),
    ("Binding core", "", "The 9 residues of a class II peptide that sit in the MHC groove; the rest overhangs both ends."),
    ("Anchor", "", "A peptide position that fits a pocket of the MHC and mainly decides binding, not T-cell recognition."),
    ("Pseudosequence", "", "The ~34 MHC residues that contact the peptide; prediction networks use it to cover alleles they were never trained on."),
    ("Pan-specific", "", "One model for all alleles instead of one per allele."),
    ("EL", "Eluted Ligand", "A peptide that was actually found on cells by mass spectrometry. EL models learn presentation, including processing."),
    ("BA", "Binding Affinity", "How tightly a peptide binds an MHC in the test tube. BA models learn binding only."),
    ("IC50", "Half-maximal inhibitory concentration", "The peptide concentration that displaces half of a reference peptide. Lower = tighter. Under 50 nM is strong, under 500 nM weak by convention."),
    ("%Rank", "Percentile rank", "Share of random natural peptides that score better for the same allele. 1 % means the peptide is in the best 1 %. It makes tools comparable across alleles."),
    ("SB / WB", "Strong / Weak Binder", "Labels by %Rank thresholds, e.g. class I EL: SB below 0.5 %, WB below 2 %."),
    ("Presentation", "", "The peptide is produced, transported and displayed by an MHC molecule."),
    ("Immunogenicity", "", "The displayed peptide actually triggers a T-cell response. Much harder to predict than binding."),
    ("Proteasome", "", "The cell's protein-cutting machine; makes the C-terminus of class I peptides."),
    ("TAP", "Transporter associated with Antigen Processing", "Carries peptides into the endoplasmic reticulum where class I molecules wait."),
    ("Processing", "", "Everything that makes a peptide available to MHC: cleavage, transport, trimming."),
    ("Consensus", "", "A combined call of several tools."),
    ("Vote / agreement", "", "Each tool says yes or no; agreement is the share of yes among the tools that scored the pair."),
    ("Gate / chain / funnel", "", "Running slow tools only on pairs that can still pass."),
    ("Cluster", "", "Overlapping peptides that count as one epitope."),
    ("Representative", "", "The best peptide of a cluster."),
    ("Minimum panel", "", "The smallest set of representatives reaching the coverage goal."),
    ("Population coverage", "", "The share of the population carrying at least one allele that presents the peptide."),
    ("Hardy-Weinberg", "", "A genetics rule turning allele frequencies into the share of people who carry them."),
    ("Promiscuous epitope", "", "A peptide presented by many different alleles."),
    ("Tolerance", "", "The immune system's learned ignoring of the body's own proteins."),
    ("Autoimmunity", "", "An immune attack on the body's own tissue."),
    ("Reference proteome", "", "A curated, non-redundant set of the proteins of an organism, one sequence per gene (UniProt)."),
    ("k-mer", "", "A stretch of k residues; a 9-mer is nine in a row."),
    ("GRAVY", "Grand Average of Hydropathy", "Mean Kyte-Doolittle hydropathy of a peptide; positive = hydrophobic."),
    ("pI", "Isoelectric point", "The pH at which the peptide has no net charge."),
    ("BLOSUM62", "Blocks Substitution Matrix 62", "A table of how similar amino acids are; used to encode peptides for neural networks."),
    ("Transfer learning", "", "Reusing a network trained on a large related task (binding) as the starting point for a small one (immunogenicity)."),
    ("Multiple-instance learning", "", "Training when labels belong to a group (a sample with several alleles) instead of to each item."),
    ("Ensemble", "", "Several networks whose outputs are averaged."),
    ("Watchdog", "", "The engine kills a tool that shows no activity for the stall time instead of hanging."),
    ("Prediction cache", "", "Stored tool results; a pair already scored is not scored again."),
]

SECTIONS_OVERVIEW = [
    ("What the engine does",
     ["You give a protein (or peptides) and HLA alleles. The engine cuts the protein into peptides, asks many prediction programs how well each peptide is presented by each allele "
      "and how immunogenic it looks, combines the answers, weighs them by how common each allele is in the chosen population, checks each candidate against the host's own proteins, "
      "and returns a ranked list of candidate T-cell epitopes, grouped into regions of the protein and reduced to the smallest peptide set that covers the population.",
      "MHC class I (HLA-A, -B, -C) presents peptides of 8-14 residues to CD8 T cells. MHC class II (HLA-DR, -DQ, -DP) presents peptides of 11-25 residues to CD4 T cells. "
      "A region where both kinds of candidate overlap (CTL+HTL) is usually the first choice for a vaccine.",
      "Everything runs on your computer. Each prediction program runs as a separate process; the web server never loads torch or TensorFlow.",
      "Not covered: antibody (B-cell) epitopes, class II proteasomal processing, and the safety filters (antigenicity, allergenicity, toxicity, conservation across variants), which are separate steps."]),
    ("The chain of a run, in order",
     ["1. Input: the protein or peptides are cut into windows.",
      "2. Quick binding tools score every peptide-allele pair (class I: NetMHCpan, MixMHCpred, MHCflurry, PRIME, IEDB immunogenicity; class II: MixMHC2pred, NetMHCIIpan).",
      "3. The gate keeps only the pairs that can still pass your rule.",
      "4. Slow tools score the survivors (NetCTLpan, NetChop, BigMHC; class II: CD4episcore, TLimmuno2, ImmuScope).",
      "5. Consensus: votes, agreement, binder and candidate calls.",
      "6. Peptide table: one row per peptide over its alleles, population coverage, priority score.",
      "7. Host proteome identity: each candidate is compared with the host's proteins; identical or near ones are flagged or removed.",
      "8. Clusters, representatives, minimum panel, then regions that combine class I and class II.",
      "9. Report, tables and a zip. The predictions are kept, so another strategy can re-rank them in seconds."]),
]


# ------------------------------------------------------------------------------------------------ lookups for the web page
def column_help():
    """{column or prefix*: text} for the '?' on table headers."""
    return {k: v for k, v in COLUMNS}


def lookup(c):
    """The help text for a column name, with the family and region rules the web page applies; None when undocumented."""
    import re
    m = column_help()
    if c in m:
        return m[c]
    for p in ("Vote_", "TopPct_", "Pass_"):
        if c.startswith(p):
            return m[p + "*"]
    if re.match(r"^S\d+_(In|Pass)$", c):
        return m["S*"]
    mm = re.match(r"^(Coverage_(?:Binder|Candidate))_(.+)$", c)
    return m[mm.group(1)] + f" ({mm.group(2)} region only)" if mm and mm.group(1) in m else None


def api():
    return {"columns": column_help(), "tiles": TILES, "version": ""}


# ------------------------------------------------------------------------------------------------ page
CSS = """
:root{--bg:#fff;--fg:#1d2430;--mu:#5b6678;--card:#f6f8fb;--line:#d8dee8;--ac:#2f6fdb}
@media(prefers-color-scheme:dark){:root{--bg:#10151f;--fg:#e5e9f1;--mu:#98a4ba;--card:#182030;--line:#2b384e;--ac:#7aa5ff}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);font:15px/1.6 system-ui,Segoe UI,Roboto,sans-serif}
header{padding:28px 24px 8px;max-width:1100px;margin:0 auto}header h1{margin:0 0 4px;font-size:28px}header p{margin:0;color:var(--mu)}
.wrap{display:flex;gap:28px;max-width:1100px;margin:0 auto;padding:12px 24px 80px}
nav{position:sticky;top:12px;align-self:flex-start;width:230px;flex:none;max-height:calc(100vh - 24px);overflow:auto;font-size:13.5px}
nav a{display:block;color:var(--fg);text-decoration:none;padding:3px 8px;border-radius:6px}nav a:hover{background:var(--card)}
nav .g{margin:12px 8px 2px;color:var(--mu);font-size:11.5px;text-transform:uppercase;letter-spacing:.06em}
main{min-width:0;flex:1}
h2{font-size:21px;margin:40px 0 10px;padding-bottom:6px;border-bottom:1px solid var(--line)}
h3{font-size:17px;margin:26px 0 6px}h4{font-size:12.5px;margin:14px 0 4px;color:var(--mu);text-transform:uppercase;letter-spacing:.05em}
p{margin:6px 0}ul{margin:6px 0 6px 20px;padding:0}li{margin:3px 0}
table{border-collapse:collapse;width:100%;font-size:13.5px;margin:6px 0}th,td{border:1px solid var(--line);padding:5px 9px;text-align:left;vertical-align:top}
th{background:var(--card);font-weight:600}code,.mono{font-family:ui-monospace,Consolas,monospace;font-size:12.5px}
.model{border:1px solid var(--line);border-radius:10px;padding:4px 18px 14px;margin:18px 0;background:var(--bg)}
.model>h3{margin-top:14px}.sub{color:var(--mu);margin:-2px 0 6px}
.badge{display:inline-block;font-size:11.5px;padding:1px 9px;border-radius:10px;background:var(--ac);color:#fff;margin-left:8px;vertical-align:middle}
dl{display:grid;grid-template-columns:150px 1fr;gap:6px 14px;margin:8px 0}dt{color:var(--mu);font-weight:600}dd{margin:0}
.ref{font-size:13.5px}.note{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px 12px;color:var(--mu);font-size:13.5px}
@media(max-width:820px){.wrap{display:block}nav{position:static;width:auto;max-height:none;margin-bottom:16px}dl{grid-template-columns:1fr}}
@media print{nav{display:none}.wrap{display:block}}
"""


def _cls_badge(m):
    return f'<span class="badge">MHC class {m["cls"]}</span>'


def page():
    out = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
           f"<title>Epitope Filter Engine: documentation</title><style>{CSS}</style></head><body>"]
    out.append("<header><h1>Epitope Filter Engine</h1><p>Documentation: what each model is, what goes in and out, how every number is calculated, "
               "what every column and term means, and where the methods come from.</p></header><div class='wrap'><nav>")
    nav = [("Overview", [("overview", "What it does"), ("chain", "The chain of a run")]),
           ("MHC class I models", [(f"m-{m['key']}", m["name"]) for m in MODELS if m["cls"] == "I"]),
           ("MHC class II models", [(f"m-{m['key']}", m["name"]) for m in MODELS if m["cls"] == "II"]),
           ("Engine methods", [(f"x{i}", t) for i, (t, _) in enumerate(METHODS)]),
           ("Reference", [("columns", "Column headers"), ("tiles", "Result tiles"), ("terms", "Glossary and full forms"), ("sources", "Reference list")])]
    for g, items in nav:
        out.append(f"<div class='g'>{E(g)}</div>" + "".join(f"<a href='#{i}'>{E(t)}</a>" for i, t in items))
    out.append("</nav><main>")

    out.append("<h2 id='overview'>What the engine does</h2>" + "".join(f"<p>{E(p)}</p>" for p in SECTIONS_OVERVIEW[0][1]))
    out.append("<h2 id='chain'>The chain of a run, in order</h2><ul>" + "".join(f"<li>{E(p)}</li>" for p in SECTIONS_OVERVIEW[1][1]) + "</ul>")

    for cls, title in (("I", "MHC class I models"), ("II", "MHC class II models")):
        out.append(f"<h2>{title}</h2>")
        for m in (x for x in MODELS if x["cls"] == cls):
            out.append(f"<div class='model' id='m-{m['key']}'><h3>{E(m['name'])} {E(m['version'])}{_cls_badge(m)}</h3>"
                       f"<div class='sub'>{E(m['full'])}. Role in the engine: {E(m['role'])}.</div><dl>"
                       f"<dt>Developed by</dt><dd>{E(m['by'])}</dd><dt>Input</dt><dd>{E(m['inputs'])}</dd></dl>"
                       "<h4>Output columns</h4><table><tr><th>Column</th><th>Meaning</th></tr>"
                       + "".join(f"<tr><td class='mono'>{E(c)}</td><td>{E(t)}</td></tr>" for c, t in m["outputs"]) + "</table>"
                       f"<h4>How it is calculated</h4><p>{E(m['method'])}</p>"
                       f"<h4>How the engine uses it</h4><p>{E(m['use'])}</p>"
                       f"<h4>Limits</h4><p>{E(m['limits'])}</p>"
                       f"<h4>Licence</h4><p>{E(m['licence'])}</p>"
                       "<h4>References</h4><ul class='ref'>" + "".join(f"<li>{E(r)}</li>" for r in m["refs"]) + "</ul></div>")

    out.append("<h2>Engine methods</h2>")
    for i, (t, paras) in enumerate(METHODS):
        out.append(f"<h3 id='x{i}'>{E(t)}</h3><ul>" + "".join(f"<li>{E(p)}</li>" for p in paras) + "</ul>")

    out.append("<h2 id='columns'>Column headers</h2><p class='note'>The same texts appear as a <b>?</b> on every column header of the result tables in the web page. "
               "A name ending in * stands for a family (for example Vote_NII_EL, Vote_MIX2 ...).</p><table><tr><th>Column</th><th>Meaning</th></tr>"
               + "".join(f"<tr><td class='mono'>{E(c)}</td><td>{E(t)}</td></tr>" for c, t in COLUMNS) + "</table>")
    out.append("<h2 id='tiles'>Result tiles</h2><table><tr><th>Tile</th><th>Meaning</th></tr>"
               + "".join(f"<tr><td>{E(k)}</td><td>{E(v)}</td></tr>" for k, v in TILES.items()) + "</table>")
    out.append("<h2 id='terms'>Glossary and full forms</h2><table><tr><th>Term</th><th>Full form</th><th>Meaning</th></tr>"
               + "".join(f"<tr><td><b>{E(a)}</b></td><td>{E(b)}</td><td>{E(c)}</td></tr>" for a, b, c in TERMS) + "</table>")
    out.append("<h2 id='sources'>Reference list</h2><p class='note'>" + E(REF_NOTE) + "</p><ul class='ref'>"
               + "".join(f"<li>{E(r)}</li>" for r in sorted({r for m in MODELS for r in m["refs"]}))
               + "<li>Bui HH, Sidney J, Dinh K, Southwood S, Newman MJ, Sette A. Predicting population coverage of T-cell epitope-based diagnostics and vaccines. BMC Bioinformatics 2006;7:153.</li>"
               "<li>Kyte J, Doolittle RF. A simple method for displaying the hydropathic character of a protein. J Mol Biol 1982;157:105-132.</li>"
               "<li>The UniProt Consortium. UniProt: the Universal Protein Knowledgebase in 2023. Nucleic Acids Res 2023;51:D523-D531.</li>"
               "<li>Vita R, Mahajan S, Overton JA, et al. The Immune Epitope Database (IEDB): 2018 update. Nucleic Acids Res 2019;47:D339-D343.</li></ul>")
    out.append("</main></div></body></html>")
    return "".join(out)
