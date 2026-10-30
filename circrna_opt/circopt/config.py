"""Parameters for a design run. Everything user-tunable lives here (JSON round-trippable)."""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields

SPECIES = ("human", "mouse")

# Common restriction sites (DNA alphabet; converted to RNA + both strands at compile time).
DEFAULT_RE_SITES = {
    "BsaI": "GGTCTC", "BsmBI": "CGTCTC", "BbsI": "GAAGAC", "SapI": "GCTCTTC",
    "EcoRI": "GAATTC", "BamHI": "GGATCC", "HindIII": "AAGCTT", "XhoI": "CTCGAG",
    "NotI": "GCGGCCGC", "NcoI": "CCATGG", "NheI": "GCTAGC", "XbaI": "TCTAGA",
}

DEFAULT_WEIGHTS = {
    # codon fitness
    "codon": 2.0,          # (1 - weighted geometric-mean CAI)
    "cai_floor": 6.0,      # per unit shortfall of the weaker species below cai_floor
    # structure (cheap, MFE based)
    "dsrna": 4.0,          # nt sitting in helices >= dsrna_min_bp (PKR/MDA5 proxy) -- main driver
    "mfe": 4.0,            # MFE per nt, direction set by structure_goal
    "paired": 2.0,         # paired fraction of MFE structure
    "stem": 1.5,           # longest continuous stem beyond stem_target
    "long_range": 1.0,     # long-range pairs
    "start": 3.0,          # start-codon context paired (in the MFE structure)
    "start_access": 10.0,  # AUG accessibility from the fast RNAplfold pass (this drives the GA)
    "start_context": 4.0,  # mean accessibility of the start context
    "early_access": 2.0,   # mean accessibility of the early ORF
    "early": 1.5,          # early ORF paired
    "min_paired": 2.0,     # penalty for being *too* unstructured (hydrolysis caveat)
    # ensemble (partition function, elites only)
    "entropy": 1.0, "diversity": 1.0, "centroid": 1.0, "mfe_freq": 0.5,
    "red": 1.0, "start_prob": 3.0,
    # sequence-level
    "gc": 3.0, "gc_low": 3.0, "motif": 2.0, "inverted_repeat": 0.5, "ggg": 0.2,
}


@dataclass
class Params:
    name: str = "circ_design"
    protein: str = ""                 # amino acids (optional if cds_start/cds_end given)
    rna: str = ""                     # full circular molecule as given (IRES, UTRs, spacers, CDS)
    cds_start: int | None = None      # 0-based; auto-detected from `protein` when None
    cds_end: int | None = None        # half-open
    protected: list = field(default_factory=list)   # [{"name","start","end"}], 0-based half-open

    # species
    primary_species: str = "human"
    secondary_species: str = "mouse"
    primary_weight: float = 0.6       # weight in the weighted geometric mean (secondary = 1 - this)
    cai_floor: float = 0.75           # hard floor on the weaker species' CAI
    # Codon fitness SATURATES here: once the weighted geometric mean reaches cai_target the
    # codon term stops pushing, so the optimiser spends the remaining freedom on structure
    # instead of trading a lot of base-pairing for a few points of CAI. Set it to 1.0 for
    # plain "maximise the geometric mean" behaviour.
    cai_target: float = 0.92
    codon_metric: str = "cai"         # "cai" | "tai" (tai needs trna_copy_file)
    trna_copy_file: str | None = None # CSV: species,anticodon,copies
    min_codon_freq: float = 0.0       # drop codons rarer than this (per-aa freq) in either species

    # sequence composition
    gc_min: float = 0.40
    gc_max: float = 0.55
    gc_window: int = 60
    gc_window_min: float = 0.30
    gc_window_max: float = 0.65
    max_poly_u: int = 4               # runs of U longer than this are forbidden
    max_homopolymer: int = 6          # A/C/G runs longer than this are forbidden
    avoid_slippery: bool = True       # X XXY YYZ ribosomal frameshift heptamers
    avoid_g4: bool = True             # G-quadruplex G3+N1-7 x4
    avoid_cryptic_splice: bool = True
    avoid_re_sites: bool = True
    re_sites: dict = field(default_factory=lambda: dict(DEFAULT_RE_SITES))
    extra_forbidden: list = field(default_factory=list)   # extra IUPAC motifs (DNA or RNA alphabet)
    ir_kmer: int = 9                  # inverted-repeat (complementary repeat) k-mer length

    # folding (ViennaRNA, matches RNAfold web-server defaults)
    temperature: float = 37.0
    dangles: int = 2
    no_lonely_pairs: bool = False
    no_gu: bool = False
    no_gu_closure: bool = False       # "no GU at helix ends"
    gquad: bool = False
    inner_span: int = 200             # max_bp_span of the fast in-loop fold (0 = unlimited)
    stress_test: bool = True          # also fold top variants with noLP + noGUclosure
    linear_compare: bool = True

    # structure objective direction.
    #   "minimize_pairing" : push MFE per nt toward 0 (least structured)
    #   "minimize_mfe"     : push MFE per nt below mfe_per_nt_target (LinearDesign/DERNA direction)
    #   "target"           : hold MFE per nt near mfe_per_nt_target
    # Total pairing is only a soft term in every mode; the hard immunogenicity proxy is
    # `dsrna_min_bp` (long contiguous duplex), which is what PKR/MDA5 actually sense.
    structure_goal: str = "minimize_pairing"
    mfe_per_nt_target: float = -0.45
    dsrna_min_bp: int = 30          # helices at/above this length count as dsRNA

    # structure targets
    stem_target: int = 20
    stem_cap: int = 30
    long_range_nt: int = 150
    red_prob: float = 0.8
    start_upstream: int = 12
    start_downstream: int = 18
    early_orf_nt: int = 45
    # Initiation-region accessibility. The hard requirements are on the start codon itself and
    # on the mean of its context; demanding every position of a 30+ nt window be >=0.85 unpaired
    # is not achievable in a real molecule, so it is a soft push instead.
    start_unpaired_min_prob: float = 0.70     # min p(unpaired) across the AUG triplet (hard)
    start_context_mean_min: float = 0.60      # mean p(unpaired) over the start context (hard)
    # RNAplfold settings for the fast accessibility screen. 200/150 is ~15x faster than 600/300
    # and ranks candidates about as well (benchmarked against the exact circular PF: mean
    # absolute error 0.17 vs 0.12 on the AUG, 0.03 on the context mean). The exact partition
    # function re-scores the elites and decides feasibility, so this only has to rank.
    access_window: int = 200
    access_span: int = 150
    min_paired_fraction: float = 0.10

    weights: dict = field(default_factory=lambda: dict(DEFAULT_WEIGHTS))

    # search
    population: int = 40
    generations_per_round: int = 12
    rounds: int = 4
    elite_k: int = 6
    top_n: int = 8
    mutation_rate: float = 0.02
    polish_codons: int = 8
    repair: bool = True               # directed repair of GC-window / motif violations
    seed: int = 1
    workers: int = 0                  # 0 = all cores
    seed_sequences: list = field(default_factory=list)  # extra full-RNA seeds (e.g. from mRNAid)

    # DERNA (bioconda `derna`): Pareto MFE/CAI designer, used as a seed + comparison arm.
    # It optimises toward *minimal* MFE (maximally structured) with its built-in human table,
    # is O(n^3) single-threaded, and its custom-table option (-c) is unusable (see INSTALL.md).
    use_derna: bool = False
    derna_bin: str | None = None          # None = look for `derna` on PATH
    derna_lambdas: list = field(default_factory=lambda: [0.0, 0.3])  # 0 = max CAI, 1 = min MFE
    derna_max_aa: int = 400               # skip DERNA above this protein length (runtime)
    derna_timeout_s: int = 7200

    @classmethod
    def from_dict(cls, d: dict) -> "Params":
        names = {f.name for f in fields(cls)}
        unknown = set(d) - names
        if unknown:
            raise ValueError(f"unknown parameter(s): {sorted(unknown)}")
        d = dict(d)
        w = dict(DEFAULT_WEIGHTS)
        w.update(d.pop("weights", {}) or {})
        unk_w = set(w) - set(DEFAULT_WEIGHTS)
        if unk_w:
            raise ValueError(f"unknown weight(s): {sorted(unk_w)}")
        p = cls(**d)
        p.weights = w
        p.validate()
        return p

    @classmethod
    def from_json(cls, path: str) -> "Params":
        with open(path) as fh:
            return cls.from_dict(json.load(fh))

    def to_dict(self) -> dict:
        return asdict(self)

    def validate(self) -> None:
        if not self.protein and (self.cds_start is None or self.cds_end is None):
            raise ValueError(
                "Cannot locate the CDS: give the protein (amino-acid sequence), or give both "
                "cds_start and cds_end (0-based, end exclusive). The protein is also what "
                "guarantees the redesign stays synonymous, so prefer it when you have it.")
        if self.primary_species not in SPECIES or self.secondary_species not in SPECIES:
            raise ValueError(f"species must be in {SPECIES}")
        if self.primary_species == self.secondary_species:
            raise ValueError("primary and secondary species must differ")
        if not 0 < self.cai_target <= 1.0:
            raise ValueError("cai_target must be in (0,1]")
        if self.cai_floor > self.cai_target:
            raise ValueError("cai_floor must be <= cai_target")
        if not 0 < self.primary_weight < 1:
            raise ValueError("primary_weight must be in (0,1)")
        if self.structure_goal not in ("minimize_pairing", "minimize_mfe", "target"):
            raise ValueError("structure_goal must be minimize_pairing | minimize_mfe | target")
        if self.dsrna_min_bp < 5:
            raise ValueError("dsrna_min_bp must be >= 5")
        if self.codon_metric not in ("cai", "tai"):
            raise ValueError("codon_metric must be 'cai' or 'tai'")
        if not 0 <= self.gc_min < self.gc_max <= 1:
            raise ValueError("need 0 <= gc_min < gc_max <= 1")
        if self.dangles not in (0, 1, 2, 3):
            raise ValueError("dangles must be 0..3")
        if self.stem_target >= self.stem_cap:
            raise ValueError("stem_target must be < stem_cap")
        if self.use_derna and any(not 0.0 <= x <= 1.0 for x in self.derna_lambdas):
            raise ValueError("derna_lambdas must be in [0,1]")
        if self.population < 4 or self.rounds < 1 or self.generations_per_round < 1:
            raise ValueError("population>=4, rounds>=1, generations_per_round>=1 required")
        for r in self.protected:
            if not {"start", "end"} <= set(r) or r["start"] >= r["end"]:
                raise ValueError(f"bad protected region: {r}")
