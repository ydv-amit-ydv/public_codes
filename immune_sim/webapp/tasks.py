#!/usr/bin/env python3
"""The operations the GUI may run, and what each one needs and produces.

This is an ALLOWLIST. The browser posts a task name and a dict of arguments;
every argument is coerced and clamped here, and anything not listed is
rejected. No part of the request reaches a shell.

Each task also declares `requires` and `produces` so the GUI can tell the
user what is missing before they press a button, and can mark a result as
STALE when something it was computed from has changed since. That staleness
check is the whole reason this file knows about file paths: a GUI that shows
a validation number computed from a model that has since been retrained is
worse than one that shows nothing.
"""
import time
from pathlib import Path

import jobs

ROOT = Path(__file__).resolve().parent.parent
PY = "python3"


def _int(args, key, default, lo, hi):
    try:
        v = int(args.get(key, default))
    except (TypeError, ValueError):
        raise ValueError(f"{key} must be a whole number")
    if not lo <= v <= hi:
        raise ValueError(f"{key} must be between {lo} and {hi}")
    return v


def _choice(args, key, default, allowed):
    v = str(args.get(key, default))
    if v not in allowed:
        raise ValueError(f"{key} must be one of {', '.join(allowed)}")
    return v


def _flag(args, key):
    return bool(args.get(key))


def _params_file(args, key="params", default="models/calibrated_wetlab.json"):
    """A parameter file, restricted to models/*.json inside the project."""
    v = str(args.get(key) or default)
    p = (ROOT / v).resolve()
    if p.parent != (ROOT / "models").resolve() or p.suffix != ".json":
        raise ValueError("parameter file must be a models/*.json in this project")
    if not p.exists():
        raise ValueError(f"{v} does not exist")
    return str(p.relative_to(ROOT))


# ---------------------------------------------------------------------------
# argv builders
# ---------------------------------------------------------------------------
def _harvest_iedb(a):
    return [PY, "-u", "harvest_iedb.py"], None


def _harvest_cmipb(a):
    return [PY, "-u", "harvest_cmipb.py"], None


def _harvest_eluted(a):
    return [PY, "-u", "harvest_eluted.py",
            "--decoy-ratio", str(_int(a, "decoy_ratio", 19, 1, 99)),
            "--max-records", str(_int(a, "max_records", 60000, 500, 500000)),
            "--min-per-allele", str(_int(a, "min_per_allele", 120, 10, 5000))], None


def _harvest_trials(a):
    return [PY, "-u", "harvest_trials.py",
            "--max-studies", str(_int(a, "max_studies", 1500, 50, 20000)),
            "--page-size", str(_int(a, "page_size", 50, 10, 200))], None


def _train_allele(a):
    return [PY, "-u", "model.py"], None


def _train_pan(a):
    argv = [PY, "-u", "-c", "import model_pan; model_pan.main()",
            "--realign", str(_int(a, "realign", 4, 1, 12)),
            "--classes", _choice(a, "classes", "I,II", ("I", "II", "I,II"))]
    if _flag(a, "quick"):
        argv.append("--quick")
    return argv, None


def _train_immuno(a):
    return [PY, "-u", "immunogenicity.py"], None


def _calibrate(a):
    """Calibration NEVER writes the shipped parameter set.

    It writes a timestamped candidate instead. Overwriting
    models/calibrated_wetlab.json from a button press is how you lose a
    converged fit to a 20-second accident: the search writes a checkpoint on
    every improvement, so even a run that is still at score 9.9 has already
    replaced a run that reached 2.5e-05. Promotion is a separate, explicit
    action, after the candidate has been scored on the held-out split.
    """
    iters = _int(a, "iters", 2000, 50, 200000)
    out = f"models/calibrated_candidate_{time.strftime('%Y%m%d_%H%M%S')}.json"
    argv = [PY, "-u", "calibrate.py", "--wetlab",
            "--iters", str(iters),
            "--workers", str(_int(a, "workers", 0, 0, 64)),
            "--restarts", str(_int(a, "restarts", 1, 1, 20)),
            "--seed", str(_int(a, "seed", 0, 0, 10**6)),
            "--out", out]
    if a.get("init"):
        argv += ["--init", _params_file(a, "init")]
    return argv, jobs.calibrate_parser(iters)


def _benchmarks(a):
    return [PY, "-u", "test_platforms.py"], jobs.counting_parser(
        [("PASS", "running checks"), ("passed,", "done")])


def _wetlab(a):
    return [PY, "-u", "wetlab_validation.py",
            "--split", _choice(a, "split", "test", ("test", "train", "all")),
            "--params", _params_file(a)], None


def _lopo(a):
    return [PY, "-u", "external_validation.py"], None


def _lopo_pan(a):
    argv = [PY, "-u", "external_validation_pan.py",
            "--realign", str(_int(a, "realign", 4, 1, 12))]
    if _flag(a, "quick"):
        argv.append("--quick")
    return argv, None


def _eluted(a):
    return [PY, "-u", "validate_eluted.py",
            "--max-per-allele", str(_int(a, "max_per_allele", 60000, 100, 500000))], None


def _trials(a):
    return [PY, "-u", "validate_trials.py"], None


def _scorecard(a):
    return [PY, "-u", "scorecard.py",
            "--robustness", str(_int(a, "robustness", 100, 0, 2000))], None


# ---------------------------------------------------------------------------
# the registry
# ---------------------------------------------------------------------------
TASKS = {
    # -- data ---------------------------------------------------------------
    "harvest_iedb": dict(
        group="data", label="Harvest IEDB binding + immunogenicity",
        build=_harvest_iedb, requires=[], network=True, est="3-10 min",
        produces=["data/binding_affinity.csv", "data/tcell_immunogenicity.csv"],
        about="Measured peptide-MHC affinities and T cell assay outcomes. "
              "~23 MB. IEDB has its own citation and licence terms, which is "
              "why this is downloaded rather than redistributed."),
    "harvest_cmipb": dict(
        group="data", label="Harvest CMI-PB human boost data",
        build=_harvest_cmipb, requires=[], network=True, est="under 1 min",
        produces=["data/cmipb_ab_titer.csv", "data/cmipb_cell_frequency.csv"],
        about="166 adults given a Tdap boost, antibody titers at days "
              "0/1/3/7/14/30. This is the dataset the dynamic model is "
              "calibrated and held-out tested against."),
    "harvest_eluted": dict(
        fields=[
            dict(name='decoy_ratio', label='Decoys per positive', type='int', default=19, min=1, max=99, help='Proteome-sampled non-ligands per observed ligand. 19 gives 5% prevalence. Changing it changes ROC-AUC very little and PPV a lot, which is why both are reported.'),
            dict(name='max_records', label='Max records to fetch', type='int', default=60000, min=500, max=500000, help='Cap on IEDB rows pulled.'),
            dict(name='min_per_allele', label='Min ligands per allele', type='int', default=120, min=10, max=5000, help='Alleles with fewer are dropped.'),
        ],
        group="data", label="Harvest mass-spec eluted ligands",
        build=_harvest_eluted, requires=[], network=True, est="5-15 min",
        produces=["data/eluted_ligands.csv"],
        about="What cells actually present, from mass spectrometry, plus "
              "proteome-sampled decoys. An independent assay modality: no "
              "affinity measurement anywhere in it, so a model trained on "
              "IC50 has never seen these labels."),
    "harvest_trials": dict(
        fields=[
            dict(name='max_studies', label='Max studies', type='int', default=1500, min=50, max=20000, help='Registered trials with posted results to scan.'),
            dict(name='page_size', label='Page size', type='int', default=50, min=10, max=200, help='ClinicalTrials.gov API page size.'),
        ],
        group="data", label="Harvest ClinicalTrials.gov fold-rises",
        build=_harvest_trials, requires=[], network=True, est="3-10 min",
        produces=["data/trial_foldrise.csv"],
        about="Posted antibody fold-rise results from registered vaccine "
              "trials, with platform and adjuvant parsed from the arm "
              "descriptions. Observational: platform is not randomised."),

    # -- training -----------------------------------------------------------
    "train_allele": dict(
        group="train", label="Train per-allele binding models (model.py)",
        build=_train_allele, requires=["data/binding_affinity.csv"],
        est="about 10 s", produces=["models/report.json"],
        about="One gradient-boosted model per allele on fixed-length "
              "peptides. Fast, and the baseline every bigger model is "
              "measured against."),
    "train_pan": dict(
        fields=[
            dict(name='realign', label='Realignment passes', type='int', default=4, min=1, max=12, help='Iterations of binding-register re-estimation. Passes 3-4 changed held-out AUC by under 0.003, so 2 is usually enough.'),
            dict(name='classes', label='MHC class', type='choice', default='I,II', options=['I', 'II', 'I,II'], help='Train one class or both.'),
            dict(name='quick', label='Quick mode', type='bool', default=False, help='Small models, for checking the plumbing. The metrics it prints are NOT the real ones.'),
        ],
        group="train", label="Train pan-allele core-scan ensemble (model_pan.py)",
        build=_train_pan, requires=["data/binding_affinity.csv"],
        est="about 6 min at realign 4", produces=["models/pan_report.json"],
        about="One model across all alleles, with the 9-mer binding register "
              "chosen by iterative realignment and a hierarchical allele "
              "encoding so an untrained allele falls back to its locus "
              "rather than being scored as a different allele. Use --quick "
              "first to check the plumbing."),
    "train_immuno": dict(
        group="train", label="Train immunogenicity model",
        build=_train_immuno, requires=["data/tcell_immunogenicity.csv"],
        est="1-3 min", produces=["models/immuno_report.json"],
        about="Whether a bound peptide actually elicits a T cell response, "
              "which is a different question from whether it binds."),

    # -- calibration --------------------------------------------------------
    "calibrate": dict(
        fields=[
            dict(name='iters', label='Evaluations', type='int', default=2000, min=50, max=200000, help='Each evaluation solves the ODE for every vaccine configuration. 12,000 with 4 restarts converged this model.'),
            dict(name='workers', label='Worker processes', type='int', default=0, min=0, max=64, help='0 uses every core. Leave one free if you want the machine to stay responsive.'),
            dict(name='restarts', label='Restarts', type='int', default=1, min=1, max=20, help='Re-heats the search from the best point found, to escape a local optimum.'),
            dict(name='seed', label='Seed', type='int', default=0, min=0, max=1000000, help='Changing it changes the search path, not the targets.'),
            dict(name='init', label='Start from', type='params', default='models/calibrated_wetlab.json', help='Resume from an existing fit. The shipped set is the sensible start: from the literature defaults a short run lands far from it.'),
        ],
        group="calibrate", label="Calibrate the dynamic model against CMI-PB",
        build=_calibrate, requires=["data/cmipb_ab_titer.csv"],
        est="minutes to hours, set by iterations",
        produces=[],
        about="Parallel annealed search over 36 ODE parameters against 13 "
              "measured targets and 7 cross-configuration constraints. "
              "Writes a checkpoint on every improvement, so cancelling keeps "
              "the best parameters found so far. The result is saved as a "
              "timestamped CANDIDATE and is not used by anything until you "
              "promote it below \u2014 a short run reaches a far worse score "
              "than the shipped fit, and overwriting on a button press would "
              "lose that fit silently."),

    # -- validation ---------------------------------------------------------
    "benchmarks": dict(
        group="validate", label="Platform benchmarks (13 checks)",
        build=_benchmarks, requires=[], est="1-3 min", produces=[],
        about="Documented wet-lab directions the model must reproduce: mRNA "
              "beats protein on CD8, alum is Th2, a boost is faster than a "
              "prime, cell counts stay physiological. No data download "
              "needed."),
    "wetlab": dict(
        fields=[
            dict(name='split', label='Subject split', type='choice', default='test', options=['test', 'train', 'all'], help='test = subjects never used in calibration. This is the only one that is an out-of-sample claim.'),
            dict(name='params', label='Parameters', type='params', default='models/calibrated_wetlab.json', help='Which fitted parameter set to score.'),
        ],
        group="validate", label="Boost kinetics vs held-out human subjects",
        build=_wetlab, requires=["data/cmipb_ab_titer.csv"], est="1-2 min",
        produces=["results/wetlab_validation_test.csv",
                  "results/wetlab_validation_test.json"],
        about="The out-of-sample test that matters for the dynamic model. "
              "Subjects are split by ID; calibration sees one half and this "
              "scores the other. Includes negative-control antigens that are "
              "not in the vaccine."),
    "lopo": dict(
        group="validate", label="Leave-one-pathogen-out (model.py)",
        build=_lopo, requires=["data/binding_affinity.csv"], est="10-30 min",
        produces=["results/external_validation.csv"],
        about="Retrains with one organism entirely removed, then ranks its "
              "epitopes. The real situation when a novel pathogen arrives."),
    "lopo_pan": dict(
        fields=[
            dict(name='realign', label='Realignment passes', type='int', default=4, min=1, max=12, help='Must match how the shipped model was trained to be comparable.'),
            dict(name='quick', label='Quick mode', type='bool', default=False, help='Plumbing check only.'),
        ],
        group="validate", label="Leave-one-pathogen-out (pan vs per-allele)",
        build=_lopo_pan, requires=["data/binding_affinity.csv"],
        est="about 90 min, 16 folds",
        produces=["results/external_validation_pan.csv"],
        about="The same test for the pan-allele model, paired against "
              "model.py on identical rows with a DeLong test and "
              "Benjamini-Hochberg correction, plus an amino-acid-composition "
              "baseline so a model that only learned composition cannot pass."),
    "eluted": dict(
        fields=[
            dict(name='max_per_allele', label='Max rows per allele', type='int', default=60000, min=100, max=500000, help='Caps the biggest alleles so one of them cannot dominate the pooled metric.'),
        ],
        group="validate", label="Validate against eluted ligands",
        build=_eluted, requires=["data/eluted_ligands.csv"], est="about 10 min",
        produces=["results/eluted_validation.csv"],
        about="An independent assay modality, stratified by whether the "
              "allele was trained on, matched only at its locus or first "
              "field, or never seen. The stratification is the point: it is "
              "how the allele-encoding bug was found."),
    "trials": dict(
        group="validate", label="Compare against clinical trial fold-rises",
        build=_trials, requires=["data/trial_foldrise.csv"], est="1-2 min",
        produces=["results/trial_validation.csv"],
        about="Clustered by trial, not by measurement, because one trial "
              "posting many rows is pseudo-replication. Reports the "
              "confound first: platform is not randomised across trials."),
    "scorecard": dict(
        fields=[
            dict(name='robustness', label='Parameter draws', type='int', default=100, min=0, max=2000, help='Resamples every load-bearing parameter by +-25%. 0 skips resampling and reports grades only.'),
        ],
        group="validate", label="Evidence scorecard + robustness resampling",
        build=_scorecard, requires=[], est="2-5 min at 100 draws",
        produces=["results/scorecard.csv"],
        about="Every claim the project makes, graded A-D, with the share of "
              "plausible parameter sets that reproduce each direction. That "
              "share is NOT a probability that the biology is true."),
}

GROUPS = [
    ("data", "Data", "Nothing is redistributed in the clone. Each dataset "
                     "is downloaded from its source, which is also what "
                     "keeps the licence and citation with the provider."),
    ("train", "Models", "Train the epitope predictors. The dynamic model "
                        "needs none of this: it ships calibrated."),
    ("calibrate", "Calibration", "Refit the dynamic model's ODE parameters "
                                 "against measured human data."),
    ("validate", "Validation", "Every number in VALIDATION.md, re-runnable. "
                               "These can fail, and that is the point."),
]


# ---------------------------------------------------------------------------
# what each output was computed FROM, so the GUI can mark it stale
# ---------------------------------------------------------------------------
# A validation number computed from a model that has since been retrained is
# worse than no number at all, because it looks current. The GUI compares
# modification times against this graph and says so in the open.
DEPENDS = {
    "models/report.json":        ["data/binding_affinity.csv"],
    "models/pan_report.json":    ["data/binding_affinity.csv"],
    "models/immuno_report.json": ["data/tcell_immunogenicity.csv"],
    "models/calibrated_wetlab.json": ["data/cmipb_ab_titer.csv"],
    "results/external_validation.csv":     ["data/binding_affinity.csv"],
    "results/external_validation_pan.csv": ["data/binding_affinity.csv"],
    "results/eluted_validation.csv": ["data/eluted_ligands.csv",
                                      "models/pan_report.json"],
    "results/wetlab_validation_test.csv": ["data/cmipb_ab_titer.csv",
                                           "models/calibrated_wetlab.json"],
    "results/trial_validation.csv": ["data/trial_foldrise.csv",
                                     "models/calibrated_wetlab.json"],
    "results/scorecard.csv": ["models/calibrated_wetlab.json"],
}

# Everything the inventory shows, in display order.
INVENTORY = [
    ("data/binding_affinity.csv",      "IEDB binding affinities"),
    ("data/tcell_immunogenicity.csv",  "IEDB T cell assays"),
    ("data/cmipb_ab_titer.csv",        "CMI-PB antibody titers"),
    ("data/cmipb_cell_frequency.csv",  "CMI-PB cell frequencies"),
    ("data/eluted_ligands.csv",        "Mass-spec eluted ligands"),
    ("data/trial_foldrise.csv",        "Clinical trial fold-rises"),
    ("models/report.json",             "Per-allele model report"),
    ("models/pan_report.json",         "Pan-allele model report"),
    ("models/immuno_report.json",      "Immunogenicity model report"),
    ("models/calibrated_wetlab.json",  "Calibrated ODE parameters"),
    ("results/external_validation.csv",     "LOPO, per-allele"),
    ("results/external_validation_pan.csv", "LOPO, pan-allele"),
    ("results/eluted_validation.csv",       "Eluted-ligand validation"),
    ("results/wetlab_validation_test.csv",  "Held-out human subjects"),
    ("results/trial_validation.csv",        "Clinical trial comparison"),
    ("results/scorecard.csv",               "Evidence scorecard"),
]


def stale_because(rel):
    """Inputs modified more recently than this output, newest first."""
    out = ROOT / rel
    if not out.exists():
        return []
    t = out.stat().st_mtime
    newer = [(d, (ROOT / d).stat().st_mtime)
             for d in DEPENDS.get(rel, [])
             if (ROOT / d).exists() and (ROOT / d).stat().st_mtime > t + 1]
    return [d for d, _ in sorted(newer, key=lambda x: -x[1])]


def build(name, args):
    """Return (label, argv, produces, parser) for an allowlisted task."""
    spec = TASKS.get(name)
    if spec is None:
        raise ValueError(f"unknown task {name!r}")
    argv, parser = spec["build"](args or {})
    return spec["label"], argv, spec["produces"], parser


def missing(name):
    """Inputs a task needs that are not present."""
    spec = TASKS.get(name)
    if spec is None:
        return []
    return [r for r in spec["requires"] if not (ROOT / r).exists()]
