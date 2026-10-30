"""Unit tests: run with `pytest -q` from the package root."""
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from circopt.codon import AA_TO_CODONS, CODON_TO_AA, CodonModel, translate  # noqa: E402
from circopt.config import Params  # noqa: E402
from circopt.constraints import Rules, gc_stats, inverted_repeat_count, revcomp  # noqa: E402
from circopt.design import Design  # noqa: E402
from circopt.score import score  # noqa: E402
from circopt.structure import (circular_accessibility, fold, helices,  # noqa: E402
                               longest_stem, pair_table)

PROTEIN = "MFVNQHLCGSHLVEALYLVCGERGFFYTPKTGIVEQCCTSICSLYQLENYCN"


def make_params(**kw):
    cds = "".join(sorted(AA_TO_CODONS[a])[0] for a in PROTEIN) + "UAA"
    ires = "GGGAGACCCAAGCUGGCUAGCGUUUAAACUUAAGCUUGGUACCGAGCUCGGAUCC"
    rna = ires + "AAACAAACAAA" + cds + "AAAGAAACAAA"
    base = dict(name="t", protein=PROTEIN, rna=rna,
                protected=[{"name": "ires", "start": 0, "end": len(ires)}],
                population=6, generations_per_round=2, rounds=1, elite_k=2, top_n=2,
                polish_codons=0, workers=1, inner_span=0, stress_test=False,
                linear_compare=False, access_window=200, access_span=120)
    base.update(kw)
    return Params.from_dict(base)


# ---- codon layer ----------------------------------------------------------
def test_genetic_code_and_translate():
    assert CODON_TO_AA["AUG"] == "M" and CODON_TO_AA["UGA"] == "*"
    assert translate("AUGUUUUAA") == "MF*"
    assert len(CODON_TO_AA) == 64


def test_cai_weights_and_index():
    m = CodonModel()
    for sp in ("human", "mouse"):
        assert max(m.w[sp][c] for c in AA_TO_CODONS["L"]) == pytest.approx(1.0)
        # the most-used Leu codon in both species is CUG
        assert m.w[sp]["CUG"] == pytest.approx(1.0)
    best = [max(AA_TO_CODONS[a], key=lambda c: m.w["human"][c]) for a in "LVSPTAGR"]
    assert m.index(best, "human") == pytest.approx(1.0)
    worst = [min(AA_TO_CODONS[a], key=lambda c: m.w["human"][c]) for a in "LVSPTAGR"]
    assert m.index(worst, "human") < m.index(best, "human")


def test_min_codon_freq_filters_rare_codons():
    m = CodonModel(min_freq=0.15)
    assert "UUA" not in m.options("L")        # ~0.08 in both species
    assert "CUG" in m.options("L")
    assert m.options("M") == ["AUG"]          # never empty for single-codon aa


def test_tai_requires_trna_file():
    with pytest.raises(ValueError):
        CodonModel(metric="tai")


# ---- structure layer ------------------------------------------------------
def test_pair_table_and_helices():
    assert pair_table("((..))") == [5, 4, -1, -1, 1, 0]
    assert [h[2] for h in helices(pair_table("((((....))))"))] == [4]
    assert sorted(h[2] for h in helices(pair_table("(((...)))..((..))"))) == [2, 3]
    assert longest_stem(pair_table(".....")) == 0
    for bad in ("(", ")", "(()"):
        with pytest.raises(ValueError):
            pair_table(bad)


def test_fold_is_circular_and_metrics_present():
    p = make_params()
    d = Design(p)
    m = fold(d.rna, p, d.cds, circ=True, full=True, span=0)
    assert m.mfe < 0 and 0 <= m.paired_fraction <= 1
    assert m.structure.count("(") == m.structure.count(")")
    assert m.efe <= m.mfe + 1e-6              # ensemble energy is never above the MFE
    assert len(m.unpaired_prob) == len(m.entropy) == d.n
    assert m.start_codon_unpaired is not None
    # Circular mode really is in effect: closing the circle turns the exterior loop into a
    # real loop, so the energy differs from the linear fold of the same sequence. (Direction
    # is not asserted -- the circular MFE can be either side of the linear one depending on
    # the construct; here it is less stable, -38.2 vs -44.0.)
    lin = fold(d.rna, p, d.cds, circ=False, full=False, span=0)
    assert lin.mfe < 0 and m.mfe < 0
    assert abs(lin.mfe - m.mfe) > 1e-6


def test_dsrna_counts_only_long_duplexes():
    p = make_params(dsrna_min_bp=8)
    d = Design(p)
    # a planted 12 bp perfect duplex
    arm = "GGCCGGCCGGCC"
    seq = arm + "AAAAAAAA" + revcomp(arm) + "AAAAAAAAAAAA"
    m = fold(seq, p, (0, 12), circ=True, full=False, span=0)
    assert m.longest_stem >= 8
    assert m.dsrna_nt >= 16 and m.dsrna_helices >= 1


def test_accessibility_profile_shape():
    p = make_params()
    d = Design(p)
    acc = circular_accessibility(d.rna, p)
    assert len(acc) == d.n
    assert all(0.0 <= a <= 1.0 for a in acc)
    # an A-rich stretch is more accessible than a GC-rich one
    assert sum(acc[55:66]) / 11 > 0.5


# ---- constraints ----------------------------------------------------------
def test_gc_stats_circular_window():
    assert gc_stats("GCGC", 2)["gc"] == 1.0
    s = gc_stats("GGGGAAAA", 4)
    assert s["gc"] == 0.5 and s["gc_win_max"] == 1.0 and s["gc_win_min"] == 0.0
    # the window wraps: GC-rich block split across the origin is still found
    s2 = gc_stats("GGAAAAGG", 4)
    assert s2["gc_win_max"] == 1.0


def test_gc_window_ignores_unfixable_protected_windows():
    # a GC block that is entirely protected must not be blamed on the design
    seq = "G" * 10 + "A" * 10
    mask = [False] * 10 + [True] * 10
    st = gc_stats(seq, 5, mask)
    assert st["gc_win_max"] == 1.0                 # the block is really there
    assert st["gc_win_max_free"] < 1.0             # but no designable window reaches it
    assert gc_stats(seq, 5, [True] * 20)["gc_win_max_free"] == 1.0


def test_gc_violation_is_reported(tmp_path):
    p = make_params(gc_min=0.0, gc_max=0.01)       # impossible on purpose
    d = Design(p)
    m = fold(d.rna, p, d.cds, circ=True, full=False, span=0)
    s = score(d, d.genome_of(), m, Rules(p))
    assert not s.feasible
    assert any(v.startswith("GC=") for v in s.violations)


def test_rules_detect_motifs_and_wrap():
    p = make_params(max_poly_u=3, extra_forbidden=["GAAUUC"])
    r = Rules(p)
    assert "polyU" in r.counts("AAAUUUUAAA")
    assert "polyU" not in r.counts("AAAUUUAAA")
    assert "motif:GAAUUC" in r.counts("AAGAAUUCAA")
    # a motif spanning the circular junction is still caught
    assert "motif:GAAUUC" in r.counts("UUCAAAAAAAAAGAA")
    assert "RE:EcoRI" in r.counts("AAAGAAUUCAAA")


def test_re_sites_check_both_strands():
    p = make_params(re_sites={"BsaI": "GGTCTC"})
    r = Rules(p)
    assert "RE:BsaI" in r.counts("AAGGUCUCAA")
    assert "RE:BsaI" in r.counts("AAGAGACCAA")      # reverse complement


def test_fixable_ignores_protected_violations():
    p = make_params(extra_forbidden=["GCUGGCUAGC"])   # sits inside the protected IRES
    d = Design(p)
    r = Rules(p)
    assert "motif:GCUGGCUAGC" in r.counts(d.rna)
    assert "motif:GCUGGCUAGC" not in r.fixable(d.rna, d.free_mask)


def test_inverted_repeat_detection():
    arm = "GGCCAAGGCCA"
    seq = arm + "AAAAAAAAAA" + revcomp(arm) + "AAAA"
    mask = [True] * len(seq)
    assert inverted_repeat_count(seq, 9, mask) > 0
    assert inverted_repeat_count("A" * 60, 9, [True] * 60) == 0


# ---- design layer ---------------------------------------------------------
def test_design_locates_cds_and_respects_protection():
    p = make_params()
    d = Design(p)
    assert d.rna[d.cds_start:d.cds_start + 3] == "AUG"
    assert d.protein.rstrip("*") == PROTEIN
    assert all(d.locked[i] for i in range(0, 55))            # the IRES
    assert d.build(d.genome_of()) == d.rna                   # round trip
    # no free slot overlaps a protected position
    for pos in d.slots:
        assert not any(d.locked[pos + k] for k in range(3))


def test_build_keeps_protein_identical():
    p = make_params()
    d = Design(p)
    import random
    rng = random.Random(0)
    for _ in range(25):
        g = [rng.randrange(len(o)) for o in d.options]
        seq = d.build(g)
        assert translate(seq[d.cds_start:d.cds_end]) == d.protein
        assert seq[:d.cds_start] == d.rna[:d.cds_start]       # scaffold untouched
        assert seq[d.cds_end:] == d.rna[d.cds_end:]
        assert len(seq) == d.n


def test_protein_mismatch_is_rejected():
    p = make_params(protein="MFVNQHLCGSHLVEALYLVCGERGFFYTPKTGIVEQCCTSICSLYQLENYCX")
    with pytest.raises(ValueError):
        Design(p)


def test_bad_rna_rejected():
    with pytest.raises(ValueError):
        Design(make_params(rna="AUGCNNNUAA", protein="M", cds_start=0, cds_end=6))


def test_cds_must_be_locatable():
    """The web form used to accept a blank protein and then crash inside the worker thread."""
    cds = "".join(sorted(AA_TO_CODONS[a])[0] for a in PROTEIN) + "UAA"
    rna = "GGGAGACCCAAGCUGGCUAGC" + cds
    with pytest.raises(ValueError, match="locate the CDS"):
        Params.from_dict({"name": "t", "rna": rna})          # no protein, no coordinates
    # coordinates alone are enough
    p = Params.from_dict({"name": "t", "rna": rna, "cds_start": 21, "cds_end": 21 + len(cds)})
    d = Design(p)
    assert d.protein.rstrip("*") == PROTEIN
    # protein alone is enough
    assert Design(Params.from_dict({"name": "t", "rna": rna, "protein": PROTEIN})).cds[0] == 21


def test_orf_suggestion_reports_real_coordinates():
    from circopt.app import _suggest_orf
    cds = "".join(sorted(AA_TO_CODONS[a])[0] for a in PROTEIN) + "UAA"
    rna = "GGGAGACCCAAGCUGGCUAGC" + cds
    msg = _suggest_orf(rna)
    assert "cds_start=21" in msg and f"cds_end={21 + len(cds)}" in msg
    assert _suggest_orf("not a sequence") == ""
    assert _suggest_orf("ACGUACGU") == ""                    # no ORF -> no bogus suggestion


def test_config_validation():
    with pytest.raises(ValueError):
        Params.from_dict({"primary_species": "rat", "secondary_species": "human"})
    with pytest.raises(ValueError):
        Params.from_dict({"primary_species": "human", "secondary_species": "human"})
    with pytest.raises(ValueError):
        Params.from_dict({"gc_min": 0.8, "gc_max": 0.3})
    with pytest.raises(ValueError):
        Params.from_dict({"structure_goal": "nonsense"})
    with pytest.raises(ValueError):
        Params.from_dict({"unknown_key": 1})
    with pytest.raises(ValueError):
        Params.from_dict({"weights": {"not_a_weight": 1.0}})


# ---- scoring --------------------------------------------------------------
def test_score_is_sum_of_parts_and_flags_violations():
    p = make_params()
    d = Design(p)
    m = fold(d.rna, p, d.cds, circ=True, full=True, span=0)
    s = score(d, d.genome_of(), m, Rules(p))
    assert s.total == pytest.approx(sum(s.parts.values()))
    assert set(s.cai) == {"human", "mouse"}
    assert all(v >= 0 for v in s.parts.values())
    row = s.row(p)
    assert "cai_human" in row and "mfe_per_nt" in row


def test_structure_goal_changes_mfe_penalty_direction():
    d_min = Design(make_params(structure_goal="minimize_pairing"))
    p1 = d_min.p
    m = fold(d_min.rna, p1, d_min.cds, circ=True, full=False, span=0)
    s1 = score(d_min, d_min.genome_of(), m, Rules(p1))
    p2 = make_params(structure_goal="minimize_mfe", mfe_per_nt_target=-0.45)
    d2 = Design(p2)
    s2 = score(d2, d2.genome_of(), m, Rules(p2))
    # the same structure is penalised under one goal and not the other
    assert s1.parts["mfe"] > 0
    assert s2.parts["mfe"] != s1.parts["mfe"]


def test_structure_goals_rank_mfe_in_opposite_orders():
    """The MFE term must be monotonic and opposite under the two goals.

    Checked on the term itself, not on the winning variant: with default weights the codon
    and GC terms can dominate, so flipping the goal does not necessarily flip which candidate
    wins -- raise weights["mfe"] to explore the stability-first direction.
    """
    import random
    base = Design(make_params())
    rng = random.Random(21)
    genomes = [base.genome_of()] + [[rng.randrange(len(o)) for o in base.options]
                                    for _ in range(7)]
    folded = []
    for g in genomes:
        seq = base.build(g)
        folded.append((g, seq, fold(seq, base.p, base.cds, circ=True, full=False, span=0)))
    pens = {}
    for goal in ("minimize_pairing", "minimize_mfe"):
        p = make_params(structure_goal=goal, mfe_per_nt_target=-0.45)
        d, r = Design(p), Rules(p)
        pens[goal] = [score(d, g, m, r, seq=s).parts["mfe"] for g, s, m in folded]
    mfes = [m.mfe_per_nt for _, _, m in folded]
    least = max(range(len(mfes)), key=lambda i: mfes[i])     # least negative = least paired
    most = min(range(len(mfes)), key=lambda i: mfes[i])
    assert pens["minimize_pairing"][least] == min(pens["minimize_pairing"])
    assert pens["minimize_mfe"][most] == min(pens["minimize_mfe"])


def test_dual_species_geometric_mean_between_the_two():
    p = make_params()
    d = Design(p)
    m = fold(d.rna, p, d.cds, circ=True, full=False, span=0)
    s = score(d, d.genome_of(), m, Rules(p))
    assert min(s.cai.values()) <= s.cai_gm <= max(s.cai.values())


# ---- DERNA integration (graft only; does not need the binary) --------------
def test_derna_graft_preserves_protein_and_scaffold():
    from circopt.derna import _parse, graft
    p = make_params()
    d = Design(p)
    prot = d.protein.rstrip("*")
    model = d.model
    # a max-CAI CDS without a stop codon, which is what DERNA emits
    cds = "".join(max(AA_TO_CODONS[a], key=lambda c: model.w["human"][c]) for a in prot)
    full, adopted = graft(d, cds)
    assert adopted > 0
    assert translate(full[d.cds_start:d.cds_end]).rstrip("*") == prot
    assert full[:d.cds_start] == d.rna[:d.cds_start]
    assert full[d.cds_end:] == d.rna[d.cds_end:]
    assert len(full) == d.n
    # a CDS of the wrong length is rejected rather than silently grafted
    assert graft(d, cds[:-3]) == (None, 0)
    # undetermined codons (X) in DERNA's MFE-only output are skipped
    assert _parse("zuker rna: AUGCUGXUG\nzuker cai rna: AUGCUGCUG\n") == ["AUGCUGCUG"]


# ---- directed repair -------------------------------------------------------
def test_repair_clears_motifs_without_breaking_the_protein():
    from circopt.optimize import repair
    p = make_params(extra_forbidden=["GAAUUC"], max_poly_u=3)
    d = Design(p)
    r = Rules(p)
    g = repair(d, d.genome_of(), r)
    seq = d.build(g)
    assert sum(r.fixable(seq, d.free_mask).values()) <= \
           sum(r.fixable(d.rna, d.free_mask).values())
    assert translate(seq[d.cds_start:d.cds_end]) == d.protein     # still synonymous
    assert seq[:d.cds_start] == d.rna[:d.cds_start]               # scaffold untouched


def test_repair_pulls_gc_windows_into_band_and_respects_the_cai_floor():
    from circopt.optimize import gc_excess, repair
    p = make_params(gc_window=30, gc_window_max=0.55, gc_window_min=0.2, cai_floor=0.5)
    d = Design(p)
    r = Rules(p)
    model = d.model
    # start from the GC-heaviest synonymous choice available
    g0 = [max(range(len(o)), key=lambda i, o=o: sum(1 for b in o[i] if b in "GC"))
          for o in d.options]
    before = gc_excess(d.build(g0), p, d.free_mask)
    g1 = repair(d, g0, r)
    after = gc_excess(d.build(g1), p, d.free_mask)
    assert after < before, "repair must reduce total GC-window excess"
    weak = min(model.index(d.cds_codons(d.build(g1)), sp) for sp in ("human", "mouse"))
    start_weak = min(model.index(d.cds_codons(d.build(g0)), sp) for sp in ("human", "mouse"))
    assert weak >= min(start_weak, p.cai_floor) - 1e-9
    assert translate(d.build(g1)[d.cds_start:d.cds_end]) == d.protein


# ---- end to end -----------------------------------------------------------
def test_full_run_end_to_end(tmp_path):
    from circopt.pipeline import run_design
    p = make_params(rounds=2, generations_per_round=2, population=6, top_n=2, elite_k=2)
    res = run_design(p, str(tmp_path), log=lambda *a: None)
    best = res["best"]
    # the protein is unchanged -- the one inviolable constraint
    assert translate(best.seq[best_cds(res)[0]:best_cds(res)[1]]).rstrip("*") == PROTEIN
    assert len(best.seq) == len(p.rna)
    assert len(res["rounds"]) == 2
    assert res["rounds"][0].ensemble, "each cycle records a sampled ensemble"
    for f in ("ranked_csv", "fasta", "comparison_csv", "run_json", "ensembles"):
        assert os.path.exists(res["paths"][f])
    with open(res["paths"]["run_json"]) as fh:
        rec = json.load(fh)
    assert rec["best"]["seq"] == best.seq
    # the scaffold is byte-identical outside the CDS
    s, e = best_cds(res)
    assert best.seq[:s] == p.rna[:s] and best.seq[e:] == p.rna[e:]


def best_cds(res):
    return tuple(res["design"]["cds"])
