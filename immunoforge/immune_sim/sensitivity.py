#!/usr/bin/env python3
"""Is the simulation actually sensitive to the antigen you give it?

C-ImmSim was reported returning the same response for an antigenic protein
and for GAPDH under identical parameters. This model had the same defect
until central tolerance and absolute epitope weights went in. This module is
the standing check that it stays fixed, and that the sensitivity is real
rather than an artefact of one lucky comparison.

It runs four kinds of arm and reports them against each other:

  positive controls  proteins that should provoke a response
  negative controls  proteins and sequences that should not, each failing for
                     a DIFFERENT reason, so that passing them all means
                     something
  scrambles          the same amino acids in a random order. A model that
                     scored composition rather than sequence would give the
                     same answer here as for the real protein. Several
                     replicates, so the comparison is against a distribution
                     rather than one draw.
  mutants            single substitutions, separated into those that hit an
                     MHC anchor position of a dominant epitope and those that
                     do not. A model that understands binding should care
                     more about the first.

Replicates matter because the ODE is deterministic: two runs of the same
epitope set are identical to the last decimal. Any spread you see here comes
from the ANTIGEN, which is the point. The scramble and mutation arms are
randomised and repeated, so each gives a distribution that a real effect has
to stand outside of.
"""
import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import predict as PR                                         # noqa: E402
from dynamics import Sim, load_epitopes                      # noqa: E402

RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

# Anchor positions of the bound 9-mer core, counting from 1. P2 and P9 are
# the primary anchors for most class I alleles; P1, P4, P6 and P9 for class
# II. Mutating one of these should cost binding; mutating P5 usually should
# not, because it points at the T cell receptor rather than into the groove.
ANCHORS = {"I": (2, 9), "II": (1, 4, 6, 9)}
NON_ANCHOR = {"I": (5, 7), "II": (5, 8)}
AAS = "ACDEFGHIKLMNPQRSTVWY"


# ---------------------------------------------------------------------------
def simulate(seq, models, params, host_ix=None, platform="mrna_lnp",
             adjuvant="intrinsic", dose_ug=30.0, n_steps=400):
    """One antigen through the whole pipeline. Returns None if it cannot respond."""
    pred = PR.scan(seq, models)
    b = pred[pred.call != "non-binder"].sort_values("pred_IC50_nM").copy()
    if b.empty:
        return {"n_binders": 0, "peak_IgG": 0.0, "peak_CTL": 0.0,
                "reason": "no predicted binders"}
    if host_ix is not None:
        from common import hostself
        b["host_self"] = [hostself.check(p, host_ix)["Host_Self"]
                          for p in b.peptide]
    try:
        ep = load_epitopes(b)
    except ValueError as e:
        return {"n_binders": int(len(b)), "peak_IgG": 0.0, "peak_CTL": 0.0,
                "reason": str(e)[:90]}
    sim = Sim(ep, p=params, injections=[(2, None, adjuvant)],
              platform=platform, dose_ug=dose_ug)
    df = sim.run(n_steps=n_steps)
    return {"n_binders": int(len(b)), "n_clones": int(len(ep)),
            "sum_weight": round(float(ep.weight.sum()), 3),
            "peak_IgG": round(float(df.IgG.max()), 2),
            "peak_CTL": round(float(df.TC_active.max()), 1),
            "peak_B": round(float(df.B_active.max()), 1),
            "reason": ""}


def scramble(seq, rng):
    """Same residues, random order: destroys epitopes, keeps composition."""
    a = list(seq)
    rng.shuffle(a)
    return "".join(a)


def mutate_at(seq, pos0, rng):
    """Substitute one residue for a different one, chosen at random."""
    a = list(seq)
    old = a[pos0]
    a[pos0] = rng.choice([c for c in AAS if c != old])
    return "".join(a), old, a[pos0]


def dominant_epitopes(seq, models, host_ix, top=6):
    """The epitopes that actually drive the response, strongest first."""
    pred = PR.scan(seq, models)
    b = pred[pred.call != "non-binder"].sort_values("pred_IC50_nM").copy()
    if b.empty:
        return b
    if host_ix is not None:
        from common import hostself
        b["host_self"] = [hostself.check(p, host_ix)["Host_Self"]
                          for p in b.peptide]
        b = b[b.host_self != "identical"]
    return b.head(top)


# ---------------------------------------------------------------------------
# the arms
# ---------------------------------------------------------------------------
# Each negative control fails for a DIFFERENT reason. Passing all of them
# means more than passing one three times: a model could be right about self
# proteins and still be scoring composition, or right about composition and
# still give a response to a sequence with no binders at all.
POSITIVES = [
    ("P0DTC2", "SARS-CoV-2 spike", "a real vaccine antigen"),
    ("P01012", "ovalbumin", "the standard experimental immunogen"),
]
NEGATIVES = [
    ("P04406", "GAPDH (human)", "self: central tolerance should delete it"),
    ("P01579", "IFN-gamma (human)", "self, and a short one"),
]
SYNTHETIC_NEGATIVES = [
    ("poly_A", "A" * 300, "one residue repeated: no real epitope structure"),
    ("poly_GS", "GS" * 150, "a flexible linker, not an antigen"),
]


def epitope_level_anchor_test(models, host_ix, seqs, rng, n=40):
    """Does changing an anchor residue destroy binding, on the peptide itself?

    The whole-protein response saturates: a 1,273-residue antigen has
    thousands of binders, so losing one changes almost nothing, and a
    substitution test at that level measures noise. Measured on the PEPTIDE,
    the question is sharp - swap the residue that sits in the MHC groove and
    the predicted affinity should collapse; swap the one pointing at the T
    cell receptor and it should not.

    Class I only. A class II 15-mer's binding core sits somewhere inside it,
    so position 2 of the peptide is not position 2 of the core, and the first
    version of this test was mutating the wrong residues for class II.
    """
    rows = []
    for label, seq in seqs.items():
        dom = dominant_epitopes(seq, models, host_ix, top=60)
        dom = dom[(dom.mhc_class == "I") & (dom.peptide.str.len() == 9)]
        if dom.empty:
            continue
        for _, e in dom.head(n).iterrows():
            for kind, positions in (("anchor", ANCHORS["I"]),
                                    ("non_anchor", NON_ANCHOR["I"])):
                for pk in positions:
                    mutpep, old, new_aa = mutate_at(e.peptide, pk - 1, rng)
                    rows.append({"parent": label, "allele": e.allele,
                                 "peptide": e.peptide, "mutant": mutpep,
                                 "kind": kind, "pos": pk,
                                 "from_aa": old, "to_aa": new_aa,
                                 "base_nM": float(e.pred_IC50_nM)})
    if not rows:
        return {}
    df = pd.DataFrame(rows)
    # score every mutant peptide with the same per-allele models
    parts = []
    for allele, sub in df.groupby("allele"):
        sub = sub.copy()
        try:
            sub["mut_nM"] = PR.score_peptides(list(sub.mutant), allele, models)
        except ValueError:
            sub["mut_nM"] = np.nan
        parts.append(sub)
    df = pd.concat(parts, ignore_index=True)
    df = df.dropna(subset=["mut_nM"])
    df["fold_worse"] = df.mut_nM / df.base_nM.clip(lower=1e-6)
    a = df[df.kind == "anchor"].fold_worse
    nn = df[df.kind == "non_anchor"].fold_worse
    return {"anchor_loss": round(float(a.median()), 3) if len(a) else 0.0,
            "non_anchor_loss": round(float(nn.median()), 3) if len(nn) else 0.0,
            "n_anchor": int(len(a)), "n_non_anchor": int(len(nn)),
            "table": df.head(200).to_dict("records")}


def _load_user_antigen(spec):
    """A UniProt accession, a FASTA file, or a bare sequence."""
    spec = str(spec).strip()
    f = Path(spec)
    if f.exists():
        txt = f.read_text()
        seq = "".join(l.strip() for l in txt.splitlines()
                      if l.strip() and not l.startswith(">"))
        name = f.stem
    elif len(spec) <= 20 and spec.isalnum():
        name, seq = PR.fetch_uniprot(spec)
        name = name.split("|")[-1].strip()[:40] or spec
    else:
        seq, name = spec, "your antigen"
    seq = "".join(c for c in seq.upper() if c.isalpha())
    if len(seq) < 9:
        raise SystemExit(f"{spec!r}: no usable protein sequence found")
    return name, seq


def run(reps=8, seed=0, quick=False, out="sensitivity", antigen=None):
    import json as _json
    rng = np.random.default_rng(seed)
    params = _json.loads((HERE / "models" / "calibrated_wetlab.json")
                         .read_text())["params"]
    models = PR.load_models()
    host_ix = None
    try:
        sys.path.insert(0, str(HERE.parent))
        from common import hostself
        if hostself.available():
            host_ix = hostself.index()
    except Exception:                                        # noqa: BLE001
        pass
    if host_ix is None:
        print("  WARNING: no host proteome, so self proteins will NOT be "
              "recognised as self. The negative controls below are then not "
              "a fair test. Settings -> Host proteome downloads it.")

    seqs, rows = {}, []

    ref_epitopes = {}

    def arm(kind, label, seq, why, extra=None):
        r = simulate(seq, models, params, host_ix)
        if kind == "positive":
            ref_epitopes[label] = set(
                dominant_epitopes(seq, models, host_ix, top=400).peptide)
        if kind == "scramble" and extra and extra.get("parent") in ref_epitopes:
            got = set(dominant_epitopes(seq, models, host_ix, top=400).peptide)
            ref = ref_epitopes[extra["parent"]]
            r["epitope_overlap"] = (len(got & ref) / max(len(ref), 1))
        r.update({"arm": kind, "label": label, "why": why,
                  "length": len(seq)})
        r.update(extra or {})
        rows.append(r)
        note = r["reason"] or ""
        print(f"  {kind:11s} {label:26s} binders {r['n_binders']:5d}  "
              f"IgG {r['peak_IgG']:9.1f}  CTL {r['peak_CTL']:11.1f}"
              + (f"   {note[:44]}" if note else ""), flush=True)
        return r

    print("\nPOSITIVE CONTROLS - these should respond")
    if antigen:
        # Your own antigen, run through the same controls: its own shuffled
        # replicates and its own anchor substitutions, beside the standard
        # positives and negatives so there is something to compare it to.
        uname, useq = _load_user_antigen(antigen)
        seqs[uname] = useq
        arm("positive", uname, useq, "the antigen you asked about")
    for acc, label, why in POSITIVES:
        _n, s = PR.fetch_uniprot(acc)
        seqs[label] = s
        arm("positive", label, s, why)

    print("\nNEGATIVE CONTROLS - each should fail, each for its own reason")
    for acc, label, why in NEGATIVES:
        _n, s = PR.fetch_uniprot(acc)
        arm("negative", label, s, why)
    for label, s, why in SYNTHETIC_NEGATIVES:
        arm("negative", label, s, why)

    print(f"\nSCRAMBLES - same residues, random order, {reps} replicates each."
          "\n  A model scoring composition rather than sequence cannot tell "
          "these from the real protein.")
    for label, s in seqs.items():
        for i in range(reps if not quick else 2):
            arm("scramble", f"{label} scrambled #{i + 1}", scramble(s, rng),
                "composition preserved, epitopes destroyed",
                {"parent": label, "rep": i + 1})

    print(f"\nPOINT MUTATIONS - {reps} replicates per class."
          "\n  An anchor position holds the peptide in the MHC groove; P5 "
          "points at the T cell.\n  Losing an anchor should cost more than "
          "losing a non-anchor.")
    for label, s in seqs.items():
        dom = dominant_epitopes(s, models, host_ix, top=max(3, reps // 2))
        if dom.empty:
            continue
        for kind, positions in (("anchor", ANCHORS), ("non_anchor", NON_ANCHOR)):
            done = 0
            for _, e in dom.iterrows():
                if done >= (reps if not quick else 2):
                    break
                cls = e.mhc_class
                core_start = int(e.start)
                for pk in positions[cls]:
                    if done >= (reps if not quick else 2):
                        break
                    pos0 = core_start + pk - 1
                    if pos0 >= len(s):
                        continue
                    mut, old, new = mutate_at(s, pos0, rng)
                    arm(f"mut_{kind}", f"{label} {old}{pos0 + 1}{new} (P{pk})",
                        mut, f"{kind} position of a dominant epitope",
                        {"parent": label, "epitope": e.peptide,
                         "position": pos0 + 1, "anchor_pos": pk,
                         "from_aa": old, "to_aa": new})
                    done += 1
        for i in range(reps if not quick else 2):
            pos0 = int(rng.integers(0, len(s)))
            mut, old, new = mutate_at(s, pos0, rng)
            arm("mut_random", f"{label} {old}{pos0 + 1}{new} (random)", mut,
                "a substitution anywhere in the protein",
                {"parent": label, "position": pos0 + 1,
                 "from_aa": old, "to_aa": new})

    print("\nEPITOPE-LEVEL SUBSTITUTIONS"
          "\n  The whole-protein response saturates, so a substitution test "
          "at that level\n  measures noise. On the peptide itself the "
          "question is sharp.")
    ep_level = epitope_level_anchor_test(models, host_ix, seqs, rng,
                                         n=6 if quick else 30)
    if ep_level:
        print(f"  anchor substitution     : affinity "
              f"{ep_level['anchor_loss']:.1f}x worse "
              f"(n={ep_level['n_anchor']})")
        print(f"  non-anchor substitution : affinity "
              f"{ep_level['non_anchor_loss']:.1f}x worse "
              f"(n={ep_level['n_non_anchor']})")

    tab = pd.DataFrame(rows)
    tab.to_csv(RESULTS / f"{out}.csv", index=False)
    summary = verdict(tab)
    summary["epitope_level"] = {k: v for k, v in ep_level.items()
                                if k != "table"}
    summary["checks"] = []
    summary = verdict(tab, ep_level)
    (RESULTS / f"{out}.json").write_text(_json.dumps(summary, indent=1))
    report(summary)
    print(f"\n  wrote {RESULTS / out}.[csv|json]")
    return tab, summary


# ---------------------------------------------------------------------------
# what the arms mean
# ---------------------------------------------------------------------------
def _med(tab, arm, parent=None):
    s = tab[tab.arm == arm]
    if parent is not None and "parent" in s.columns:
        s = s[s.parent == parent]
    return float(s.peak_IgG.median()) if len(s) else float("nan")


def verdict(tab, ep_level=None):
    """Turn the arms into pass/fail statements with the numbers behind them."""
    out = {"checks": [], "n_arms": int(len(tab)),
           "epitope_level": {k: v for k, v in (ep_level or {}).items()
                             if k != "table"}}

    def add(name, ok, detail, why):
        out["checks"].append({"name": name, "ok": bool(ok),
                              "detail": detail, "why": why})

    pos = tab[tab.arm == "positive"]
    neg = tab[tab.arm == "negative"]
    pos_med = float(pos.peak_IgG.median()) if len(pos) else float("nan")
    neg_max = float(neg.peak_IgG.max()) if len(neg) else float("nan")

    add("positive controls respond", len(pos) and pos.peak_IgG.min() > 0,
        f"median peak IgG {pos_med:,.0f} across {len(pos)} antigens",
        "a vaccine antigen the model cannot respond to would be useless")

    add("negative controls do not", len(neg) and neg_max <= 0,
        f"every one of the {len(neg)} negative controls gives exactly 0"
        if neg_max <= 0 else
        f"the strongest negative control still reaches {neg_max:,.0f}",
        "self proteins and sequences with no epitope structure must give "
        "nothing, each for its own reason")

    add("positive and negative are separated",
        neg_max <= 0 or pos_med > 10 * max(neg_max, 1e-9),
        f"positives {pos_med:,.0f} against the worst negative {neg_max:,.0f}",
        "this is the comparison C-ImmSim failed: an antigen and a "
        "housekeeping protein must not come out the same")

    # Scrambles. The first version of this check asserted that the real
    # protein out-responds its shuffled replicates, and that was the wrong
    # test: a shuffled protein is still a foreign protein of the same size,
    # and it SHOULD be immunogenic. What must change is WHICH epitopes are
    # found. So the check is on epitope identity, and the magnitudes are
    # reported beside it without being asserted on.
    scr_rows = []
    for parent in tab[tab.arm == "scramble"].get("parent", pd.Series(dtype=object)).dropna().unique():
        real = float(tab[(tab.arm == "positive") & (tab.label == parent)].peak_IgG.iloc[0])
        sc = tab[(tab.arm == "scramble") & (tab.parent == parent)]
        ov = sc.get("epitope_overlap", pd.Series(dtype=float)).dropna()
        scr_rows.append({"parent": parent, "real": real,
                         "scramble_median": float(sc.peak_IgG.median()),
                         "scramble_max": float(sc.peak_IgG.max()),
                         "n": int(len(sc)),
                         "epitope_overlap_median":
                             round(float(ov.median()), 4) if len(ov) else None})
    out["scrambles"] = scr_rows
    ovs = [r["epitope_overlap_median"] for r in scr_rows
           if r["epitope_overlap_median"] is not None]
    add("the epitopes come from the sequence, not the composition",
        bool(ovs) and max(ovs) < 0.02,
        "; ".join(f"{r['parent']}: {r['epitope_overlap_median']:.1%} of "
                  f"epitopes survive shuffling; response {r['real']:,.0f} vs "
                  f"{r['scramble_median']:,.0f}" for r in scr_rows) or "not run",
        "a shuffled protein has identical composition, so if the same "
        "peptides came back the model would be counting residues. It is "
        "EXPECTED that a shuffled protein is still immunogenic - it is still "
        "foreign - so the magnitudes are reported, not asserted on")

    # mutations: anchor should cost more than non-anchor
    mut = {}
    for parent in tab[tab.arm.astype(str).str.startswith("mut_")].get(
            "parent", pd.Series(dtype=object)).dropna().unique():
        real = float(tab[(tab.arm == "positive") & (tab.label == parent)].peak_IgG.iloc[0])
        row = {"parent": parent, "real": real}
        for k in ("mut_anchor", "mut_non_anchor", "mut_random"):
            m = _med(tab, k, parent)
            row[k] = m
            row[k + "_drop"] = (round(100 * (real - m) / real, 2)
                                if real > 0 and m == m else None)
        mut[parent] = row
    out["mutations"] = list(mut.values())

    ep = out.get("epitope_level") or {}
    if ep:
        add("an anchor substitution destroys binding, a non-anchor one does not",
            ep.get("anchor_loss", 0) > 3 * max(ep.get("non_anchor_loss", 0), 1e-9),
            f"predicted affinity worsens {ep['anchor_loss']:.1f}-fold when an "
            f"anchor residue is changed, against {ep['non_anchor_loss']:.1f}-fold "
            f"for a residue pointing at the T cell "
            f"(n={ep['n_anchor']} and {ep['n_non_anchor']} substitutions)",
            "measured on the peptide itself, where it cannot be hidden by a "
            "whole-protein response that saturates. This is the check that "
            "the model has learned where the MHC groove grips")

    drops_a = [r["mut_anchor_drop"] for r in mut.values()
               if r.get("mut_anchor_drop") is not None]
    drops_n = [r["mut_non_anchor_drop"] for r in mut.values()
               if r.get("mut_non_anchor_drop") is not None]
    # NOT a pass/fail. The first version of this asserted that a single
    # substitution moves the whole-protein response, and it "passed" on
    # drops of 0.09% and -0.49% - noise, including a negative one. A check
    # that passes on noise is worse than no check. The honest statement is
    # that at the whole-protein level one substitution is NOT detectable,
    # and the epitope-level check above is where the question is answerable.
    out["limitations"] = [{
        "name": "a single substitution does not move the whole-protein response",
        "detail": "; ".join(
            f"{r['parent']}: anchor {r['mut_anchor_drop']}%, non-anchor "
            f"{r['mut_non_anchor_drop']}%, random {r['mut_random_drop']}%"
            for r in mut.values()) or "not run",
        "why": "an antigen with thousands of predicted binders does not "
               "notice losing one, so escape has to be studied at the "
               "epitope, not by re-simulating the mutated protein and "
               "comparing totals",
    }]
    out["whole_protein_mutation_effect"] = {
        "anchor_mean_drop_pct": round(float(np.mean(drops_a)), 3) if drops_a else None,
        "non_anchor_mean_drop_pct": round(float(np.mean(drops_n)), 3) if drops_n else None,
    }

    out["passed"] = sum(1 for c in out["checks"] if c["ok"])
    out["total"] = len(out["checks"])
    return out


def report(summary):
    print("\n" + "=" * 74)
    print("VERDICT")
    print("=" * 74)
    for c in summary["checks"]:
        print(f"  {'PASS' if c['ok'] else 'FAIL'}  {c['name']}")
        print(f"        {c['detail']}")
    for l in summary.get("limitations", []):
        print(f"  NOTE  {l['name']}")
        print(f"        {l['detail']}")
        print(f"        {l['why']}")
    print(f"\n  {summary['passed']} of {summary['total']} checks passed")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=8,
                    help="replicates per randomised arm")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="sensitivity")
    ap.add_argument("--antigen", default=None,
                    help="your own protein: a UniProt accession, a FASTA "
                         "file, or a bare sequence. It is put through the "
                         "same controls as the built-in positives.")
    a = ap.parse_args()
    run(reps=a.reps, seed=a.seed, quick=a.quick, out=a.out,
        antigen=a.antigen)
