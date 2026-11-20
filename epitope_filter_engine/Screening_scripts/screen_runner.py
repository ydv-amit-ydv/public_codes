"""Safety and function screens for peptides and whole constructs (Raghava lab tools).

These are not MHC predictors and they do not vote in the binding consensus. They answer the
questions that decide whether a designed construct is usable at all:

    toxicity        is this peptide or construct toxic?         ToxinPred3
    allergenicity   will it be seen as an allergen?             AlgPred2
    IFN-gamma       will this HTL epitope induce IFN-gamma?     IFNepitope2
    haemolysis      will it lyse red blood cells? (optional)    HemoPI2

They all share one CLI shape (-i input, -o output, -t threshold, -m model, -d display) and
all write a small CSV, so one runner drives them all.

Two things this runner has to handle, both found by running them:

* **ToxinPred3 writes a file called `Sequence_1` into the current directory** and reads it
  back, so two runs in one directory would race. Every call gets its own scratch directory.
* They are pinned to **scikit-learn 1.2.x**: their models are pickled trees from before
  scikit-learn 1.3 added `missing_go_to_left` to the node dtype, and a newer scikit-learn
  refuses to unpickle them. The environment the installer builds pins that; `health()`
  reports it plainly rather than letting a run die with a dtype traceback.
"""
import csv
import os
import subprocess
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))

# key -> console script, extra flags, the column holding the score, and what a positive call
# looks like in the Prediction column
TOOLS = {
    "toxinpred3": {"cmd": "toxinpred3", "label": "ToxinPred3 (toxicity)",
                   "flags": ["-d", "2"], "score": ("Hybrid Score", "ML Score", "ML_Score"), "positive": "toxin",
                   "column": "Toxicity", "default_threshold": 0.38, "adverse": True,
                   "cite": "Rathore AS et al. (2024) Brief Bioinform 25(5):bbae208"},
    "algpred2": {"cmd": "algpred2", "label": "AlgPred2 (allergenicity)",
                 "flags": ["-d", "2"], "score": ("Hybrid Score", "ML_Score", "ML Score"), "positive": "allergen",
                 "column": "Allergenicity", "default_threshold": 0.3, "adverse": True,
                 "cite": "Sharma N et al. (2021) Brief Bioinform 22(4):bbaa294"},
    "ifnepitope2": {"cmd": "ifnepitope2", "label": "IFNepitope2 (IFN-gamma induction)",
                    "flags": ["-d", "2", "-j", "1"], "score": ("Total_Score", "ML_Score"),
                    "positive": "inducer", "column": "IFN_gamma", "default_threshold": 0.5, "adverse": False,
                    "cite": "Dhall A et al. (2023) Brief Bioinform 24(6):bbad354"},
    "ntxpred2": {"cmd": "ntxpred2", "label": "NTxPred2 (neurotoxicity)",
                 "flags": ["-d", "2"], "score": ("ML_Score", "Score", "ML Score"), "positive": "neurotoxin",
                 "column": "Neurotoxicity", "default_threshold": 0.5, "optional": True, "adverse": True,
                 "needs": "transformers and torch (it is an ESM model)",
                 "cite": "Rathore AS et al. (2024) NTxPred2, Raghava lab"},
    "hemopi2": {"cmd": "hemopi2_classification", "label": "HemoPI2 (haemolytic activity)",
                "flags": ["-d", "2"], "score": ("ML_Score", "Score", "ML Score"), "positive": "hemolytic",
                "column": "Haemolysis", "default_threshold": 0.5, "optional": True, "adverse": True,
                "cite": "Rathore AS et al. (2025) Nucleic Acids Res, HemoPI2.0"},
}

# padding for the batch-of-one bug above; any valid peptide does, it is never reported
FILLER = "AAAAAAAAAA"
FILLER2 = "CCCCCCCCCC"

HINT = ("Install the screening tools: bash start.sh --redo-setup, or press "
        "'Repair / install predictors' in Settings. They need scikit-learn 1.2.x, which the "
        "installer pins, because their models cannot be unpickled by newer versions.")


class ScreenError(RuntimeError):
    pass


def _python(cfg):
    """Interpreter of the screening environment."""
    return (cfg or {}).get("python") or ""


def available(cfg):
    """Which tools can actually run, as {key: (ok, detail)}."""
    py = _python(cfg)
    out = {}
    for key, spec in TOOLS.items():
        if not py:
            out[key] = (False, "no screening environment found")
            continue
        exe = os.path.join(os.path.dirname(py), spec["cmd"])
        out[key] = (os.path.exists(exe), exe if os.path.exists(exe) else f"{spec['cmd']} is not installed")
    return out


def _fasta(peptides):
    lines = []
    for i, p in enumerate(peptides, 1):
        lines.append(f">s{i}")
        lines.append(str(p).strip().upper())
    return "\n".join(lines) + "\n"


def _read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def predict(tool, peptides, cfg, threshold=None, timeout=1800):
    """Run one screen over `peptides`. Returns DataFrame[Peptide, <Column>, <Column>_Score].

    Peptides are sent as s1..sN and mapped back by that id, not by the sequence the tool
    echoes: IFNepitope2 can return a window rather than the whole peptide.
    """
    spec = TOOLS.get(tool)
    if spec is None:
        raise ScreenError(f"unknown screen: {tool}")
    peptides = [str(p).strip().upper() for p in peptides if str(p).strip()]
    col = spec["column"]
    if not peptides:
        return pd.DataFrame(columns=["Peptide", col, f"{col}_Score"])

    py = _python(cfg)
    if not py:
        raise ScreenError(f"{spec['label']} is not installed. " + HINT)
    exe = os.path.join(os.path.dirname(py), spec["cmd"])
    if not os.path.exists(exe):
        raise ScreenError(f"{spec['label']} is not installed. " + HINT)

    # ToxinPred3 and AlgPred2 both fail on a batch of ONE: their feature matrix
    # collapses to a single row and numpy/scikit-learn reject it ("axis 1 is out of bounds",
    # "Reshape your data... if it contains a single sample"). Send a filler alongside and
    # drop it afterwards; the mapping is by id, so the filler is never looked up.
    sent = list(peptides)
    if len(sent) < 2:
        sent = sent + [FILLER if sent[0] != FILLER else FILLER2]

    # its own directory: ToxinPred3 drops a file named Sequence_1 in the working directory
    with tempfile.TemporaryDirectory(prefix=f"screen_{tool}_") as d:
        fin, fout = os.path.join(d, "in.fa"), os.path.join(d, "out.csv")
        with open(fin, "w") as f:
            f.write(_fasta(sent))
        cmd = [exe, "-i", fin, "-o", fout, *spec["flags"]]
        if threshold is not None:
            cmd += ["-t", str(threshold)]
        env = dict(os.environ, PYTHONUNBUFFERED="1")
        try:
            r = subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                               stdin=subprocess.DEVNULL, timeout=timeout, env=env)
        except subprocess.TimeoutExpired:
            raise ScreenError(f"{spec['label']} took longer than {timeout}s")
        if not os.path.exists(fout):
            tail = (r.stderr or r.stdout or "").strip().splitlines()[-3:]
            raise ScreenError(f"{spec['label']} produced no result: " + " / ".join(tail))
        rows = _read_csv(fout)
    if not rows:
        raise ScreenError(f"{spec['label']} returned an empty table for {len(peptides)} peptide(s)")

    by_id = {f"s{i}": p for i, p in enumerate(sent, 1)}
    score_key = next((k for k in spec["score"] if rows and k in rows[0]), None)
    # "Subject" is what ToxinPred3 calls it under the hybrid model, where there is also no
    # Sequence column, so the id is the only way back to the peptide
    id_key = next((k for k in ("ID", "Seq_ID", "SeqID", "Id", "Subject") if rows and k in rows[0]), None)
    pred_key = next((k for k in ("Prediction", "Predicted", "Result") if rows and k in rows[0]), None)

    out = {}
    for row in rows:
        pid = (row.get(id_key) or "").strip()
        pep = by_id.get(pid) or (row.get("Sequence") or "").strip().upper()
        if not pep or pep not in set(peptides):          # the padding row, if one was added
            continue
        try:
            score = float(row.get(score_key)) if score_key and row.get(score_key) not in (None, "") else None
        except ValueError:
            score = None
        call = (row.get(pred_key) or "").strip()
        # a tool that reports several windows for one peptide: keep the worst (highest score)
        prev = out.get(pep)
        if prev is None or (score is not None and prev[1] is not None and score > prev[1]):
            out[pep] = (call, score)

    data = [{"Peptide": p, col: out.get(p, ("", None))[0], f"{col}_Score": out.get(p, ("", None))[1]}
            for p in peptides]
    return pd.DataFrame(data, columns=["Peptide", col, f"{col}_Score"])


# ---------------------------------------------------------------- allergen homology (no ML)
#
# A second, independent allergenicity check: sequence similarity to known allergens, which is
# the regulatory criterion rather than a model's opinion. FAO/WHO Codex alimentarius asks for
#
#     >= 35 % identity over a window of 80 amino acids, or
#     >= 6 contiguous identical amino acids
#
# with a known allergen. The first limb only means anything for construct-length sequences;
# the second is the one that bites for epitopes.
#
# It needs no download: AlgPred2 ships both a blastp binary and its allergen database. That
# database also holds AlgPred2's NON-allergen training sequences, named N_*, so only P_*
# subjects may be counted - treating an N_ hit as allergen homology would invert the answer.
#
# Known limitation, and it is the rule's rather than this code's: the 6-contiguous limb
# over-calls on low-complexity sequence. A poly-aspartate stretch matches an acidic run in
# some allergen at 100 % and is flagged, which is why this is reported next to AlgPred2's
# opinion rather than instead of it. Two methods agreeing means something; either alone is
# a prompt to look, not a verdict.

def _algpred_dir(cfg):
    py = _python(cfg)
    if not py:
        return None
    import glob
    root = os.path.dirname(os.path.dirname(py))
    for pat in (os.path.join(root, "lib", "python*", "site-packages", "algpred2"),
                os.path.join(root, "lib", "site-packages", "algpred2")):
        hit = glob.glob(pat)
        if hit:
            return hit[0]
    return None


def _blastp(cfg):
    d = _algpred_dir(cfg)
    if not d:
        return None, None
    import platform
    plat = {"Linux": "linux", "Darwin": "mac", "Windows": "windows"}.get(platform.system(), "linux")
    exe = os.path.join(d, "blast_binaries", plat, "blastp" + (".exe" if plat == "windows" else ""))
    db = os.path.join(d, "blast_db", "data")
    return (exe if os.path.exists(exe) else None), (db if os.path.exists(db + ".phr") else None)


def allergen_homology(peptides, cfg, min_identity=35.0, min_window=80, contiguous=6, evalue=10.0):
    """Best allergen match per peptide, and whether it meets the FAO/WHO criterion.

    Returns DataFrame[Peptide, Allergen_Homology, Allergen_Identity, Allergen_Match,
    Allergen_Align_Len, Allergen_Evalue].
    """
    cols = ["Peptide", "Allergen_Homology", "Allergen_Identity", "Allergen_Match",
            "Allergen_Align_Len", "Allergen_Evalue"]
    peptides = [str(p).strip().upper() for p in peptides if str(p).strip()]
    if not peptides:
        return pd.DataFrame(columns=cols)
    exe, db = _blastp(cfg)
    if not exe or not db:
        raise ScreenError("allergen homology needs AlgPred2's blastp and database. " + HINT)

    with tempfile.TemporaryDirectory(prefix="screen_homology_") as d:
        fin = os.path.join(d, "in.fa")
        with open(fin, "w") as f:
            f.write(_fasta(peptides))
        # short queries need the short-word-size task, or BLAST finds nothing at all
        task = "blastp-short" if max(len(p) for p in peptides) < 30 else "blastp"
        cmd = [exe, "-query", fin, "-db", db, "-task", task, "-evalue", str(evalue),
               "-max_target_seqs", "50", "-outfmt", "6 qseqid sseqid pident length evalue"]
        try:
            r = subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                               stdin=subprocess.DEVNULL, timeout=1800)
        except subprocess.TimeoutExpired:
            raise ScreenError("the allergen homology search timed out")
        if r.returncode != 0 and not r.stdout:
            raise ScreenError("blastp failed: " + (r.stderr or "").strip().splitlines()[-1:][0]
                              if (r.stderr or "").strip() else "blastp failed")

    by_id = {f"s{i}": p for i, p in enumerate(peptides, 1)}
    best = {}
    for line in r.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) < 5:
            continue
        qid, sid, pident, length, ev = parts[0], parts[1], parts[2], parts[3], parts[4]
        if not sid.startswith("P_"):          # N_* are AlgPred2's non-allergens
            continue
        pep = by_id.get(qid)
        if not pep:
            continue
        try:
            pident, length, ev = float(pident), int(length), float(ev)
        except ValueError:
            continue
        hit = (pident >= min_identity and length >= min_window) or (pident >= 99.9 and length >= contiguous)
        prev = best.get(pep)
        rank = (hit, pident * length)         # a qualifying hit first, then the strongest
        if prev is None or rank > prev[0]:
            best[pep] = (rank, {"Allergen_Homology": bool(hit), "Allergen_Identity": round(pident, 2),
                                "Allergen_Match": sid, "Allergen_Align_Len": length,
                                "Allergen_Evalue": ev})
    rows = []
    for p in peptides:
        v = best.get(p)
        rows.append({"Peptide": p, **(v[1] if v else {"Allergen_Homology": False, "Allergen_Identity": None,
                                                      "Allergen_Match": "", "Allergen_Align_Len": None,
                                                      "Allergen_Evalue": None})})
    return pd.DataFrame(rows, columns=cols)


def screen(peptides, cfg, tools=None, thresholds=None, homology=True):
    """Run several screens and return one table, plus per-tool notes.

    A tool that is missing or fails does not stop the others: its columns are left out and the
    reason is returned, the same way a failed model does not end a prediction run.
    """
    peptides = [str(p).strip().upper() for p in peptides if str(p).strip()]
    tools = [t for t in (tools or ["toxinpred3", "algpred2", "ifnepitope2"]) if t in TOOLS]
    table = pd.DataFrame({"Peptide": peptides})
    notes = []
    for t in tools:
        try:
            df = predict(t, peptides, cfg, (thresholds or {}).get(t))
            table = table.merge(df, on="Peptide", how="left")
        except ScreenError as e:
            msg = str(e)
            notes.append(msg if TOOLS[t]["label"] in msg else f"{TOOLS[t]['label']}: {msg}")
        except Exception as e:                            # noqa: BLE001 - one screen must not sink the rest
            notes.append(f"{TOOLS[t]['label']} failed: {type(e).__name__}: {e}")
    if homology:
        try:
            table = table.merge(allergen_homology(peptides, cfg), on="Peptide", how="left")
        except ScreenError as e:
            notes.append(str(e))
        except Exception as e:                            # noqa: BLE001
            notes.append(f"allergen homology failed: {type(e).__name__}: {e}")
    return table, notes


def is_positive(call, token):
    """Did the tool say yes?

    The verdicts are worded "Toxin"/"Non-Toxin", "Allergen"/"Non-Allergen",
    "IFN-gamma inducer"/"Non-inducer", so a plain substring test matches BOTH: "Non-Allergen"
    contains "allergen". The negatives have to be excluded first.
    """
    v = str(call or "").strip().lower().replace("\u2013", "-").replace("\u2014", "-")
    if not v:
        return False
    if v.startswith(("non-", "non ", "nontoxic", "nonallergen", "noninducer", "nonhemolytic", "no")):
        return False
    return token in v


def flags(table, adverse_only=True):
    """Peptides a screen called a problem: {peptide: [reasons]}.

    Only the adverse screens count. A positive IFNepitope2 call is a reason to KEEP an HTL
    epitope, not to drop it, so treating "positive" as "bad" across the board would throw
    away the good ones.
    """
    bad = {}
    for _, row in table.iterrows():
        why = [spec["column"] for spec in TOOLS.values()
               if (spec["adverse"] or not adverse_only)
               and spec["column"] in table.columns
               and is_positive(row.get(spec["column"]), spec["positive"])]
        # the homology check is not an ML tool, but it is adverse evidence all the same
        if "Allergen_Homology" in table.columns and bool(row.get("Allergen_Homology")):
            why.append("Allergen_Homology")
        if why:
            bad[row["Peptide"]] = why
    return bad


def benefits(table):
    """Peptides a non-adverse screen called positive, e.g. the IFN-gamma inducers."""
    good = {}
    for _, row in table.iterrows():
        why = [spec["column"] for spec in TOOLS.values()
               if not spec["adverse"] and spec["column"] in table.columns
               and is_positive(row.get(spec["column"]), spec["positive"])]
        if why:
            good[row["Peptide"]] = why
    return good
