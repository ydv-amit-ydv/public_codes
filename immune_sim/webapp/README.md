# immune_sim web GUI

    pip install -r ../requirements.txt
    python3 app.py
    # http://127.0.0.1:5000

Everything this project can do is reachable from here. A fresh clone carries
no data and no trained weights, so the GUI's first job is getting them: the
**Data** page downloads each dataset from its source, the **Models** page
trains the predictors, and **Calibrate** refits the dynamic model. Nothing
needs a terminal.

## Pages

| Page | What it does | Needs |
|---|---|---|
| **Simulate** | Platform, adjuvant, dose and schedule → the full T/B/antibody course, with every compartment plotted and the frame downloadable as CSV. | nothing |
| **Compare** | Several formulations on one axis. | nothing |
| **Epitopes** | Fetch a UniProt accession or paste a sequence → ranked epitopes → feed them straight to a simulation. | trained models |
| **Data** | Download IEDB affinities, IEDB T cell assays, CMI-PB human boost data, mass-spec eluted ligands, and ClinicalTrials.gov fold-rises. | network |
| **Models** | Train the per-allele models, the pan-allele core-scan ensemble, and the immunogenicity model; shows each one's held-out metrics. | data |
| **Calibrate** | Refit the 36 ODE parameters against measured human data, with a live score. | CMI-PB data |
| **Validation** | Every claim graded A–D, the measured evidence behind each, and a button to re-run any of it. | varies |

The Simulate, Compare and Validation pages work in a clean clone, because the
calibrated parameters and a demo epitope set are committed. Only the epitope
predictions need the harvest-and-train path.

## How the long operations work

Each button starts the project's **own CLI script** as a subprocess and
streams its output into the page: live log, a progress bar where the script
reports progress, elapsed time, and a Cancel that kills the whole process
group. The GUI does not reimplement any analysis. That is deliberate — a
figure in the browser and a number in `VALIDATION.md` come from the same code
path, so the GUI cannot drift from the validated pipeline, and anything you
see here is reproducible from a terminal.

The browser never sends a command. It sends a **task name**; `tasks.py` turns
that into an argument vector and clamps every parameter. That allowlist is the
only reason this is safe to bind to a port.

## Two safeguards worth knowing about

**Calibration never overwrites the parameter set in use.** It writes a
timestamped `models/calibrated_candidate_*.json`, which nothing reads until
you promote it. This exists because the first version of this GUI did
overwrite it, and a twenty-second test run replaced a converged fit
(score 2.5e-05) with a search that had reached 9.9 — the search writes a
checkpoint on every improvement, so even an immediately-cancelled run has
already clobbered the file. Promotion is a separate click, it keeps the
previous set as `calibrated_wetlab_prev.json`, and it does **not** refuse a
candidate with a worse training score: a lower training score is not a better
model, and the held-out split is what decides.

**Results are marked stale, not silently served.** Every output knows what it
was computed from. Retrain a model and the validation numbers computed from
the old one are flagged `stale — models/pan_report.json is newer`, in the
card and in the "what is on disk" table. A validation number that looks
current but was computed from a model that has since been retrained is worse
than no number at all.

## Caveats

This is Flask's development server bound to `127.0.0.1`. It is a single-user
local tool, not a deployment: there is no authentication, and a job started
from one browser tab is visible to any other. Jobs live in memory, so
restarting the server loses the log of a running job — the subprocess itself
keeps going, and its output file is still written.
