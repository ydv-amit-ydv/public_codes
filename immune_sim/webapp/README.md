# immune_sim GUI

```bash
pip install -r ../requirements.txt
python3 app.py                 # http://127.0.0.1:5000
python3 app.py --port 8080 --host 0.0.0.0
```

Four pages, each calling the same functions as the CLI so the GUI cannot
drift away from the validated pipeline.

| Page | What it does |
|---|---|
| **Simulate** | Platform, adjuvant, dose, schedule and epitope set, then the full time course: antigen supply, antibody, T cells, B cells and plasma cells, innate activation, immunodominance per clone. Downloads the time course as CSV. |
| **Compare** | Up to six formulations on one axis, with a summary table. Carries the warning that platform rankings are a structured prior, not a validated prediction. |
| **Epitopes** | Fetch a UniProt accession or paste a sequence, scan it with the trained models, and write the result where the Simulate page can pick it up. |
| **Validation** | The evidence behind every claim, read live from `../results/` — the scorecard, leave-one-pathogen-out binding results, and the held-out CMI-PB comparison. |

## Notes

- **Charts** are Plotly, themed from the CSS variables, so they follow the
  browser's light/dark setting rather than sitting on a white rectangle in a
  dark page. Plotly is loaded from a CDN; everything else is local.
- **The parameter set in use is printed in the header.** A figure cannot
  outlive the parameters that produced it, because nothing is cached: every
  request re-integrates.
- **Epitope files are read and written only inside `results/`.** The path is
  reduced to its basename before use, so a crafted `epitope_source` cannot
  reach the rest of the filesystem.
- `python3 app.py` runs Flask's development server, which is single-process
  and not hardened. It binds to localhost by default. For anything beyond one
  person on one machine, put it behind gunicorn or uwsgi.
