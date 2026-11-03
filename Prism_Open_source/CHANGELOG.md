# Changelog

## 1.0.0

First stable release.

### Fluid interaction (no blocking pop-ups)
* Deleting a table, result, graph or layout happens at once; a small notice offers **Undo** (8 s) instead of an
  "Are you sure?" dialog.
* Renaming and other quick text inputs open as a small field next to where you clicked (Enter accepts, Esc cancels,
  clicking elsewhere accepts) instead of a dimmed full-screen dialog.
* When the server is momentarily busy, calculations wait and retry quietly instead of showing an error; the default
  per-user allowance was raised from 120 to 600 calculations per minute and from 2 to 4 at once.
* Confirmations remain only for actions that cannot be undone (deleting a project from the server, discarding unsaved
  work, resetting encryption keys).

### Robustness (found by systematic testing)
* Every analysis was run on every table type in 11 data variants (empty, one row, two rows, identical values, missing
  values, 1e±150 magnitudes, negatives, zeros, an empty column) — 2,453 runs — and with every dialog option.
* Fixed: results containing "not a number" (e.g. a t test on identical values) were sent as invalid JSON, so the
  browser could not show them at all. All server responses are now valid JSON (not-a-number values are shown blank).
* Fixed: the frequency distribution failed whenever *Bin width* or *First bin* was left blank (the default). Blank
  optional fields now mean "automatic" in every analysis.
* Degenerate data (all values identical, too few rows, a variable that is a combination of others) now give a clear
  explanation instead of "the analysis failed"; the mixed-model message explains when rows do not behave like
  matched subjects.
* A second sweep fed every analysis invalid option values (empty, null, text, negative, huge, wrong type; 1,903
  cases). Fixed: an empty option could crash the server; an unknown option value (e.g. from an old project file) gave
  "the analysis failed"; a negative column number silently analysed the *last* column instead of being refused
  (column choices are now checked centrally); absurd permutation or lag counts could run for minutes (now capped).
* A polynomial order too high for the number of different X values is explained instead of failing.
* Opening the app with `http://` on the HTTPS port now redirects instead of silently dropping the connection; a
  port that is already in use (macOS AirPlay uses 5000) is reported with the fix.

## 0.9.0
* Gray's test for cumulative incidence curves and Fine-Gray competing-risks regression (validated against cmprsk).

## 0.8.0
* Strict number format, Satterthwaite mixed models, summary-data tables, encrypted sharing, images on layouts,
  calculated columns, more power analyses, competing-risks cumulative incidence, continuous integration.

## 0.7.0
* Random-slope mixed models, GraphPad Prism .pzfx import/export, EPS and EMF export, shapes and templates, activity
  log, client-side graph summaries and faster start-up.
