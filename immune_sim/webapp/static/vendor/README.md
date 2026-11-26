# Vendored, so the GUI works with no internet

`plotly-basic.min.js` — Plotly.js 2.35.2, basic bundle (1.0 MB), MIT licence,
from `cdn.jsdelivr.net/npm/plotly.js-basic-dist-min@2.35.2`.

It is committed rather than loaded from a CDN because this project is meant to
run offline: the datasets ship in the repository for the same reason. Loading
charts from a CDN meant that on a machine without internet every plot was a
silent empty rectangle — which is how this was found.

The *basic* bundle is used because the only trace type here is `scatter`. The
full bundle is roughly four times the size and adds nothing this GUI draws.

To update: download the same path at a new version, replace the file, and
change the `<script src>` in `templates/base.html` if the filename changes.
