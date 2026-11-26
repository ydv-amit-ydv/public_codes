/* immune_sim GUI front end.
   Plotly is themed from the CSS variables so the charts follow light/dark
   instead of sitting on a white rectangle in a dark page. */
"use strict";

const css = (n) => getComputedStyle(document.documentElement)
  .getPropertyValue(n).trim();

const PALETTE = ["#2b6cb0", "#1a7f4b", "#b4530a", "#7a4fbf",
                 "#0d8b8b", "#b03060", "#6b7280", "#8a6d00"];

function layout(title, yTitle, opts) {
  const ink = css("--fg"), muted = css("--muted"), grid = css("--grid");
  return Object.assign({
    title: { text: title, font: { size: 13.5, color: ink }, x: 0, xanchor: "left" },
    paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
    font: { color: muted, size: 11 },
    margin: { l: 58, r: 14, t: 34, b: 42 },
    xaxis: { title: "day", gridcolor: grid, zerolinecolor: grid,
             linecolor: grid },
    yaxis: { title: yTitle, gridcolor: grid, zerolinecolor: grid,
             linecolor: grid },
    legend: { orientation: "h", y: -0.22, font: { size: 10 } },
    showlegend: true, hovermode: "x unified",
  }, opts || {});
}
const CONFIG = { displaylogo: false, responsive: true,
                 modeBarButtonsToRemove: ["lasso2d", "select2d"] };

function draw(el, traces, lay) {
  // If the chart library is missing, say so in the box rather than leaving an
  // empty rectangle. The numbers above the plots are the real content; a
  // missing chart should not look like a missing result.
  if (typeof Plotly === "undefined") {
    el.innerHTML = '<p class="muted sm" style="padding:12px">The chart library '
      + 'did not load, so this plot is missing. The numbers above are '
      + 'unaffected, and the CSV download has the full time course. '
      + 'Check that <code>webapp/static/vendor/plotly-basic.min.js</code> '
      + 'is present.</p>';
    return;
  }
  traces.forEach((t, i) => {
    t.type = t.type || "scatter"; t.mode = t.mode || "lines";
    t.line = Object.assign({ width: 2, color: PALETTE[i % PALETTE.length] },
                           t.line || {});
  });
  Plotly.react(el, traces, lay, CONFIG);
}

function setStatus(msg, isErr) {
  const s = document.getElementById("status");
  if (!s) return;
  s.textContent = msg || "";
  s.className = "status" + (isErr ? " err" : "");
}

async function post(url, body) {
  const r = await fetch(url, {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const j = await r.json().catch(() => ({ error: `HTTP ${r.status}` }));
  if (!r.ok || j.error) throw new Error(j.error || `HTTP ${r.status}`);
  return j;
}

const fmt = (v) => {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  const a = Math.abs(v);
  if (a === 0) return "0";
  if (a >= 1e5 || a < 1e-3) return v.toExponential(2);
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  return v.toFixed(2);
};

function fillSelect(el, items, valueKey, labelKey, keep) {
  if (!el) return;
  const prev = keep ? el.value : null;
  for (const it of items) {
    const o = document.createElement("option");
    o.value = it[valueKey]; o.textContent = it[labelKey];
    el.appendChild(o);
  }
  if (prev) el.value = prev;
}

function epitopeOptions(el) {
  if (!el) return;
  el.innerHTML = "";
  const o = document.createElement("option");
  o.value = "demo"; o.textContent = "demo (SARS-CoV-2 spike)";
  el.appendChild(o);
  for (const f of (META.epitope_files || [])) {
    if (f === "demo_epitopes.csv") continue;
    const x = document.createElement("option");
    x.value = f; x.textContent = f.replace(/_epitopes\.csv$/, "");
    el.appendChild(x);
  }
}

function days(id) {
  return (document.getElementById(id).value || "0")
    .split(",").map((s) => parseFloat(s.trim()))
    .filter((v) => Number.isFinite(v) && v >= 0);
}

/* ------------------------------------------------------------------ simulate */
function initSimulate() {
  const $ = (id) => document.getElementById(id);
  fillSelect($("platform"), META.platforms, "key", "label");
  fillSelect($("adjuvant"), META.adjuvants, "key", "label");
  fillSelect($("heterologous"), META.platforms, "key", "label");
  epitopeOptions($("epitope_source"));
  $("platform").value = "protein";
  $("adjuvant").value = "alum";

  const notes = () => {
    const p = META.platforms.find((x) => x.key === $("platform").value);
    const a = META.adjuvants.find((x) => x.key === $("adjuvant").value);
    if (p) {
      $("platform-note").textContent =
        `${p.notes} — MHC-I routing ${p.mhc1}, ` +
        (p.template ? `template t½ ${p.half_life_h} h` : "bolus antigen");
      $("dose-range").textContent = `typical ${p.dose[0]}–${p.dose[1]}`;
    }
    if (a) {
      $("adjuvant-note").textContent =
        `${a.notes} — Th1 bias ${a.th1_bias} ` +
        `(${a.th1_bias > 0.55 ? "CTL-leaning" :
            a.th1_bias < 0.45 ? "antibody-leaning" : "balanced"})`;
    }
  };
  $("platform").addEventListener("change", notes);
  $("adjuvant").addEventListener("change", notes);
  notes();

  const payload = () => ({
    platform: $("platform").value,
    adjuvant: $("adjuvant").value,
    heterologous: $("heterologous").value || null,
    dose_ug: $("dose_ug").value || null,
    inject_days: days("inject_days"),
    days: parseFloat($("days").value) || 240,
    epitope_source: $("epitope_source").value,
    calibration: $("calibration").value,
  });

  $("run").addEventListener("click", async () => {
    const btn = $("run"); btn.disabled = true;
    setStatus("integrating…");
    try {
      const j = await post("/api/simulate", payload());
      renderRun(j);
      setStatus(`${j.info.n_clones} epitope clones · ` +
                `injections day ${j.info.injection_days.join(", ")}`);
    } catch (e) { setStatus(e.message, true); }
    btn.disabled = false;
  });

  $("dl").addEventListener("click", async () => {
    const r = await fetch("/api/timecourse.csv", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload()),
    });
    if (!r.ok) { setStatus("download failed", true); return; }
    const url = URL.createObjectURL(await r.blob());
    const a = document.createElement("a");
    a.href = url; a.download = "timecourse.csv"; a.click();
    URL.revokeObjectURL(url);
  });

  $("run").click();
}

const CARDS = [
  ["peak_IgG", "Peak IgG"], ["peak_IgG_day", "at day"],
  ["boost_fold_over_baseline", "Boost fold-rise"],
  ["boost_peak_day_after", "Boost peak (d)"],
  ["peak_CTL", "Peak CTL"], ["peak_CD4", "Peak CD4"],
  ["affinity_fold", "Affinity ×"], ["IgG_over_IgM", "IgG:IgM"],
  ["CTL_memory_end", "CTL memory"], ["B_memory_end", "B memory"],
  ["peak_B", "Peak B cells"], ["LLPC_end", "Long-lived PC"],
  ["igg_durability", "IgG retained"],
];

function renderRun(j) {
  const s = j.summary;
  document.getElementById("summary").innerHTML = CARDS
    .filter(([k]) => s[k] !== undefined)
    .map(([k, lab]) =>
      `<div class="card"><div class="k">${lab}</div>` +
      `<div class="v">${fmt(s[k])}</div></div>`).join("");

  const g = (keys, titleY) => keys.filter((k) => j.series[k])
    .map((k) => ({ x: j.day, y: j.series[k],
                   name: (META.series.find((m) => m.key === k) || {}).label || k }));

  const logy = { yaxis: { type: "log", gridcolor: css("--grid"),
                          linecolor: css("--grid") } };
  draw("c-antigen",
       g(["antigen", "antigen_depot", "template", "antigen_endogenous"]),
       layout("Antigen supply", "arbitrary units", logy));
  draw("c-antibody", g(["IgM", "IgG"]),
       layout("Antibody", "arbitrary units", logy));
  draw("c-tcell", g(["TH_active", "TH_memory", "TC_active", "TC_memory"]),
       layout("T cells", "cells", logy));
  draw("c-bcell", g(["B_active", "B_memory", "plasma",
                     "plasma_longlived"]),
       layout("B cells and plasma cells", "cells", logy));
  draw("c-innate", g(["danger", "DC_mature", "IL2"]),
       layout("Innate activation and IL-2", "arbitrary units", logy));

  const cl = Object.keys(j.clones || {});
  if (cl.length) {
    draw("c-clones",
         cl.map((k) => ({ x: j.day, y: j.clones[k],
                          name: k.replace("TH_clone_", "CD4 ")
                                 .replace("TC_clone_", "CD8 ") })),
         layout("Immunodominance by clone", "cells", logy));
  }
  document.getElementById("describe").textContent = j.describe || "";
}

/* ------------------------------------------------------------------- compare */
function armRow(i) {
  const d = document.createElement("div");
  d.className = "arm";
  d.innerHTML =
    `<div class="row"><select class="p"></select><select class="a"></select>` +
    `<button class="x" title="remove">&times;</button></div>` +
    `<div class="row"><input class="l" placeholder="label (optional)">` +
    `<input class="d" type="number" step="0.5" placeholder="dose µg"></div>`;
  fillSelect(d.querySelector(".p"), META.platforms, "key", "label");
  fillSelect(d.querySelector(".a"), META.adjuvants, "key", "label");
  const defaults = [["protein", "alum"], ["mrna_lnp", "intrinsic"],
                    ["protein", "as01"], ["adenovirus", "intrinsic"],
                    ["circrna_lnp", "intrinsic"], ["protein", "none"]];
  const [p, a] = defaults[i % defaults.length];
  d.querySelector(".p").value = p; d.querySelector(".a").value = a;
  d.querySelector(".x").addEventListener("click", () => d.remove());
  return d;
}

function initCompare() {
  const $ = (id) => document.getElementById(id);
  epitopeOptions($("epitope_source"));
  const box = $("arms");
  [0, 1, 2].forEach((i) => box.appendChild(armRow(i)));
  $("add-arm").addEventListener("click", () => {
    if (box.children.length >= 6) { setStatus("six arms maximum"); return; }
    box.appendChild(armRow(box.children.length));
  });

  $("run").addEventListener("click", async () => {
    const arms = [...box.children].map((d) => ({
      platform: d.querySelector(".p").value,
      adjuvant: d.querySelector(".a").value,
      dose_ug: d.querySelector(".d").value || null,
      label: d.querySelector(".l").value || null,
    }));
    if (!arms.length) { setStatus("add at least one arm", true); return; }
    $("run").disabled = true;
    setStatus(`integrating ${arms.length} arms…`);
    try {
      const j = await post("/api/compare", {
        arms, inject_days: days("inject_days"),
        days: parseFloat($("days").value) || 240,
        epitope_source: $("epitope_source").value,
      });
      renderCompare(j);
      setStatus("");
    } catch (e) { setStatus(e.message, true); }
    $("run").disabled = false;
  });
  $("run").click();
}

function renderCompare(j) {
  const ok = j.arms.filter((a) => !a.error);
  const logy = { yaxis: { type: "log", gridcolor: css("--grid"),
                          linecolor: css("--grid") } };
  draw("c-igg", ok.map((a) => ({ x: a.day, y: a.IgG, name: a.label })),
       layout("IgG", "arbitrary units", logy));
  draw("c-ctl", ok.map((a) => ({ x: a.day, y: a.TC_active, name: a.label })),
       layout("CD8 T effector (CTL)", "cells", logy));

  const cols = [["peak_IgG", "peak IgG"], ["peak_IgG_day", "IgG peak day"],
                ["boost_fold_over_baseline", "boost fold-rise"],
                ["peak_CTL", "peak CTL"],
                ["CTL_over_IgG", "CTL:IgG"],
                ["CTL_memory_end", "CTL memory"],
                ["B_memory_end", "B memory"], ["dose_ug", "dose µg"]];
  document.getElementById("table").innerHTML =
    `<div class="tablewrap"><table><thead><tr><th>arm</th>` +
    cols.map(([, l]) => `<th>${l}</th>`).join("") +
    `</tr></thead><tbody>` +
    j.arms.map((a) => a.error
      ? `<tr><td>${a.label}</td><td colspan="${cols.length}">${a.error}</td></tr>`
      : `<tr><td>${a.label}</td>` +
        cols.map(([k]) => `<td>${fmt(a.summary[k])}</td>`).join("") + `</tr>`)
      .join("") + `</tbody></table></div>`;
}

/* ------------------------------------------------------------------ epitopes */
function initEpitopes() {
  const $ = (id) => document.getElementById(id);

  $("fetch").addEventListener("click", async () => {
    const acc = $("acc").value.trim();
    if (!acc) return;
    setStatus(`fetching ${acc}…`);
    try {
      const r = await fetch(`/api/uniprot/${encodeURIComponent(acc)}`);
      const j = await r.json();
      if (j.error) throw new Error(j.error);
      $("sequence").value = j.sequence;
      if (!$("name").value) $("name").value = acc;
      $("seqinfo").textContent = `${j.name} — ${j.length} residues`;
      setStatus("");
    } catch (e) { setStatus(e.message, true); }
  });

  $("run").addEventListener("click", async () => {
    $("run").disabled = true;
    setStatus("scanning every register against each allele…");
    try {
      const j = await post("/api/predict", {
        sequence: $("sequence").value, name: $("name").value, limit: 150,
      });
      renderEpitopes(j);
      setStatus(`saved as ${j.saved_as} — now selectable on the ` +
                `Simulate page`);
    } catch (e) { setStatus(e.message, true); }
    $("run").disabled = false;
  });
}

function renderEpitopes(j) {
  document.getElementById("summary").innerHTML = [
    ["Binders", j.n_binders], ["Strong", j.n_strong],
    ["Class I", j.by_class.I || 0], ["Class II", j.by_class.II || 0],
  ].map(([k, v]) => `<div class="card"><div class="k">${k}</div>` +
                    `<div class="v">${v}</div></div>`).join("");

  const ci = j.columns.indexOf("mhc_class");
  const si = j.columns.indexOf("start");
  const ii = j.columns.indexOf("pred_IC50_nM");
  const pi = j.columns.indexOf("peptide");
  const ai = j.columns.indexOf("allele");
  const traces = ["I", "II"].map((c, n) => {
    const rows = j.rows.filter((r) => r[ci] === c);
    return {
      x: rows.map((r) => r[si]), y: rows.map((r) => r[ii]),
      text: rows.map((r) => `${r[pi]} · ${r[ai]}`),
      name: `class ${c}`, mode: "markers", type: "scatter",
      marker: { size: 7, color: PALETTE[n], opacity: 0.75 },
    };
  }).filter((t) => t.x.length);
  Plotly.react("c-map", traces,
    layout("Predicted binders along the sequence", "IC50 (nM, log)", {
      xaxis: { title: "residue position", gridcolor: css("--grid"),
               linecolor: css("--grid") },
      yaxis: { title: "IC50 nM", type: "log", autorange: "reversed",
               gridcolor: css("--grid"), linecolor: css("--grid") },
    }), CONFIG);

  document.getElementById("table").innerHTML =
    `<div class="tablewrap"><table><thead><tr>` +
    j.columns.map((c) => `<th>${c}</th>`).join("") +
    `</tr></thead><tbody>` +
    j.rows.map((r) => `<tr>${r.map((v) => `<td>${v}</td>`).join("")}</tr>`)
      .join("") + `</tbody></table></div>`;
}

/* ------------------------------------------------------------------ tips
   The small circled "?" next to a field, same idiom as the engine: hover or
   focus shows one short paragraph explaining what the setting does and, where
   it matters, what the measured consequence of getting it wrong is. */
(function tips() {
  const tip = document.getElementById("tip");
  if (!tip) return;
  const show = (el) => {
    tip.innerHTML = el.dataset.tip || "";
    tip.style.display = "block";
    const r = el.getBoundingClientRect();
    const w = tip.offsetWidth, h = tip.offsetHeight;
    let x = r.left, y = r.bottom + 7;
    if (x + w > innerWidth - 10) x = Math.max(10, innerWidth - w - 10);
    if (y + h > innerHeight - 10) y = Math.max(10, r.top - h - 7);
    tip.style.left = x + "px";
    tip.style.top = y + "px";
  };
  const hide = () => { tip.style.display = "none"; };
  document.addEventListener("mouseover", (e) => {
    const q = e.target.closest("span.q[data-tip]");
    if (q) show(q); else if (!e.target.closest("#tip")) hide();
  });
  document.addEventListener("focusin", (e) => {
    const q = e.target.closest("span.q[data-tip]");
    if (q) show(q);
  });
  document.addEventListener("focusout", hide);
  addEventListener("scroll", hide, true);
  // Scoped to span.q[data-tip]: a bare ".q" also matched .verdict.q on the
  // validation page and replaced four whole verdicts with a "?".
  document.querySelectorAll("span.q[data-tip]").forEach((q) => {
    q.tabIndex = 0;
    q.textContent = "?";
  });
})();

/* ---------------------------------------------------------------- design
   The merged pipeline: a sequence goes in, epitopes and a simulated response
   come out of one request. Before this page the two halves of the project
   were joined only by a CSV on disk. */
function initDesign() {
  const $ = (id) => document.getElementById(id);
  fillSelect($("platform"), META.platforms, "key", "label");
  fillSelect($("adjuvant"), META.adjuvants, "key", "label");
  const het = $("heterologous");
  het.innerHTML = '<option value="">— same platform —</option>';
  META.platforms.forEach((p) => {
    const o = document.createElement("option");
    o.value = p.key; o.textContent = p.label; het.appendChild(o);
  });
  $("platform").value = "mrna_lnp";
  $("adjuvant").value = "intrinsic";

  // allele chips, from whichever models are trained
  fetch("/api/alleles").then((r) => r.json()).then((d) => {
    const box = $("alleles");
    if (!d.alleles || !d.alleles.length) {
      box.innerHTML = '<span class="muted">No trained models yet. '
        + '<a href="/setup">Set them up</a> — no network needed.</span>';
      return;
    }
    box.innerHTML = d.alleles.map((a) =>
      `<label class="chip"><input type="checkbox" value="${a.allele}"> `
      + `${a.allele} <span class="muted">${a.mhc_class}</span></label>`).join("");
  }).catch(() => { $("alleles").textContent = "could not list alleles"; });

  document.querySelectorAll('input[name="asrc"]').forEach((r) =>
    r.addEventListener("change", () => {
      $("box-sequence").hidden = r.value !== "sequence" || !r.checked
        ? document.querySelector('input[name="asrc"]:checked').value !== "sequence"
        : false;
      $("box-uniprot").hidden =
        document.querySelector('input[name="asrc"]:checked').value !== "uniprot";
    }));

  $("sequence").addEventListener("input", () => {
    const raw = $("sequence").value;
    const seq = raw.split("\n").filter((l) => !l.startsWith(">")).join("")
      .replace(/[^A-Za-z]/g, "").toUpperCase();
    const bad = [...new Set(seq.split(""))]
      .filter((c) => !"ACDEFGHIKLMNPQRSTVWY".includes(c));
    $("seqinfo").innerHTML = !seq.length ? "&nbsp;"
      : `${seq.length.toLocaleString()} residues`
        + (bad.length ? ` <span class="bad">— not amino acids: ${bad.join(", ")}`
           + (/^[ACGTU]+$/.test(seq) ? " (this looks like a DNA or RNA sequence)" : "")
           + "</span>" : "");
  });

  $("run").addEventListener("click", runDesign);
  $("dlcsv").addEventListener("click", () => downloadCsv(designPayload()));
}

function designPayload() {
  const $ = (id) => document.getElementById(id);
  const src = document.querySelector('input[name="asrc"]:checked').value;
  const alleles = [...document.querySelectorAll("#alleles input:checked")]
    .map((c) => c.value);
  return {
    antigen_source: src,
    sequence: $("sequence").value,
    accession: $("accession").value,
    name: $("name").value,
    alleles: alleles.length ? alleles : null,
    platform: $("platform").value,
    adjuvant: $("adjuvant").value,
    dose_ug: Number($("dose").value),
    heterologous: $("heterologous").value || null,
    inject_days: $("inject").value.split(",").map((s) => Number(s.trim()))
      .filter((n) => !Number.isNaN(n)),
    days: Number($("days").value),
    resolution_h: Number($("resolution").value),
    calibration: $("calibration").value,
    replicating: $("replicating").checked,
  };
}

async function runDesign() {
  const $ = (id) => document.getElementById(id);
  const btn = $("run");
  btn.disabled = true;
  $("empty").hidden = true;
  $("results").hidden = true;
  $("progress").hidden = false;
  $("pstate").textContent = "Working…";
  $("pstate").className = "";
  $("pbar").classList.add("ind");
  $("pstages").innerHTML =
    '<div class="stage"><span>&#9654; Ranking epitopes</span>'
    + '<div class="bar thin ind"><i></i></div></div>'
    + '<div class="stage"><span>&middot; Simulating the response</span>'
    + '<div class="bar thin"><i></i></div></div>';
  const t0 = Date.now();
  const timer = setInterval(() => {
    $("ptime").textContent = Math.round((Date.now() - t0) / 1000) + " s elapsed";
  }, 500);
  try {
    const j = await post("/api/design", designPayload());
    clearInterval(timer);
    $("progress").hidden = true;
    renderDesign(j);
  } catch (e) {
    clearInterval(timer);
    $("pbar").classList.remove("ind");
    $("pstate").textContent = "Could not run";
    $("pstate").className = "bad";
    $("pstages").innerHTML = `<div class="bad">${e.message}</div>`;
  } finally {
    btn.disabled = false;
  }
}

function renderDesign(j) {
  const $ = (id) => document.getElementById(id);
  $("results").hidden = false;
  const a = j.antigen;
  $("rtitle").textContent = a.name + (a.length ? ` · ${a.length} residues` : "");
  $("describe").textContent = j.describe || "";

  const s = j.summary;
  // Keys must match _summary() in app.py exactly. They did not the first
  // time: this read s.boost_fold, which does not exist, and so displayed a
  // real 1.41x boost as "0.00x".
  const K = [
    ["Peak IgG", fmt(s.peak_IgG), "arbitrary units, comparable only to another run"],
    ["Day of peak", Math.round(s.peak_IgG_day) + " d", "after the first dose"],
    ["Boost rise", s.boost_fold_over_baseline != null
      ? s.boost_fold_over_baseline.toFixed(2) + "×" : "—",
     "over the pre-boost level"],
    ["Peak CD8", fmt(s.peak_CTL), "effector T cells"],
    ["Memory B", fmt(s.B_memory_end), "cells at the end"],
    ["IgG still up", s.igg_durability != null
      ? Math.round(s.igg_durability * 100) + "%" : "—", "of its peak, at the end"],
    ["Affinity gain", s.affinity_fold != null
      ? s.affinity_fold.toFixed(2) + "×" : "—", "antibody affinity maturation"],
  ];
  $("kpis").innerHTML = K.map(([k, v, sub]) =>
    `<div class="kpi"><b>${v}</b><span>${k}</span>`
    + `<span class="sm" style="display:block">${sub}</span></div>`).join("");

  const d = j.day, S = j.series;
  const ser = (k) => S[k] || [];
  draw($("plot-ab"), [
    { x: d, y: ser("IgG"), name: "IgG" },
    { x: d, y: ser("IgM"), name: "IgM" },
    { x: d, y: ser("affinity"), name: "affinity", yaxis: "y2" },
  ], Object.assign(layout("Antibody", "arbitrary units"), {
    yaxis2: { overlaying: "y", side: "right", title: "affinity",
              showgrid: false },
  }));
  draw($("plot-cell"), [
    { x: d, y: ser("TC_active"), name: "CD8 effector" },
    { x: d, y: ser("TH_active"), name: "CD4 effector" },
    { x: d, y: ser("B_active"), name: "B cell" },
    { x: d, y: ser("TC_memory"), name: "CD8 memory", line: { dash: "dot" } },
    { x: d, y: ser("B_memory"), name: "B memory", line: { dash: "dot" } },
  ], layout("Cells", "cells"));
  draw($("plot-ag"), [
    { x: d, y: ser("antigen"), name: "antigen (free)" },
    { x: d, y: ser("antigen_depot"), name: "antigen on FDC" },
    { x: d, y: ser("template"), name: "template (mRNA/vector)" },
    { x: d, y: ser("danger"), name: "danger signal" },
  ], layout("Antigen and danger signal", "arbitrary units"));

  const ep = j.epitopes;
  if (!ep.rows || !ep.rows.length) {
    $("epcard").hidden = true;
    return;
  }
  $("epcard").hidden = false;
  const byClass = Object.entries(ep.by_class || {})
    .map(([k, v]) => `${v.toLocaleString()} class ${k}`).join(", ");
  $("epnote").textContent =
    `${(ep.n_binders || 0).toLocaleString()} binders (${byClass}) over `
    + `${ep.n_alleles} alleles · saved as results/${ep.saved_as}`;
  const head = ep.columns.map((c) => `<th>${c}</th>`).join("");
  const body = ep.rows.map((r) => "<tr>" + r.map((v, i) =>
    `<td class="${ep.columns[i] === "peptide" ? "seq" : (typeof v === "number" ? "num" : "")}">`
    + `${v === null ? "" : v}</td>`).join("") + "</tr>").join("");
  $("eptable").innerHTML = `<table><thead><tr>${head}</tr></thead><tbody>${body}</tbody></table>`;
}

async function downloadCsv(payload) {
  const r = await fetch("/api/timecourse.csv", {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify(payload),
  });
  if (!r.ok) return;
  const blob = await r.blob();
  const u = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = u; a.download = "timecourse.csv"; a.click();
  URL.revokeObjectURL(u);
}

window.initDesign = initDesign;
