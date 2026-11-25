/* immune_sim GUI front end.
   Plotly is themed from the CSS variables so the charts follow light/dark
   instead of sitting on a white rectangle in a dark page. */
"use strict";

const css = (n) => getComputedStyle(document.documentElement)
  .getPropertyValue(n).trim();

const PALETTE = ["#2b6cb0", "#1a7f4b", "#b4530a", "#7a4fbf",
                 "#0d8b8b", "#b03060", "#6b7280", "#8a6d00"];

function layout(title, yTitle, opts) {
  const ink = css("--ink"), muted = css("--muted"), grid = css("--grid");
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
