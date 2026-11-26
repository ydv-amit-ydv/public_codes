/* Running the project's own scripts from the browser.
 *
 * Every card here starts a subprocess on the server and streams its output
 * back. Nothing is simulated in JavaScript: what you read in the log is what
 * the script printed, so a number on screen is reproducible from a terminal.
 */
"use strict";

let STATUS = null;
let GROUP = null;
const POLLING = new Map();          // jobId -> interval handle

const fmtBytes = (n) => {
  if (n == null) return "";
  const u = ["B", "KB", "MB", "GB"];
  let i = 0, v = n;
  while (v >= 1024 && i < u.length - 1) { v /= 1024; i += 1; }
  return `${v < 10 && i > 0 ? v.toFixed(1) : Math.round(v)} ${u[i]}`;
};

const fmtAgo = (t) => {
  if (!t) return "";
  const s = Date.now() / 1000 - t;
  if (s < 90) return "just now";
  if (s < 5400) return `${Math.round(s / 60)} min ago`;
  if (s < 172800) return `${Math.round(s / 3600)} h ago`;
  return `${Math.round(s / 86400)} d ago`;
};

const el = (tag, cls, text) => {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (text != null) n.textContent = text;
  return n;
};

async function getStatus() {
  const r = await fetch("/api/status");
  if (!r.ok) throw new Error("status failed");
  STATUS = await r.json();
  return STATUS;
}

/* ---------------------------------------------------------------- forms */
function fieldInput(task, f) {
  const wrap = el("label", "field");
  wrap.appendChild(el("span", "flabel", f.label));
  let input;
  if (f.type === "bool") {
    input = el("input");
    input.type = "checkbox";
    input.checked = !!f.default;
    wrap.classList.add("checkbox");
  } else if (f.type === "choice") {
    input = el("select");
    (f.options || []).forEach((o) => {
      const opt = el("option", null, o);
      opt.value = o;
      if (o === f.default) opt.selected = true;
      input.appendChild(opt);
    });
  } else if (f.type === "params") {
    input = el("select");
    const blank = el("option", null, "— literature defaults —");
    blank.value = "";
    input.appendChild(blank);
    (STATUS.param_files || []).forEach((o) => {
      const opt = el("option", null, o.replace("models/", ""));
      opt.value = o;
      if (o === f.default) opt.selected = true;
      input.appendChild(opt);
    });
  } else {
    input = el("input");
    input.type = "number";
    input.value = f.default;
    if (f.min != null) input.min = f.min;
    if (f.max != null) input.max = f.max;
  }
  input.dataset.field = f.name;
  input.dataset.ftype = f.type;
  wrap.appendChild(input);
  if (f.help) wrap.appendChild(el("span", "fhelp", f.help));
  return wrap;
}

function readArgs(card) {
  const args = {};
  card.querySelectorAll("[data-field]").forEach((i) => {
    if (i.dataset.ftype === "bool") args[i.dataset.field] = i.checked;
    else if (i.dataset.ftype === "int") args[i.dataset.field] = Number(i.value);
    else args[i.dataset.field] = i.value;
  });
  return args;
}

/* ----------------------------------------------------------------- cards */
function outputLine(o) {
  const row = el("div", "out");
  if (!o.exists) {
    row.appendChild(el("span", "pill none", "not built"));
    row.appendChild(el("code", null, o.path));
    return row;
  }
  const stale = o.stale && o.stale.length;
  row.appendChild(el("span", `pill ${stale ? "stale" : "ok"}`,
                     stale ? "stale" : "ready"));
  row.appendChild(el("code", null, o.path));
  const bits = [fmtBytes(o.bytes)];
  if (o.rows != null) bits.push(`${o.rows.toLocaleString()} rows`);
  bits.push(fmtAgo(o.mtime));
  row.appendChild(el("span", "dim", bits.join(" · ")));
  if (stale) {
    row.appendChild(el("span", "warn",
      `rebuilt since: ${o.stale.join(", ")} — re-run to make this current`));
  }
  return row;
}

function progressBlock() {
  /* The engine's progress idiom: state + elapsed/ETA on one line, an overall
     bar, a row per step with its own bar, then the script's own log. */
  const wrap = el("div", "progress");
  wrap.hidden = true;
  const head = el("div", "thead");
  head.appendChild(el("b", "pstate", "Working…"));
  head.appendChild(el("span", "muted ptime"));
  wrap.appendChild(head);
  const bar = el("div", "bar");
  bar.style.margin = "8px 0";
  bar.appendChild(el("i"));
  wrap.appendChild(bar);
  wrap.appendChild(el("div", "pstages"));
  const det = el("details");
  det.appendChild(el("summary", "muted sm", "Log"));
  det.appendChild(el("pre", "log"));
  wrap.appendChild(det);
  return wrap;
}

function taskCard(spec) {
  const card = el("section", "card task");
  card.dataset.task = spec.name;

  const head = el("div", "thead");
  head.appendChild(el("h3", "t", spec.label));
  const tags = el("div");
  tags.appendChild(el("span", "tag plain", spec.est));
  if (spec.network) tags.appendChild(el("span", "tag net", "downloads"));
  head.appendChild(tags);
  card.appendChild(head);
  card.appendChild(el("p", "sm muted", spec.about));

  if (spec.fields && spec.fields.length) {
    const form = el("div", "grid2");
    form.style.margin = "8px 0";
    spec.fields.forEach((f) => fieldRow(form, f));
    card.appendChild(form);
  }

  if (spec.outputs && spec.outputs.length) {
    const outs = el("div");
    spec.outputs.forEach((o) => outs.appendChild(outputLine(o)));
    card.appendChild(outs);
  }

  card.appendChild(actionBar(card, spec, () => startJob(spec.name, card)));
  card.appendChild(progressBlock());
  return card;
}

function planCard(spec) {
  /* A plan is several scripts in order. Its steps are listed before it runs,
     so one click is still a visible plan rather than a black box. */
  const card = el("section", "card task");
  card.dataset.task = spec.name;

  const head = el("div", "thead");
  head.appendChild(el("h3", "t", spec.label));
  const tags = el("div");
  tags.appendChild(el("span", "tag plain", spec.est));
  head.appendChild(tags);
  card.appendChild(head);
  card.appendChild(el("p", "sm muted", spec.about));

  const list = el("ol", "sm muted");
  list.style.margin = "6px 0 8px 18px";
  spec.steps.forEach((st) => list.appendChild(el("li", null, st.label)));
  card.appendChild(list);

  if (spec.outputs && spec.outputs.length) {
    const outs = el("div");
    spec.outputs.forEach((o) => outs.appendChild(outputLine(o)));
    card.appendChild(outs);
  }

  card.appendChild(actionBar(card, spec, () => startPlan(spec.name, card),
                             spec.name.startsWith("setup")));
  card.appendChild(progressBlock());
  return card;
}

function actionBar(card, spec, onRun, primary) {
  const bar = el("div");
  bar.style.cssText = "display:flex;align-items:center;gap:8px;flex-wrap:wrap;margin-top:8px";
  const run = el("button", primary ? "primary" : "b", "Run");
  const cancel = el("button", "b cancel", "Cancel");
  cancel.hidden = true;
  bar.appendChild(run);
  bar.appendChild(cancel);
  bar.appendChild(el("span", "muted sm jstatus"));
  run.addEventListener("click", onRun);
  cancel.addEventListener("click", () => {
    const id = card.dataset.job;
    if (id) fetch(`/api/jobs/${id}/cancel`, { method: "POST" });
  });
  if (spec.missing && spec.missing.length) {
    run.disabled = true;
    const m = el("p", "sm warn");
    m.style.margin = "6px 0 0";
    m.textContent = `Needs ${spec.missing.join(", ")} first. `;
    const link = el("a", null, "Go to Data →");
    link.href = "/data";
    m.appendChild(link);
    bar.appendChild(m);
  }
  card.dataset.primary = primary ? "1" : "";
  return bar;
}

function fieldRow(form, f) {
  const lab = el("span", null, f.label);
  if (f.help) {
    const q = el("span", "q", "?");
    q.dataset.tip = f.help;
    q.tabIndex = 0;
    lab.appendChild(q);
  }
  form.appendChild(lab);
  let input;
  if (f.type === "bool") {
    input = el("input");
    input.type = "checkbox";
    input.checked = !!f.default;
  } else if (f.type === "choice") {
    input = el("select");
    (f.options || []).forEach((o) => {
      const opt = el("option", null, o);
      opt.value = o;
      if (o === f.default) opt.selected = true;
      input.appendChild(opt);
    });
  } else if (f.type === "params") {
    input = el("select");
    const blank = el("option", null, "— literature defaults —");
    blank.value = "";
    input.appendChild(blank);
    (STATUS.param_files || []).forEach((o) => {
      const opt = el("option", null, o.replace("models/", ""));
      opt.value = o;
      if (o === f.default) opt.selected = true;
      input.appendChild(opt);
    });
  } else {
    input = el("input");
    input.type = "number";
    input.value = f.default;
    if (f.min != null) input.min = f.min;
    if (f.max != null) input.max = f.max;
  }
  input.dataset.field = f.name;
  input.dataset.ftype = f.type;
  form.appendChild(input);
}

/* ------------------------------------------------------------ job driving */
async function startPlan(name, card) {
  const btn = card.querySelector("button.primary, button.b");
  const status = card.querySelector(".jstatus");
  btn.disabled = true;
  status.textContent = "starting…";
  let r;
  try {
    r = await fetch("/api/plans", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ plan: name }),
    });
  } catch (e) {
    status.textContent = "could not reach the server";
    btn.disabled = false;
    return;
  }
  const body = await r.json().catch(() => ({}));
  if (!r.ok) {
    status.className = "muted sm jstatus bad";
    status.textContent = body.error || `failed (${r.status})`;
    btn.disabled = false;
    return;
  }
  attach(card, body.id);
}

async function startJob(name, card) {
  const run = card.querySelector("button.primary, button.b");
  const status = card.querySelector(".jstatus");
  run.disabled = true;
  status.className = "muted sm jstatus";
  status.textContent = "starting…";
  let r;
  try {
    r = await fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ task: name, args: readArgs(card) }),
    });
  } catch (e) {
    status.className = "muted sm jstatus bad";
    status.textContent = "could not reach the server";
    run.disabled = false;
    return;
  }
  const body = await r.json().catch(() => ({}));
  if (!r.ok) {
    status.className = "muted sm jstatus bad";
    status.textContent = body.error || `failed (${r.status})`;
    run.disabled = false;
    return;
  }
  attach(card, body.id);
}

function attach(card, jobId) {
  card.dataset.job = jobId;
  card.querySelector(".progress").hidden = false;
  card.querySelector("button.cancel").hidden = false;
  if (POLLING.has(jobId)) return;
  const tick = async () => {
    let snap;
    try {
      const r = await fetch(`/api/jobs/${jobId}?log=200`);
      if (!r.ok) throw new Error("gone");
      snap = await r.json();
    } catch (e) {
      return;                                  // transient; try again next tick
    }
    render(card, snap);
    if (["done", "failed", "cancelled"].includes(snap.state)) {
      clearInterval(POLLING.get(jobId));
      POLLING.delete(jobId);
      card.querySelectorAll("button").forEach((b) => { b.disabled = false; });
      card.querySelector("button.cancel").hidden = true;
      refresh();                               // outputs and staleness changed
    }
  };
  POLLING.set(jobId, setInterval(tick, 1200));
  tick();
}

function render(card, snap) {
  const q = (sel) => card.querySelector(sel);
  const pstate = q(".pstate"), ptime = q(".ptime");
  const bar = q(".bar"), stages = q(".pstages"), log = q(".log");

  const pct = snap.progress == null ? 0 : Math.round(snap.progress * 100);
  const names = {
    queued: "Waiting to start…",
    running: snap.progress == null ? "Running" : `Running: ${pct}%`,
    done: "Finished", failed: "Failed", cancelled: "Cancelled",
  };
  pstate.textContent = names[snap.state] || snap.state;
  pstate.className = "pstate "
    + (snap.state === "failed" ? "bad" : snap.state === "done" ? "ok"
       : snap.state === "cancelled" ? "warn" : "");

  ptime.textContent = snap.state === "running"
    ? fmtT(snap.elapsed) + " elapsed"
      + (snap.eta != null ? `, about ${fmtT(snap.eta)} left` : "")
    : snap.state === "queued" ? "" : fmtT(snap.elapsed);

  bar.classList.toggle("ind", snap.state === "running" && !pct);
  bar.firstElementChild.style.width = pct + "%";

  // Per-step rows, so a ten-minute plan says what it is doing and what is
  // still to come, rather than showing one bar and a wall of log.
  if (snap.stages && snap.stages.length) {
    stages.innerHTML = snap.stages.map((st) => {
      const known = st.total ? Math.min(100, Math.round(100 * st.done / st.total))
        : (st.state === "done" ? 100 : 0);
      const ind = st.state === "running" && !st.total;
      const mark = st.state === "done" ? "&#10003; "
        : st.state === "running" ? "&#9654; "
        : st.state === "failed" ? "&#10007; " : "&middot; ";
      return `<div class="stage ${st.state}"><span>${mark}${esc(st.label)}</span>`
        + `<div class="bar thin ${ind ? "ind" : ""}"><i style="width:${known}%"></i></div>`
        + (st.msg ? `<div class="m">${esc(st.msg)}</div>` : "") + "</div>";
    }).join("");
  } else {
    stages.innerHTML = snap.note
      ? `<div class="stage"><span>&#9654; working</span>`
        + `<div class="bar thin ${pct ? "" : "ind"}"><i style="width:${pct}%"></i></div>`
        + `<div class="m">${esc(snap.note)}</div></div>`
      : "";
  }
  if (snap.error) {
    stages.innerHTML += `<div class="bad sm" style="margin-top:6px">`
      + `${esc(snap.error)} — the log below is the script's own output</div>`;
  }

  const atEnd = log.scrollTop + log.clientHeight >= log.scrollHeight - 20;
  log.textContent = (snap.log || []).join("\n");
  if (atEnd) log.scrollTop = log.scrollHeight;
}

const esc = (t) => String(t == null ? "" : t)
  .replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
const fmtT = (s) => s == null ? ""
  : s < 90 ? Math.round(s) + " s"
  : Math.floor(s / 60) + " min " + Math.round(s % 60) + " s";

/* ------------------------------------------------------------ candidates */
function renderCandidates() {
  const host = document.getElementById("candidates");
  if (!host) return;
  host.textContent = "";
  const cands = STATUS.candidates || [];
  const shipped = STATUS.shipped;

  const cur = el("p", "shipped");
  cur.textContent = shipped && shipped.score != null
    ? `In use: models/calibrated_wetlab.json, training score ${Number(shipped.score).toPrecision(4)}, fitted ${fmtAgo(shipped.mtime)}.`
    : "In use: models/calibrated_wetlab.json.";
  host.appendChild(cur);

  if (!cands.length) {
    host.appendChild(el("p", "hint",
      "No candidates yet. A calibration run saves one; nothing uses it until you promote it."));
    return;
  }
  host.appendChild(el("p", "hint",
    "A lower training score is not a better model. Score the candidate on the held-out "
    + "subject split first (Validation → boost kinetics, with this file selected as the "
    + "parameter set), then promote it."));

  const list = el("div", "cands");
  cands.forEach((c) => {
    const row = el("div", "cand");
    const better = shipped && shipped.score != null && c.score != null
      && c.score < shipped.score;
    row.appendChild(el("span", `pill ${better ? "ok" : "stale"}`,
      better ? "better fit" : "worse fit"));
    row.appendChild(el("code", null, c.path.replace("models/", "")));
    row.appendChild(el("span", "dim",
      `score ${c.score == null ? "?" : Number(c.score).toPrecision(4)} · ${fmtAgo(c.mtime)}`));
    const btn = el("button", "ghost", "Promote to the parameter set in use");
    btn.addEventListener("click", async () => {
      btn.disabled = true;
      const r = await fetch("/api/promote", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ file: c.path }),
      });
      const b = await r.json().catch(() => ({}));
      const msg = el("p", r.ok ? "note" : "blocked",
        r.ok ? [b.warning, b.note].filter(Boolean).join(" — ")
             : (b.error || "promotion failed"));
      row.appendChild(msg);
      btn.disabled = false;
      refresh();
    });
    row.appendChild(btn);
    list.appendChild(row);
  });
  host.appendChild(list);
}

/* ------------------------------------------------------------- inventory */
function renderInventory() {
  const host = document.getElementById("inventory");
  if (!host) return;
  host.textContent = "";
  const t = el("table");
  const thead = el("thead");
  const hr = el("tr");
  ["", "What", "Path", "Size", "Rows", "Built", "State"].forEach((h) =>
    hr.appendChild(el("th", null, h)));
  thead.appendChild(hr);
  t.appendChild(thead);
  const tb = el("tbody");
  (STATUS.inventory || []).forEach((o) => {
    const tr = el("tr");
    const dot = el("td");
    const stale = o.stale && o.stale.length;
    dot.appendChild(el("span", `dot ${!o.exists ? "none" : stale ? "stale" : "ok"}`));
    tr.appendChild(dot);
    tr.appendChild(el("td", null, o.label));
    const pc = el("td");
    pc.appendChild(el("code", null, o.path));
    tr.appendChild(pc);
    tr.appendChild(el("td", "num", o.exists ? fmtBytes(o.bytes) : "—"));
    tr.appendChild(el("td", "num", o.rows != null ? o.rows.toLocaleString() : "—"));
    tr.appendChild(el("td", null, o.exists ? fmtAgo(o.mtime) : "—"));
    tr.appendChild(el("td", "small",
      !o.exists ? "not built" :
      stale ? `stale — ${o.stale.join(", ")} is newer` : "current"));
    tb.appendChild(tr);
  });
  t.appendChild(tb);
  const wrap = el("div", "tablewrap");
  wrap.appendChild(t);
  host.appendChild(wrap);
}

function renderPlans() {
  const host = document.getElementById("plans");
  if (!host) return;
  const open = {};
  host.querySelectorAll(".task").forEach((c) => {
    if (c.dataset.job) open[c.dataset.task] = c.dataset.job;
  });
  host.textContent = "";
  Object.values(STATUS.plans || {})
    .filter((p) => p.group === GROUP)
    .forEach((p) => {
      const card = planCard(p);
      host.appendChild(card);
      if (open[p.name]) attach(card, open[p.name]);
    });
  (STATUS.jobs || []).forEach((j) => {
    if (j.state !== "running" && j.state !== "queued") return;
    const card = host.querySelector(`.task[data-task="${j.task}"]`);
    if (card && !card.dataset.job) attach(card, j.id);
  });
}

function renderData() {
  /* Where each dataset came from and whether you are using the copy that
     ships with the repository or one you harvested yourself. Citations are
     shown because redistributing someone's database carries their terms. */
  const host = document.getElementById("datainv");
  if (!host) return;
  host.textContent = "";
  const rows = STATUS.data || [];
  if (!rows.length) return;
  const t = el("table");
  const th = el("thead"), hr = el("tr");
  ["", "Dataset", "Copy in use", "Size", "Source", "Cite"].forEach((h) =>
    hr.appendChild(el("th", null, h)));
  th.appendChild(hr); t.appendChild(th);
  const tb = el("tbody");
  rows.forEach((r) => {
    const tr = el("tr");
    const dot = el("td");
    dot.appendChild(el("span", `dot ${r.present ? "ok" : "none"}`));
    tr.appendChild(dot);
    const name = el("td");
    name.appendChild(el("b", null, r.what));
    name.appendChild(el("br"));
    name.appendChild(el("code", null, r.path || `data/${r.name}.csv`));
    name.style.whiteSpace = "normal";
    tr.appendChild(name);
    tr.appendChild(el("td", null, r.copy === "bundled" ? "bundled with the repo"
      : r.copy === "harvested" ? "harvested by you" : "missing"));
    tr.appendChild(el("td", "num", r.bytes ? fmtBytes(r.bytes) : "—"));
    tr.appendChild(el("td", null, r.source));
    const cite = el("td", "sm muted", r.citation);
    cite.style.whiteSpace = "normal";
    cite.style.maxWidth = "320px";
    tr.appendChild(cite);
    tb.appendChild(tr);
  });
  t.appendChild(tb);
  const wrap = el("div", "tablewrap");
  wrap.appendChild(t);
  host.appendChild(wrap);
}

function renderTasks() {
  const host = document.getElementById("tasks");
  if (!host) return;
  const open = {};
  host.querySelectorAll(".task").forEach((c) => {
    if (c.dataset.job) open[c.dataset.task] = c.dataset.job;
  });
  host.textContent = "";
  Object.values(STATUS.tasks)
    .filter((s) => s.group === GROUP)
    .forEach((s) => {
      const card = taskCard(s);
      host.appendChild(card);
      if (open[s.name]) attach(card, open[s.name]);
    });
  // a job started before this page loaded, or on another tab
  (STATUS.jobs || []).forEach((j) => {
    if (j.state !== "running" && j.state !== "queued") return;
    const card = host.querySelector(`.task[data-task="${j.task}"]`);
    if (card && !card.dataset.job) attach(card, j.id);
  });
}

async function refresh() {
  try {
    await getStatus();
  } catch (e) {
    return;
  }
  renderTasks();
  renderPlans();
  renderData();
  renderCandidates();
  renderInventory();
}

async function initRunPage(group) {
  GROUP = group;
  await refresh();
  setInterval(() => { if (!POLLING.size) refresh(); }, 15000);
}

window.initRunPage = initRunPage;
window.jobsRefresh = refresh;
