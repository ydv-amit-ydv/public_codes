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

function taskCard(spec) {
  const card = el("section", "card task");
  card.dataset.task = spec.name;

  const head = el("div", "thead");
  head.appendChild(el("h3", null, spec.label));
  const tags = el("div", "tags");
  tags.appendChild(el("span", "tag", spec.est));
  if (spec.network) tags.appendChild(el("span", "tag net", "downloads"));
  head.appendChild(tags);
  card.appendChild(head);

  card.appendChild(el("p", "about", spec.about));

  if (spec.fields && spec.fields.length) {
    const form = el("div", "fields");
    spec.fields.forEach((f) => form.appendChild(fieldInput(spec.name, f)));
    card.appendChild(form);
  }

  if (spec.outputs && spec.outputs.length) {
    const outs = el("div", "outs");
    spec.outputs.forEach((o) => outs.appendChild(outputLine(o)));
    card.appendChild(outs);
  }

  const bar = el("div", "actions");
  const run = el("button", "primary", "Run");
  const cancel = el("button", "ghost", "Cancel");
  cancel.hidden = true;
  bar.appendChild(run);
  bar.appendChild(cancel);
  const status = el("span", "jstatus");
  bar.appendChild(status);
  card.appendChild(bar);

  if (spec.missing && spec.missing.length) {
    run.disabled = true;
    const m = el("p", "blocked",
      `Needs ${spec.missing.join(", ")} first.`);
    const link = el("a", null, "Go to Data →");
    link.href = "/data";
    m.appendChild(document.createTextNode(" "));
    m.appendChild(link);
    card.appendChild(m);
  }

  const prog = el("div", "progwrap");
  prog.hidden = true;
  const bar2 = el("div", "prog");
  const fill = el("div", "fill");
  bar2.appendChild(fill);
  prog.appendChild(bar2);
  prog.appendChild(el("div", "note"));
  card.appendChild(prog);

  const logBox = el("pre", "log");
  logBox.hidden = true;
  card.appendChild(logBox);

  run.addEventListener("click", () => startJob(spec.name, card));
  cancel.addEventListener("click", () => {
    const id = card.dataset.job;
    if (id) fetch(`/api/jobs/${id}/cancel`, { method: "POST" });
  });
  return card;
}

/* ------------------------------------------------------------ job driving */
async function startJob(name, card) {
  const run = card.querySelector("button.primary");
  const status = card.querySelector(".jstatus");
  run.disabled = true;
  status.className = "jstatus";
  status.textContent = "starting…";
  let r;
  try {
    r = await fetch("/api/jobs", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ task: name, args: readArgs(card) }),
    });
  } catch (e) {
    status.className = "jstatus bad";
    status.textContent = "could not reach the server";
    run.disabled = false;
    return;
  }
  const body = await r.json().catch(() => ({}));
  if (!r.ok) {
    status.className = "jstatus bad";
    status.textContent = body.error || `failed (${r.status})`;
    run.disabled = false;
    return;
  }
  attach(card, body.id);
}

function attach(card, jobId) {
  card.dataset.job = jobId;
  card.querySelector(".log").hidden = false;
  card.querySelector(".progwrap").hidden = false;
  card.querySelector("button.ghost").hidden = false;
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
      card.querySelector("button.primary").disabled = false;
      card.querySelector("button.ghost").hidden = true;
      refresh();                               // outputs and staleness changed
    }
  };
  POLLING.set(jobId, setInterval(tick, 1200));
  tick();
}

function render(card, snap) {
  const status = card.querySelector(".jstatus");
  const prog = card.querySelector(".progwrap");
  const fill = card.querySelector(".fill");
  const note = card.querySelector(".progwrap .note");
  const log = card.querySelector(".log");

  const mins = Math.floor(snap.elapsed / 60);
  const secs = Math.round(snap.elapsed % 60);
  const t = mins ? `${mins}m ${secs}s` : `${secs}s`;
  status.className = "jstatus " + (
    snap.state === "done" ? "good" :
    snap.state === "failed" ? "bad" :
    snap.state === "cancelled" ? "warnc" : "");
  status.textContent = `${snap.state} · ${t}`;

  if (snap.progress != null) {
    fill.style.width = `${Math.round(snap.progress * 100)}%`;
    fill.classList.remove("indet");
  } else if (snap.state === "running") {
    fill.style.width = "100%";
    fill.classList.add("indet");
  }
  note.textContent = snap.note || "";
  prog.hidden = false;

  const atBottom = log.scrollHeight - log.scrollTop - log.clientHeight < 40;
  log.textContent = (snap.log || []).join("\n");
  if (atBottom) log.scrollTop = log.scrollHeight;

  if (snap.state === "failed" && snap.error) {
    note.textContent = snap.error + " — the log above is the script's own output";
  }
}

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
