// Home Designer — application shell: state, panels, toolbar, keyboard, save/load.
import { CATEGORIES, ITEMS, ITEM_COUNT, ROOM_TYPES } from './catalog.js';
import {
  uid, newDesign, newFloor, normalizeDesign, History, polygonArea, polygonPerimeter, polygonCentroid,
  wallLength, wallAngle, openingGeom, clampOpeningT, fmtLen, fmtArea, itemCorners, pointInPolygon, dist, parseLen,
} from './model.js';
import { FLOOR_MATERIALS } from './textures.js';
import { PlanView, shade } from './plan.js';
import { Store, safeName, download } from './storage.js';

const $ = (s, el = document) => el.querySelector(s);
const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k.startsWith('on')) el.addEventListener(k.slice(2), v);
    else if (k === 'class') el.className = v;
    else if (k === 'html') el.innerHTML = v;
    else el.setAttribute(k, v === true ? '' : v);
  }
  for (const c of kids.flat()) if (c !== null && c !== undefined && c !== false) el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  return el;
};
const esc = (s) => String(s ?? '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);

const TOOLS = [
  ['select', '↖', 'Select / move (V)', 'v'],
  ['wall', '▭', 'Draw walls (W) — click to chain, double-click / right-click / Esc to stop, Shift = 45°', 'w'],
  ['roombox', '⬚', 'Room with walls (B) — drag a rectangle', 'b'],
  ['area', '▢', 'Room / area only (A) — drag a rectangle, no walls', 'a'],
  ['room', '⬠', 'Polygon room (P) — click corners, click the first point to close', 'p'],
  ['pin', '📍', 'Note pin (N) — annotate a place', 'n'],
  ['label', 'T', 'Text label (T)', 't'],
  ['dim', '↔', 'Dimension line (D)', 'd'],
  ['pan', '✋', 'Pan (H) — or hold Space / drag empty space', 'h'],
];

class App {
  constructor() {
    this.store = new Store();
    this.history = new History();
    this.design = null;
    this.floorIdx = 0;
    this.selection = null;
    this.tool = 'select';
    this.placing = null;
    this.fileName = null;   // name it was last saved/opened as
    this.dirty = false;
    this.view = 'split';
    this.view3d = null;
    this.lenBuf = '';
  }

  get floor() { return this.design.floors[this.floorIdx]; }

  async start() {
    this.buildShell();
    this.plan = new PlanView($('#plan'), this);
    await this.store.init();
    $('#where').textContent = this.store.server ? '💾 saving to designs/ folder' : '💾 saving in this browser';
    $('#where').title = this.store.server
      ? 'Designs are written as JSON files to home_designer/designs/ — commit them to git to keep them.'
      : 'No server detected: designs are kept in this browser. Run start.sh to save files into the repo.';

    const draft = this.store.loadDraft();
    if (draft && draft.design) {
      this.loadDesign(draft.design, draft.fileName || null, { keepDirty: draft.dirty });
      this.toast('Restored your last session');
    } else {
      let loaded = false;
      if (this.store.server) {
        try { const d = await this.store.load('sample_2bhk'); this.loadDesign(d, null); this.design.name = 'My home'; loaded = true; this.toast('Loaded the sample 2BHK — edit it or start a New design'); } catch { /* no sample */ }
      }
      if (!loaded) this.loadDesign(newDesign('My home'), null);
    }
    let v = 'split'; try { v = localStorage.getItem('homeDesigner:view') || 'split'; } catch { /* ignore */ }
    await this.setView(v);
    setTimeout(() => this.plan.fit(), 60);
    this.bindKeys();
    window.addEventListener('beforeunload', (e) => { if (this.dirty) { this.saveDraftNow(); } });
  }

  // properties for a newly drawn wall, from the current wall type (exterior / interior)
  newWallProps() {
    const s = this.design.settings, kind = s.wallKind || 'exterior';
    return { kind, thickness: kind === 'interior' ? s.intThickness : s.extThickness, height: this.floor.height || s.wallHeight, color: '#f1ede6' };
  }

  setWallKind(kind) {
    this.design.settings.wallKind = kind;
    const b = document.getElementById('wallkind');
    if (b) {
      const s = this.design.settings;
      b.textContent = kind === 'interior' ? 'Int' : 'Ext';
      b.title = `New walls: ${kind} (${Math.round((kind === 'interior' ? s.intThickness : s.extThickness) * 1000) / 10} cm) — click or press E to switch`;
      b.classList.toggle('int', kind === 'interior');
      const sub = document.getElementById('wallkind-th');
      if (sub) sub.textContent = Math.round((kind === 'interior' ? s.intThickness : s.extThickness) * 1000) / 10;
    }
    if (this.tool === 'wall' || this.tool === 'roombox') this.setTool(this.tool);
    this.plan?.requestRender();
  }

  // change a wall's length keeping its start point. Walls joined at the far end are moved
  // along with it (so a rectangular room simply gets bigger) and the walls beyond them stretch.
  setWallLength(w, L) {
    const a = wallAngle(w);
    const nx = w.x1 + Math.cos(a) * L, ny = w.y1 + Math.sin(a) * L;
    const ox = w.x2, oy = w.y2, dx = nx - ox, dy = ny - oy;
    const near = (x, y, px, py) => dist(x, y, px, py) < 1e-3;
    const moved = [[ox, oy, nx, ny]];
    const shifted = new Set([w]);
    for (const o of this.floor.walls) {
      if (o === w) continue;
      const a1 = near(o.x1, o.y1, ox, oy), a2 = near(o.x2, o.y2, ox, oy);
      if (!a1 && !a2) continue;
      const parallel = Math.abs(Math.sin(wallAngle(o) - a)) < 0.02;
      if (parallel) { if (a1) { o.x1 = nx; o.y1 = ny; } else { o.x2 = nx; o.y2 = ny; } shifted.add(o); continue; }
      const far = a1 ? [o.x2, o.y2] : [o.x1, o.y1];
      moved.push([far[0], far[1], far[0] + dx, far[1] + dy]);
      o.x1 += dx; o.y1 += dy; o.x2 += dx; o.y2 += dy;
      shifted.add(o);
    }
    for (const [px, py, qx, qy] of moved.slice(1)) {
      for (const o of this.floor.walls) {
        if (shifted.has(o)) continue;
        if (near(o.x1, o.y1, px, py)) { o.x1 = qx; o.y1 = qy; }
        if (near(o.x2, o.y2, px, py)) { o.x2 = qx; o.y2 = qy; }
      }
    }
    w.x2 = nx; w.y2 = ny;
    for (const r of this.floor.rooms) r.points = r.points.map(([x, y]) => { const m = moved.find(([px, py]) => near(x, y, px, py)); return m ? [m[2], m[3]] : [x, y]; });
    for (const o of this.floor.openings) { const ww = this.floor.walls.find((x) => x.id === o.wallId); if (ww) o.t = clampOpeningT(o, ww); }
  }

  // mark each wall exterior if one side of it is outside every indoor room
  autoWallKinds() {
    const OUT = /garden|lawn|parking|garage|driveway|sit-out|porch|terrace|balcony|verandah|courtyard|open area/i;
    const indoor = this.floor.rooms.filter((r) => !OUT.test(r.type || '') && !OUT.test(r.name || ''));
    const s = this.design.settings;
    let ext = 0, int = 0;
    for (const w of this.floor.walls) {
      const a = wallAngle(w), nx = -Math.sin(a), ny = Math.cos(a);
      const L = wallLength(w), off = (w.thickness || 0.15) / 2 + 0.1;
      let both = true;
      for (const t of [0.25, 0.5, 0.75]) {
        const mx = w.x1 + (w.x2 - w.x1) * t, my = w.y1 + (w.y2 - w.y1) * t;
        const inA = indoor.some((r) => pointInPolygon(mx + nx * off, my + ny * off, r.points));
        const inB = indoor.some((r) => pointInPolygon(mx - nx * off, my - ny * off, r.points));
        if (!(inA && inB)) both = false;
      }
      if (L < 1e-6) continue;
      w.kind = both ? 'interior' : 'exterior';
      w.thickness = both ? s.intThickness : s.extThickness;
      both ? int++ : ext++;
    }
    this.commit('Auto wall types', { props: true });
    this.toast(`${ext} exterior and ${int} interior walls`);
  }

  applyThickness(kind) {
    const s = this.design.settings;
    for (const f of this.design.floors) for (const w of f.walls) if ((w.kind || 'exterior') === kind) w.thickness = kind === 'interior' ? s.intThickness : s.extThickness;
    this.setWallKind(s.wallKind || 'exterior');
  }

  // ---------- state ----------
  loadDesign(d, fileName, { keepDirty = false } = {}) {
    this.design = normalizeDesign(JSON.parse(JSON.stringify(d)));
    this.fileName = fileName;
    this.floorIdx = 0;
    this.selection = null;
    this.history.reset(this.snap());
    this.dirty = keepDirty;
    this.refreshAll();
    if (this.plan) { this.plan._fitted = true; this.plan.fit(); }
    if (this.view3d) { this.view3d.fitted = false; this.view3d.markDirty(); }
  }

  commit(label = 'Edit', { props = false } = {}) {
    this.design.updated = new Date().toISOString();
    const snap = this.snap();
    // compare without the timestamp so no-op drags do not create history steps
    const last = this.history.past[this.history.past.length - 1];
    if (last && stripTime(last) === stripTime(snap)) return;
    this.history.push(snap);
    this.dirty = true;
    this.afterChange({ props });
    this.status(label);
  }

  liveUpdate() {
    this.plan.requestRender();
    if (this.view3d) this.view3d.markDirty();
    this.refreshPropsValues();
  }

  afterChange({ props = true } = {}) {
    this.validateSelection();
    this.plan.requestRender();
    if (this.view3d) this.view3d.markDirty();
    this.refreshLists();
    if (props) this.refreshProps();
    this.updateTitle();
    clearTimeout(this._draftT);
    this._draftT = setTimeout(() => this.saveDraftNow(), 400);
  }

  saveDraftNow() {
    const ok = this.store.saveDraft({ design: this.design, fileName: this.fileName, dirty: this.dirty });
    // browser storage is small (~5 MB): retry without embedded images if needed
    if (!ok) this.store.saveDraft({ design: { ...this.design, assets: {} }, fileName: this.fileName, dirty: this.dirty });
  }

  // ---------- trace image (underlay) ----------
  async loadUnderlayFile(file) {
    if (!file) return;
    if (!file.type.startsWith('image/')) { alert('Please choose an image (JPG, PNG, WebP…). For a PDF, take a screenshot of the page first.'); return; }
    const url = await new Promise((res, rej) => { const fr = new FileReader(); fr.onload = () => res(fr.result); fr.onerror = rej; fr.readAsDataURL(file); });
    await this.setUnderlay(url, file.name);
  }

  async loadUnderlayURL() {
    const url = prompt('Paste the web address (URL) of an image — a floor plan, sketch or photo:');
    if (!url) return;
    let data = null;
    try {
      if (this.store.server) {
        const r = await fetch('api/fetch-image?url=' + encodeURIComponent(url));
        if (!r.ok) throw new Error(await r.text());
        data = (await r.json()).dataUrl;
      } else {
        const r = await fetch(url);
        const blob = await r.blob();
        data = await new Promise((res) => { const fr = new FileReader(); fr.onload = () => res(fr.result); fr.readAsDataURL(blob); });
      }
    } catch (err) {
      alert('Could not download that image: ' + err.message + '\nTip: save the image to your computer and use “Upload image” instead.');
      return;
    }
    await this.setUnderlay(data, url.split('/').pop().slice(0, 60));
  }

  async setUnderlay(dataUrl, name) {
    // shrink big photos so design files stay small
    const img = await new Promise((res, rej) => { const i = new Image(); i.onload = () => res(i); i.onerror = () => rej(new Error('not an image')); i.src = dataUrl; });
    const maxPx = 2400, k = Math.min(1, maxPx / Math.max(img.naturalWidth, img.naturalHeight));
    const cv = document.createElement('canvas');
    cv.width = Math.round(img.naturalWidth * k); cv.height = Math.round(img.naturalHeight * k);
    cv.getContext('2d').drawImage(img, 0, 0, cv.width, cv.height);
    const small = cv.toDataURL('image/jpeg', 0.85);
    const id = uid('img');
    this.design.assets = this.design.assets || {};
    const old = this.floor.underlay;
    if (old && old.src?.startsWith('asset:') && !this.design.floors.some((f) => f !== this.floor && f.underlay?.src === old.src)) delete this.design.assets[old.src.slice(6)];
    this.design.assets[id] = small;
    // first guess: image is 15 m wide; calibrate next
    const b = this.plan.bounds();
    this.floor.underlay = { src: 'asset:' + id, name: name || 'image', x: b ? b.x0 : 0, y: b ? b.y0 : 0, mpp: 15 / cv.width, opacity: 0.5, visible: true };
    this.commit('Trace image', { props: true });
    this.plan.fit();
    const u = this.floor.underlay;
    this.plan.focusOn(u.x + (cv.width * u.mpp) / 2, u.y + (cv.height * u.mpp) / 2);
    this.setTool('calib');
    this.toast('Now set the scale: click both ends of something whose length you know (e.g. a wall), then type its length', 6000);
  }

  calibrateUnderlay(a, b) {
    const u = this.floor.underlay;
    this.setTool('select');
    if (!u) return;
    const dw = dist(a[0], a[1], b[0], b[1]);
    if (dw < 1e-6) return;
    const ans = prompt('Real length between the two points? (e.g. 4.5, 450cm, 15ft, 14\'6)', fmtLen(dw, this.design.settings.units).replace(/ m$/, ''));
    if (ans === null) return;
    const L = parseLen(ans.replace(/\s+/g, '').replace(/′/g, "'").replace(/″/g, '"'), this.design.settings.units);
    if (!(L > 0)) { alert('Could not read that length'); return; }
    const f = L / dw;
    u.mpp *= f;
    u.x = a[0] + (u.x - a[0]) * f;
    u.y = a[1] + (u.y - a[1]) * f;
    this.commit('Calibrate trace image', { props: true });
    this.plan.fit();
    this.toast('Scale set. Now draw walls over the image (W). Move it with “Move image”.');
  }

  undo() { const s = this.history.undo(); if (s) this.restore(s, 'Undo'); }
  redo() { const s = this.history.redo(); if (s) this.restore(s, 'Redo'); }
  // undo snapshots leave out embedded images (design.assets) so history stays small
  snap() { return JSON.stringify(this.design, (k, v) => (k === 'assets' ? undefined : v)); }
  assetSrc(src) { return src && src.startsWith('asset:') ? this.design.assets?.[src.slice(6)] : src; }

  restore(snap, label) {
    const assets = this.design.assets;
    this.design = JSON.parse(snap);
    if (assets) this.design.assets = assets;
    this.floorIdx = Math.min(this.floorIdx, this.design.floors.length - 1);
    this.dirty = true;
    this.refreshAll();
    this.status(label);
  }

  refreshAll() {
    this.validateSelection();
    this.refreshFloors();
    this.refreshLists();
    this.refreshProps();
    this.updateTitle();
    this.plan?.requestRender();
    this.view3d?.markDirty();
    $('#design-name').value = this.design.name;
    $('#units').value = this.design.settings.units;
    this.setWallKind(this.design.settings.wallKind || 'exterior');
  }

  validateSelection() {
    if (this.selection && !this.getSelected()) this.selection = null;
  }

  coll(type) {
    return { wall: this.floor.walls, opening: this.floor.openings, item: this.floor.items, room: this.floor.rooms, note: this.floor.notes }[type];
  }

  getSelected() {
    if (!this.selection) return null;
    const c = this.coll(this.selection.type);
    return c ? c.find((o) => o.id === this.selection.id) || null : null;
  }

  select(type, id, { focus = false } = {}) {
    this.selection = type ? { type, id } : null;
    this.refreshProps();
    this.refreshLists();
    this.plan.requestRender();
    if (focus) this.focusProps();
  }

  focusProps() {
    setTimeout(() => { const el = $('#props input[type=text], #props textarea'); if (el) { el.focus(); el.select?.(); } }, 30);
  }

  setTool(tool, placing = null) {
    this.tool = tool;
    this.placing = placing;
    if (this.plan) { this.plan.wallStart = null; this.plan.roomPts = []; this.plan.dimStart = null; this.plan.boxStart = null; this.plan.requestRender(); }
    document.querySelectorAll('#tools button').forEach((b) => b.classList.toggle('active', b.dataset.tool === tool));
    document.querySelectorAll('.cat-item').forEach((b) => b.classList.toggle('active', !!placing && b.dataset.id === placing.id));
    const hints = {
      select: 'Click to select · drag to move · double-click to edit · drag empty space to pan · wheel to zoom',
      wall: `Drawing ${this.design?.settings.wallKind || 'exterior'} walls · click to start, click for each corner · TYPE a length (e.g. 4.5) + Enter for exact walls · E = switch Ext/Int · Shift = 45° · Esc to finish`,
      roombox: 'Drag a rectangle to create a room with four walls — or type a size like 4x3.5 and press Enter',
      area: 'Drag a rectangle to mark a room / place (no walls) — e.g. garden, parking, open kitchen',
      room: 'Click corners of the room · click the first point (or Enter) to close',
      pin: 'Click anywhere to drop a numbered note pin',
      label: 'Click to place a text label',
      dim: 'Click two points to add a dimension line',
      pan: 'Drag to pan',
      calib: 'Calibrate: click the two ends of a known length on the image (e.g. a wall)',
      movebg: 'Drag the trace image to line it up · Esc when done',
      place: placing ? `Placing “${placing.name}” — click to place · R rotates · Shift+click keeps placing · Esc cancels${placing.kind === 'opening' ? ' · must be on a wall' : placing.wall ? ' · snaps to walls (Alt to place freely)' : ''}` : '',
    };
    $('#hint').textContent = hints[tool] || '';
  }

  status(t) { $('#status-left').textContent = t; }

  toast(msg, ms = 2600) {
    const t = h('div', { class: 'toast' }, msg);
    $('#toasts').append(t);
    setTimeout(() => t.classList.add('show'), 10);
    setTimeout(() => { t.classList.remove('show'); setTimeout(() => t.remove(), 300); }, ms);
  }

  updateTitle() {
    const n = this.design.name + (this.dirty ? ' •' : '');
    document.title = `${n} — Home Designer`;
    $('#save-state').textContent = this.dirty ? 'unsaved changes' : this.fileName ? `saved as ${this.fileName}` : '';
  }

  // ---------- shell ----------
  buildShell() {
    // tools
    const tools = $('#tools');
    for (const [id, icon, title] of TOOLS) {
      tools.append(h('button', { 'data-tool': id, title, onclick: () => this.setTool(id) }, h('span', {}, icon)));
    }
    tools.append(h('div', { class: 'sep' }));
    tools.append(h('button', { id: 'wallkind', class: 'wallkind', onclick: () => this.setWallKind(this.design.settings.wallKind === 'interior' ? 'exterior' : 'interior') }, 'Ext'));
    tools.append(h('div', { id: 'wallkind-th', class: 'wallkind-th' }, ''));
    tools.append(h('div', { class: 'sep' }));
    tools.append(h('button', { title: 'Zoom in (+)', onclick: () => this.plan.zoomBy(1.25) }, '+'));
    tools.append(h('button', { title: 'Zoom out (−)', onclick: () => this.plan.zoomBy(0.8) }, '−'));
    tools.append(h('button', { title: 'Fit to drawing (F)', onclick: () => this.plan.fit() }, '⤢'));

    // top bar
    $('#btn-new').onclick = () => this.newDesign();
    $('#btn-open').onclick = () => this.openDialog();
    $('#btn-save').onclick = () => this.save();
    $('#btn-saveas').onclick = () => this.save(true);
    $('#btn-undo').onclick = () => this.undo();
    $('#btn-redo').onclick = () => this.redo();
    $('#btn-help').onclick = () => $('#help').showModal();
    $('#btn-import').onclick = () => $('#file-import').click();
    $('#file-import').onchange = (e) => this.importFile(e.target.files[0]);
    $('#file-underlay').onchange = (e) => { this.loadUnderlayFile(e.target.files[0]); e.target.value = ''; };
    $('#btn-trace').onclick = () => $('#file-underlay').click();
    $('#file-model').onchange = (e) => { this.uploadModel(e.target.files[0]); e.target.value = ''; };
    $('#exp-json').onclick = () => download(safeName(this.design.name) + '.json', JSON.stringify(this.design, null, 1));
    $('#exp-png').onclick = () => download(safeName(this.design.name) + `_${safeName(this.floor.name)}_plan.png`, this.plan.exportPNG());
    $('#exp-3d').onclick = async () => { await this.ensure3D(); this.setView(this.view === '2d' ? 'split' : this.view); setTimeout(() => download(safeName(this.design.name) + '_3d.png', this.view3d.snapshot(2)), 600); };
    $('#exp-csv').onclick = () => download(safeName(this.design.name) + '_schedule.csv', this.scheduleCSV(), 'text/csv');
    $('#exp-print').onclick = () => this.printSheet();
    $('#design-name').onchange = (e) => { this.design.name = e.target.value.trim() || 'Untitled design'; this.commit('Rename design'); };
    $('#units').onchange = (e) => { this.design.settings.units = e.target.value; this.commit('Units', { props: true }); };
    document.querySelectorAll('[data-view]').forEach((b) => (b.onclick = () => this.setView(b.dataset.view)));
    $('#floor-sel').onchange = (e) => this.switchFloor(+e.target.value);
    $('#btn-addfloor').onclick = () => this.addFloor();

    // left tabs
    document.querySelectorAll('#left .tab').forEach((t) => (t.onclick = () => {
      document.querySelectorAll('#left .tab').forEach((x) => x.classList.toggle('active', x === t));
      document.querySelectorAll('#left .pane').forEach((p) => p.classList.toggle('active', p.id === t.dataset.pane));
      this.refreshLists();
    }));
    this.buildCatalog();
    $('#cat-search').oninput = () => this.filterCatalog();
    $('#places-search').oninput = () => this.refreshLists();
    $('#sched-all').onchange = () => this.refreshLists();

    // 3D controls
    $('#opt-allfloors').onchange = (e) => { this.view3d.opts.allFloors = e.target.checked; this.view3d.markDirty(); };
    $('#opt-cutaway').onchange = (e) => { this.view3d.opts.cutaway = e.target.checked; this.view3d.markDirty(); };
    $('#opt-labels').onchange = (e) => { this.view3d.opts.labels = e.target.checked; this.view3d.markDirty(); };
    $('#opt-hq').onchange = (e) => { this.view3d.opts.hq = e.target.checked; try { localStorage.setItem('homeDesigner:hq', e.target.checked ? '1' : '0'); } catch { /* ignore */ } };
    $('#opt-sun').oninput = (e) => { this.view3d.opts.sunAngle = +e.target.value; this.view3d.fitSun(); };
    $('#btn-walk').onclick = () => this.view3d.enterWalk();
    $('#btn-shot3d').onclick = () => download(safeName(this.design.name) + '_3d.png', this.view3d.snapshot(2));
    $('#btn-fit3d').onclick = () => this.view3d.fit();
    $('#btn-top3d').onclick = () => this.view3d.topView();
  }

  async ensure3D() {
    if (this.view3d) return;
    try {
      const { View3D } = await import('./view3d.js');
      this.view3d = new View3D($('#view3d'), this);
      let hq = null; try { hq = localStorage.getItem('homeDesigner:hq'); } catch { /* ignore */ }
      if (hq === '0') { this.view3d.opts.hq = false; $('#opt-hq').checked = false; }
    } catch (err) {
      console.error(err);
      $('#view3d').innerHTML = `<div class="err">3D view could not start: ${esc(err.message)}</div>`;
    }
  }

  async setView(v) {
    this.view = v;
    try { localStorage.setItem('homeDesigner:view', v); } catch { /* ignore */ }
    document.body.dataset.view = v;
    document.querySelectorAll('[data-view]').forEach((b) => b.classList.toggle('active', b.dataset.view === v));
    if (v !== '2d') await this.ensure3D();
    this.view3d?.setVisible(v !== '2d');
    setTimeout(() => { this.plan.resize(); this.view3d?.resize(); }, 30);
  }

  // ---------- catalog ----------
  buildCatalog() {
    const root = $('#catalog');
    root.innerHTML = '';
    $('#cat-count').textContent = `${ITEM_COUNT} items`;
    for (const cat of CATEGORIES) {
      const body = h('div', { class: 'cat-body' });
      const det = h('details', { class: 'cat', 'data-cat': cat.name }, h('summary', {}, `${cat.icon} ${cat.name}`, h('span', { class: 'count' }, cat.items.length)), body);
      for (const it of cat.items) {
        const btn = h('button', {
          class: 'cat-item', 'data-id': it.id, title: `${it.name} — ${Math.round(it.w * 100)}×${Math.round(it.d * 100)}×${Math.round(it.h * 100)} cm`,
          onclick: () => { if (this.placing?.id === it.id) { this.setTool('select'); return; } this.plan.ghostRot = 0; this.setTool('place', it); },
        }, h('canvas', { width: 88, height: 88 }), h('span', {}, it.name));
        body.append(btn);
      }
      det.addEventListener('toggle', () => { if (det.open) this.drawThumbs(det); });
      root.append(det);
    }
    root.firstChild.open = true;
  }

  drawThumbs(scope) {
    scope.querySelectorAll('.cat-item').forEach((b) => {
      if (b.dataset.drawn) return;
      b.dataset.drawn = '1';
      const def = ITEMS[b.dataset.id];
      const cv = b.querySelector('canvas');
      drawThumb(cv, def);
    });
  }

  filterCatalog() {
    const q = $('#cat-search').value.trim().toLowerCase();
    document.querySelectorAll('#catalog .cat').forEach((det) => {
      let n = 0;
      det.querySelectorAll('.cat-item').forEach((b) => {
        const def = ITEMS[b.dataset.id];
        const ok = !q || def.name.toLowerCase().includes(q) || def.cat.toLowerCase().includes(q) || def.id.includes(q);
        b.style.display = ok ? '' : 'none';
        if (ok) n++;
      });
      det.style.display = n ? '' : 'none';
      if (q && n) { det.open = true; this.drawThumbs(det); }
    });
  }

  // ---------- floors ----------
  refreshFloors() {
    const sel = $('#floor-sel');
    sel.innerHTML = '';
    this.design.floors.forEach((f, i) => sel.append(h('option', { value: i, selected: i === this.floorIdx }, f.name)));
  }

  switchFloor(i) {
    this.floorIdx = i;
    this.selection = null;
    this.refreshAll();
    this.plan.requestRender();
  }

  addFloor() {
    const top = this.design.floors.reduce((a, f) => Math.max(a, (f.elevation || 0) + (f.height || 3)), 0);
    const n = this.design.floors.length;
    const names = ['Ground floor', 'First floor', 'Second floor', 'Third floor', 'Fourth floor', 'Terrace'];
    const f = newFloor(names[n] || `Floor ${n}`, +(top + 0.15).toFixed(2), this.design.settings.wallHeight);
    if (confirm('Copy the outer walls of the current floor to the new floor?')) {
      f.walls = this.floor.walls.map((w) => ({ ...w, id: uid('w') }));
    }
    this.design.floors.push(f);
    this.floorIdx = this.design.floors.length - 1;
    this.commit('Add floor');
    this.refreshAll();
  }

  deleteFloor() {
    if (this.design.floors.length < 2) { this.toast('A design needs at least one floor'); return; }
    if (!confirm(`Delete “${this.floor.name}” and everything on it?`)) return;
    this.design.floors.splice(this.floorIdx, 1);
    this.floorIdx = Math.max(0, this.floorIdx - 1);
    this.commit('Delete floor');
    this.refreshAll();
  }

  // ---------- left panel lists ----------
  refreshLists() {
    clearTimeout(this._listT);
    this._listT = setTimeout(() => { this.renderPlaces(); this.renderSchedule(); }, 30);
  }

  renderPlaces() {
    const root = $('#places');
    if (!root || !$('#pane-places').classList.contains('active')) return;
    const q = $('#places-search').value.trim().toLowerCase();
    const f = this.floor, units = this.design.settings.units;
    const match = (...s) => !q || s.some((x) => String(x || '').toLowerCase().includes(q));
    const sel = this.selection;
    const row = (type, id, title, sub, color, extra = '') => h('div', {
      class: 'place-row' + (sel && sel.id === id ? ' sel' : ''),
      onclick: () => { this.select(type, id); this.focusObj(type, id); },
    }, h('span', { class: 'swatch', style: `background:${color}` }), h('div', { class: 'grow' }, h('div', { class: 'pt' }, title), sub ? h('div', { class: 'ps' }, sub) : null), extra ? h('span', { class: 'badge' }, extra) : null);
    root.innerHTML = '';

    const rooms = f.rooms.filter((r) => match(r.name, r.type, r.notes, r.floorFinish));
    const total = f.rooms.reduce((a, r) => a + polygonArea(r.points), 0);
    root.append(h('h4', {}, `Rooms & places (${f.rooms.length}) · ${fmtArea(total, units)}`));
    if (!rooms.length) root.append(h('div', { class: 'empty' }, 'No rooms yet — use ⬚ Room or ▢ Area to mark places.'));
    for (const r of rooms) {
      const items = f.items.filter((it) => pointInRoom(it, r)).length;
      root.append(row('room', r.id, r.name, [r.type !== 'Other' && r.type !== r.name ? r.type : '', r.notes].filter(Boolean).join(' · '), r.color, `${fmtArea(polygonArea(r.points), units)}${items ? ' · ' + items + ' items' : ''}`));
    }

    const pins = f.notes.filter((n) => n.kind === 'pin');
    root.append(h('h4', {}, `Notes (${pins.length})`));
    if (!pins.length) root.append(h('div', { class: 'empty' }, 'No notes — use 📍 to annotate places (e.g. “add socket here”).'));
    pins.forEach((n, i) => { if (match(n.title, n.text)) root.append(row('note', n.id, `${i + 1}. ${n.title || 'Note'}`, n.text, n.color)); });

    const labels = f.notes.filter((n) => n.kind !== 'pin' && match(n.text, n.kind));
    if (labels.length) {
      root.append(h('h4', {}, 'Labels & dimensions'));
      for (const n of labels) root.append(row('note', n.id, n.kind === 'dim' ? `↔ ${n.text || 'Dimension'}` : n.text || 'Label', n.kind === 'dim' ? fmtLen(Math.hypot(n.x2 - n.x, n.y2 - n.y), units) : '', n.color));
    }

    const annotated = f.items.filter((it) => (it.label || it.notes) && match(it.label, it.notes, ITEMS[it.catalogId]?.name));
    const annotatedOps = f.openings.filter((o) => (o.label || o.notes) && match(o.label, o.notes));
    if (annotated.length || annotatedOps.length) {
      root.append(h('h4', {}, 'Annotated items'));
      for (const it of annotated) root.append(row('item', it.id, it.label || ITEMS[it.catalogId]?.name, it.notes, it.color || ITEMS[it.catalogId]?.color));
      for (const o of annotatedOps) root.append(row('opening', o.id, o.label || ITEMS[o.catalogId]?.name, o.notes, ITEMS[o.catalogId]?.color));
    }
    if (q) {
      const found = f.items.filter((it) => !annotated.includes(it) && match(ITEMS[it.catalogId]?.name, ITEMS[it.catalogId]?.cat));
      if (found.length) {
        root.append(h('h4', {}, 'Matching items'));
        for (const it of found) root.append(row('item', it.id, ITEMS[it.catalogId]?.name, ITEMS[it.catalogId]?.cat, it.color || ITEMS[it.catalogId]?.color));
      }
    }
  }

  focusObj(type, id) {
    const o = this.coll(type).find((x) => x.id === id);
    if (!o) return;
    let x, y;
    if (type === 'room') [x, y] = polygonCentroid(o.points);
    else if (type === 'wall') { x = (o.x1 + o.x2) / 2; y = (o.y1 + o.y2) / 2; }
    else if (type === 'opening') { const w = this.floor.walls.find((w) => w.id === o.wallId); if (!w) return; const g = openingGeom(o, w); x = g.cx; y = g.cy; }
    else { x = o.x; y = o.y; }
    const [sx, sy] = this.plan.w2s(x, y);
    if (sx < 40 || sy < 40 || sx > this.plan.W - 40 || sy > this.plan.H - 40) this.plan.focusOn(x, y);
  }

  scheduleRows(allFloors) {
    const floors = allFloors ? this.design.floors : [this.floor];
    const rows = new Map();
    const add = (cat, name, w, d, hh, floor, room, label) => {
      const key = cat + '|' + name + '|' + w + '|' + d + '|' + hh;
      if (!rows.has(key)) rows.set(key, { cat, name, w, d, h: hh, qty: 0, where: new Set(), labels: new Set() });
      const r = rows.get(key);
      r.qty++;
      r.where.add(allFloors ? `${floor}${room ? ': ' + room : ''}` : room || '—');
      if (label) r.labels.add(label);
    };
    for (const f of floors) {
      for (const it of f.items) {
        const def = ITEMS[it.catalogId] || { cat: 'Other', name: it.catalogId };
        const room = f.rooms.find((r) => pointInRoom(it, r));
        add(def.cat, def.name, Math.round(it.w * 100), Math.round(it.d * 100), Math.round(it.h * 100), f.name, room?.name, it.label);
      }
      for (const o of f.openings) {
        const def = ITEMS[o.catalogId] || { cat: 'Doors & Windows', name: o.catalogId };
        add(def.cat, def.name, Math.round(o.width * 100), '', Math.round(o.height * 100), f.name, '', o.label);
      }
    }
    return [...rows.values()].sort((a, b) => a.cat.localeCompare(b.cat) || a.name.localeCompare(b.name));
  }

  renderSchedule() {
    const root = $('#schedule');
    if (!root || !$('#pane-schedule').classList.contains('active')) return;
    const all = $('#sched-all').checked;
    const rows = this.scheduleRows(all);
    root.innerHTML = '';
    const floors = all ? this.design.floors : [this.floor];
    const walls = floors.reduce((a, f) => a + f.walls.reduce((s, w) => s + wallLength(w), 0), 0);
    const area = floors.reduce((a, f) => a + f.rooms.reduce((s, r) => s + polygonArea(r.points), 0), 0);
    const units = this.design.settings.units;
    root.append(h('div', { class: 'stats' },
      h('div', {}, h('b', {}, rows.reduce((a, r) => a + r.qty, 0)), ' items & openings'),
      h('div', {}, h('b', {}, fmtArea(area, units)), ' room area'),
      h('div', {}, h('b', {}, fmtLen(walls, units)), ' of walls'),
      h('div', {}, h('b', {}, floors.reduce((a, f) => a + f.rooms.length, 0)), ' rooms')));
    let cat = null, tbody;
    for (const r of rows) {
      if (r.cat !== cat) {
        cat = r.cat;
        tbody = h('tbody');
        root.append(h('h4', {}, cat), h('table', { class: 'sched' }, tbody));
      }
      tbody.append(h('tr', { title: [...r.where].join('\n') }, h('td', { class: 'q' }, r.qty + '×'), h('td', {}, r.name, r.labels.size ? h('div', { class: 'ps' }, [...r.labels].join(', ')) : null), h('td', { class: 'dim' }, [r.w, r.d, r.h].filter((x) => x !== '').join('×'))));
    }
    if (!rows.length) root.append(h('div', { class: 'empty' }, 'Nothing placed yet. Pick items from the Catalog tab.'));
  }

  scheduleCSV() {
    const rows = this.scheduleRows(true);
    const q = (s) => `"${String(s).replace(/"/g, '""')}"`;
    const lines = [['Category', 'Item', 'Qty', 'Width cm', 'Depth cm', 'Height cm', 'Where', 'Labels'].join(',')];
    for (const r of rows) lines.push([r.cat, r.name, r.qty, r.w, r.d, r.h, [...r.where].join('; '), [...r.labels].join('; ')].map(q).join(','));
    lines.push('');
    lines.push(['Floor', 'Room', 'Type', 'Area m2', 'Perimeter m', 'Finish', 'Notes'].map(q).join(','));
    for (const f of this.design.floors) for (const r of f.rooms) lines.push([f.name, r.name, r.type, polygonArea(r.points).toFixed(2), polygonPerimeter(r.points).toFixed(2), r.floorFinish || '', r.notes || ''].map(q).join(','));
    lines.push('');
    lines.push(['Floor', 'No', 'Note', 'Text', 'x m', 'y m'].map(q).join(','));
    for (const f of this.design.floors) f.notes.filter((n) => n.kind === 'pin').forEach((n, i) => lines.push([f.name, i + 1, n.title || '', n.text || '', n.x.toFixed(2), n.y.toFixed(2)].map(q).join(',')));
    return lines.join('\n');
  }

  printSheet() {
    const units = this.design.settings.units;
    const parts = [];
    const cur = this.floorIdx;
    for (let i = 0; i < this.design.floors.length; i++) {
      this.floorIdx = i;
      const f = this.floor;
      const img = this.plan.exportPNG(1.5);
      const rooms = f.rooms.map((r) => `<tr><td>${esc(r.name)}</td><td>${esc(r.type)}</td><td>${fmtArea(polygonArea(r.points), units)}</td><td>${esc(r.floorFinish || '')}</td><td>${esc(r.notes || '')}</td></tr>`).join('');
      const pins = f.notes.filter((n) => n.kind === 'pin').map((n, k) => `<li><b>${k + 1}. ${esc(n.title || 'Note')}</b> ${esc(n.text || '')}</li>`).join('');
      parts.push(`<section><h2>${esc(f.name)}</h2><img src="${img}"><h3>Rooms</h3><table><tr><th>Room</th><th>Type</th><th>Area</th><th>Finish</th><th>Notes</th></tr>${rooms}</table>${pins ? `<h3>Notes</h3><ol>${pins}</ol>` : ''}</section>`);
    }
    this.floorIdx = cur;
    this.plan.render();
    const rows = this.scheduleRows(true).map((r) => `<tr><td>${esc(r.cat)}</td><td>${esc(r.name)}</td><td>${r.qty}</td><td>${[r.w, r.d, r.h].filter((x) => x !== '').join('×')} cm</td><td>${esc([...r.where].join('; '))}</td></tr>`).join('');
    const html = `<!doctype html><meta charset="utf-8"><title>${esc(this.design.name)}</title>
<style>body{font:13px system-ui,sans-serif;margin:24px;color:#222}img{max-width:100%;border:1px solid #ddd}table{border-collapse:collapse;width:100%;margin:8px 0 16px}td,th{border:1px solid #ccc;padding:4px 6px;text-align:left}section{page-break-after:always}h1{margin:0 0 4px}</style>
<h1>${esc(this.design.name)}</h1><div>${esc(this.design.description || '')}</div>${parts.join('')}<h2>Item schedule</h2><table><tr><th>Category</th><th>Item</th><th>Qty</th><th>Size</th><th>Where</th></tr>${rows}</table>`;
    const wnd = window.open('', '_blank');
    if (!wnd) { download(safeName(this.design.name) + '_sheet.html', html, 'text/html'); return; }
    wnd.document.write(html); wnd.document.close();
    setTimeout(() => wnd.print(), 400);
  }

  // ---------- properties panel ----------
  refreshProps() {
    const root = $('#props');
    root.innerHTML = '';
    const sel = this.selection, obj = this.getSelected();
    if (!obj) { this.designProps(root); return; }
    const fn = { wall: this.wallProps, opening: this.openingProps, item: this.itemProps, room: this.roomProps, note: this.noteProps }[sel.type];
    fn.call(this, root, obj);
  }

  // update just the numbers while dragging (keeps focus/caret intact)
  refreshPropsValues() {
    document.querySelectorAll('#props [data-live]').forEach((el) => { if (document.activeElement !== el) el._sync?.(); });
  }

  // field builders
  F = {
    row: (label, ...inputs) => h('label', { class: 'row' }, h('span', { class: 'lbl' }, label), h('span', { class: 'inp' }, ...inputs)),
    num: (obj, key, { k = 100, step = 1, min = -1e9, max = 1e9, unit = 'cm', after, label } = {}) => {
      const el = h('input', { type: 'number', step, 'data-live': 1 });
      el._sync = () => { el.value = +(obj[key] * k).toFixed(2); };
      el._sync();
      el.onchange = () => {
        const v = parseFloat(el.value);
        if (!isFinite(v)) { el._sync(); return; }
        obj[key] = Math.max(min, Math.min(max, v)) / k;
        after?.(obj[key]);
        this.commit('Edit ' + (label || key));
        el._sync();
      };
      return h('span', { class: 'numwrap' }, el, h('span', { class: 'unit' }, unit));
    },
    text: (obj, key, { placeholder = '', label } = {}) => {
      const el = h('input', { type: 'text', value: obj[key] || '', placeholder });
      el.onchange = () => { obj[key] = el.value; this.commit('Edit ' + (label || key)); };
      return el;
    },
    area: (obj, key, { placeholder = '' } = {}) => {
      const el = h('textarea', { rows: 3, placeholder });
      el.value = obj[key] || '';
      el.onchange = () => { obj[key] = el.value; this.commit('Edit notes'); };
      return el;
    },
    color: (obj, key, fallback, after) => {
      const el = h('input', { type: 'color', value: toHex(obj[key] || fallback) });
      el.oninput = () => { obj[key] = el.value; after?.(); this.liveUpdate(); };
      el.onchange = () => { obj[key] = el.value; after?.(); this.commit('Colour'); };
      return el;
    },
    check: (obj, key, text) => {
      const el = h('input', { type: 'checkbox', checked: !!obj[key] });
      el.onchange = () => { obj[key] = el.checked; this.commit(text || key); };
      return h('span', { class: 'chk' }, el, ' ', text || '');
    },
    select: (obj, key, options, after) => {
      const el = h('select', {}, options.map(([v, t]) => h('option', { value: v, selected: String(obj[key]) === String(v) }, t)));
      el.onchange = () => { obj[key] = el.value; after?.(el.value); this.commit('Edit ' + key, { props: true }); };
      return el;
    },
    btn: (text, onclick, cls = '') => h('button', { class: 'btn ' + cls, onclick }, text),
    info: (label, value) => h('div', { class: 'row info' }, h('span', { class: 'lbl' }, label), h('span', { class: 'inp' }, value)),
  };

  header(root, title, sub, color) {
    root.append(h('div', { class: 'prop-head' }, color ? h('span', { class: 'swatch big', style: `background:${color}` }) : null, h('div', {}, h('div', { class: 'pt' }, title), sub ? h('div', { class: 'ps' }, sub) : null)));
  }

  designProps(root) {
    const F = this.F, d = this.design, s = d.settings, f = this.floor;
    this.header(root, 'Design', 'Nothing selected — design & floor settings');
    root.append(
      F.row('Name', F.text(d, 'name', { label: 'name' })),
      F.row('Author', F.text(d, 'author', { placeholder: 'who drew this' })),
      F.row('Description', F.area(d, 'description', { placeholder: 'address, plot size, brief…' })),
      h('h4', {}, 'Floor: ' + f.name),
      F.row('Floor name', F.text(f, 'name', { label: 'floor name' })),
      F.row('Elevation', F.num(f, 'elevation', { min: -100 })),
      F.row('Ceiling height', F.num(f, 'height', { min: 100, max: 1000 })),
      h('div', { class: 'btns' },
        F.btn('Apply height to all walls', () => { f.walls.forEach((w) => (w.height = f.height)); this.commit('Wall heights'); }),
        F.btn('Delete floor', () => this.deleteFloor(), 'danger')),
      ...this.underlayProps(f),
      h('h4', {}, 'Drawing'),
      F.row('Grid', F.select(s, 'grid', [[0.05, '5 cm'], [0.1, '10 cm'], [0.25, '25 cm'], [0.5, '50 cm'], [1, '1 m']], (v) => (s.grid = +v))),
      F.row('', F.check(s, 'snap', 'Snap to grid')),
      F.row('', F.check(s, 'showDims', 'Show all wall lengths')),
      F.row('', F.check({ get v() { return s.showLabels !== false; }, set v(x) { s.showLabels = x; } }, 'v', 'Show item names')),
      h('h4', {}, 'Wall types'),
      F.row('Exterior walls', F.num(s, 'extThickness', { min: 2, max: 100, after: () => this.applyThickness('exterior'), label: 'exterior thickness' })),
      F.row('Interior walls', F.num(s, 'intThickness', { min: 2, max: 100, after: () => this.applyThickness('interior'), label: 'interior thickness' })),
      F.row('New walls are', F.select(s, 'wallKind', [['exterior', 'Exterior'], ['interior', 'Interior']], (v) => this.setWallKind(v))),
      F.row('New wall height', F.num(s, 'wallHeight', { min: 100, max: 1000 })),
      h('div', { class: 'btns' },
        F.btn('9″ / 4.5″ (brick)', () => { s.extThickness = 0.23; s.intThickness = 0.115; this.applyThickness('exterior'); this.applyThickness('interior'); this.commit('Wall thickness', { props: true }); }),
        F.btn('Auto-detect ext / int', () => this.autoWallKinds())),
      h('h4', {}, 'This floor'),
      F.info('Walls', `${f.walls.length} · ${fmtLen(f.walls.reduce((a, w) => a + wallLength(w), 0), s.units)}`),
      F.info('Rooms', `${f.rooms.length} · ${fmtArea(f.rooms.reduce((a, r) => a + polygonArea(r.points), 0), s.units)}`),
      F.info('Items', `${f.items.length} items · ${f.openings.length} doors/windows`),
      F.info('Notes', `${f.notes.length}`),
      h('p', { class: 'tip' }, 'Tip: pick an item in the Catalog, then click on the plan. Press R to rotate while placing. Doors and windows snap into walls.'),
    );
  }

  underlayProps(f) {
    const F = this.F, u = f.underlay;
    const out = [h('h4', {}, 'Trace image (sketch / photo / plan)')];
    const fileBtn = F.btn(u ? 'Replace image…' : 'Upload image…', () => $('#file-underlay').click());
    out.push(h('div', { class: 'btns' }, fileBtn, F.btn('From web URL…', () => this.loadUnderlayURL())));
    if (!u) { out.push(h('p', { class: 'tip' }, 'Put a floor-plan image, a hand sketch or a photo under this floor, set its scale, then trace walls over it.')); return out; }
    const op = h('input', { type: 'range', min: 0.05, max: 1, step: 0.05, value: u.opacity ?? 0.5 });
    op.oninput = () => { u.opacity = +op.value; this.plan.requestRender(); };
    op.onchange = () => this.commit('Image opacity');
    out.push(
      F.info('Image', u.name || 'image'),
      F.row('Opacity', op),
      F.row('Scale', F.num(u, 'mpp', { k: 1000, step: 0.1, unit: 'mm/px', label: 'image scale' })),
      F.row('', F.check(u, 'visible', 'Show image')),
      h('div', { class: 'btns' },
        F.btn('📏 Calibrate scale', () => this.setTool('calib')),
        F.btn('✋ Move image', () => this.setTool('movebg')),
        F.btn('Remove', () => { const src = u.src; delete f.underlay; if (src?.startsWith('asset:') && !this.design.floors.some((x) => x.underlay?.src === src)) delete this.design.assets[src.slice(6)]; this.commit('Remove trace image', { props: true }); }, 'danger')),
    );
    return out;
  }

  wallProps(root, w) {
    const F = this.F, units = this.design.settings.units;
    this.header(root, 'Wall', fmtLen(wallLength(w), units));
    const app = this;
    const len = { get v() { return wallLength(w); }, set v(L) { app.setWallLength(w, L); } };
    root.append(
      F.row('Length', F.num(len, 'v', { min: 5, label: 'length' })),
      h('p', { class: 'tip' }, 'Fast: with the wall selected just type the new length and press Enter (e.g. 3.6 or 360cm). Joined walls follow.'),
      F.row('Type', F.select(w, 'kind', [['exterior', 'Exterior'], ['interior', 'Interior']], (v) => { w.thickness = v === 'interior' ? this.design.settings.intThickness : this.design.settings.extThickness; })),
      F.row('Thickness', F.num(w, 'thickness', { min: 2, max: 100 })),
      F.row('Height', F.num(w, 'height', { min: 10, max: 2000 })),
      F.row('Colour (3D)', F.color(w, 'color', '#eeeeee')),
      F.row('Label', F.text(w, 'label', { placeholder: 'e.g. load-bearing' })),
      F.row('Notes', F.area(w, 'notes')),
      F.info('Openings', String(this.floor.openings.filter((o) => o.wallId === w.id).length)),
      h('div', { class: 'btns' },
        F.btn('Split in half', () => this.splitWall(w)),
        F.btn('Delete', () => this.deleteSelection(), 'danger')),
    );
  }

  openingProps(root, o) {
    const F = this.F, def = ITEMS[o.catalogId] || {};
    const w = this.floor.walls.find((x) => x.id === o.wallId);
    this.header(root, def.name || 'Opening', 'Door / window — drag to slide along the wall', def.color);
    const types = CATEGORIES.find((c) => c.name === 'Doors & Windows').items.map((i) => [i.id, i.name]);
    const offset = { get v() { return o.t * wallLength(w); }, set v(x) { o.t = x / wallLength(w); o.t = clampOpeningT(o, w); } };
    root.append(
      F.row('Type', F.select(o, 'catalogId', types, (v) => { const nd = ITEMS[v]; o.width = nd.w; o.height = nd.h; o.sill = nd.sill ?? 0; o.t = clampOpeningT(o, w); })),
      F.row('Width', F.num(o, 'width', { min: 20, max: 1000, after: () => (o.t = clampOpeningT(o, w)) })),
      F.row('Height', F.num(o, 'height', { min: 10, max: 1000 })),
      def.model === 'window' ? F.row('Sill height', F.num(o, 'sill', { min: 0, max: 500 })) : null,
      F.row('From wall start', F.num(offset, 'v', { label: 'position' })),
      F.row('Label', F.text(o, 'label', { placeholder: 'e.g. D1, W3' })),
      F.row('Notes', F.area(o, 'notes', { placeholder: 'material, hardware, grill…' })),
      ['door', 'door2'].includes(def.model) ? h('div', { class: 'btns' },
        F.btn('Flip swing side', () => { o.swing = -(o.swing || 1); this.commit('Flip door'); }),
        def.model === 'door' ? F.btn('Flip hinge', () => { o.hinge = o.hinge ? 0 : 1; this.commit('Flip hinge'); }) : null) : null,
      h('div', { class: 'btns' }, F.btn('Delete', () => this.deleteSelection(), 'danger')),
    );
  }

  itemProps(root, it) {
    const F = this.F, def = ITEMS[it.catalogId] || { name: it.catalogId, cat: '' };
    this.header(root, it.label || def.name, `${def.cat} · ${def.name}`, it.color || def.color);
    const pos = (k) => F.num(it, k, { k: 1, step: 0.05, unit: 'm' });
    root.append(
      F.row('Label', F.text(it, 'label', { placeholder: def.name })),
      F.row('Width', F.num(it, 'w', { min: 1 })),
      F.row('Depth', F.num(it, 'd', { min: 1 })),
      F.row('Height', F.num(it, 'h', { min: 0.5 })),
      F.row('Above floor', F.num(it, 'elev', { min: -500 })),
      F.row('Rotation', F.num(it, 'rot', { k: 1, step: 15, unit: '°', after: () => (it.rot = ((it.rot % 360) + 360) % 360) })),
      F.row('X', pos('x')), F.row('Y', pos('y')),
      F.row('Colour', F.color(it, 'color', def.color)),
      F.row('Notes', F.area(it, 'notes', { placeholder: 'brand, price, link, status…' })),
      F.row('', F.check(it, 'locked', 'Lock (can’t be moved by dragging)')),
      ...this.modelProps(it, def),
      h('div', { class: 'btns' },
        F.btn('⟳ 90°', () => { it.rot = (it.rot + 90) % 360; this.commit('Rotate'); }),
        F.btn('Duplicate', () => this.duplicate()),
        F.btn('Reset size', () => { it.w = def.w; it.d = def.d; it.h = def.h; it.elev = def.elev || 0; it.color = def.color; this.commit('Reset', { props: true }); }),
        F.btn('Delete', () => this.deleteSelection(), 'danger')),
    );
  }

  modelProps(it, def) {
    const F = this.F;
    const shared = this.design.modelMap?.[it.catalogId];
    const url = it.model3d || shared;
    const out = [h('h4', {}, '3D model')];
    out.push(F.info('Using', url ? (it.modelName || url.split('/').pop()) + (it.model3d ? '' : ' (for all ' + def.name + ')') : 'built-in shape'));
    out.push(h('div', { class: 'btns' },
      F.btn('🔎 Browse free models', () => this.openModels(it)),
      F.btn('⬆ Upload .glb', () => { this._modelTarget = it; $('#file-model').click(); })));
    if (url) {
      out.push(h('div', { class: 'btns' },
        F.btn('⟳ Turn model 90°', () => { it.modelRot = ((it.modelRot || 0) + 90) % 360; this.commit('Turn model'); }),
        F.btn('Use model’s real size', () => this.modelRealSize(it, url)),
        it.model3d ? F.btn(`Use for every “${def.name}”`, () => { this.design.modelMap = this.design.modelMap || {}; this.design.modelMap[it.catalogId] = it.model3d; this.commit('Model for all', { props: true }); this.toast('All ' + def.name + ' items now use this model'); }) : null,
        F.btn('Remove model', () => { if (it.model3d) { delete it.model3d; delete it.modelName; } else delete this.design.modelMap[it.catalogId]; this.commit('Remove model', { props: true }); }, 'danger')));
    }
    return out;
  }

  async assignModel(it, url, name) {
    it.model3d = url;
    it.modelName = name;
    it.modelRot = 0;
    this.commit('3D model', { props: true });
    await this.ensure3D();
    if (this.view === '2d') this.setView('split');
    this.modelRealSize(it, url, true);
  }

  // copy the model's true size (metres) to the item — Poly Haven models are real-world scale
  modelRealSize(it, url, quiet = false) {
    if (!this.view3d) { this.ensure3D().then(() => this.modelRealSize(it, url, quiet)); return; }
    this.view3d.getModel(url, (m) => {
      const s = this.view3d.modelSize(m, it.modelRot || 0);
      if (!(s.x > 0.02 && s.y > 0.02 && s.z > 0.02 && s.x < 50 && s.z < 50)) { if (!quiet) this.toast('This model has no usable real size'); return; }
      it.w = +s.x.toFixed(3); it.h = +s.y.toFixed(3); it.d = +s.z.toFixed(3);
      this.commit('Model size', { props: true });
    });
  }

  async uploadModel(file) {
    const it = this._modelTarget;
    if (!file || !it) return;
    if (!/\.glb$/i.test(file.name)) { alert('Please choose a .glb file (binary glTF). Most 3D sites and Blender can export it.'); return; }
    try {
      if (this.store.server) {
        const r = await fetch('api/models/upload/' + encodeURIComponent(safeName(file.name.replace(/\.glb$/i, '')) + '.glb'), { method: 'PUT', body: file });
        if (!r.ok) throw new Error(await r.text());
        const { url } = await r.json();
        await this.assignModel(it, url, file.name);
      } else {
        const data = await new Promise((res, rej) => { const fr = new FileReader(); fr.onload = () => res(fr.result); fr.onerror = rej; fr.readAsDataURL(file); });
        const id = uid('glb');
        this.design.assets = this.design.assets || {};
        this.design.assets[id] = data;
        await this.assignModel(it, 'asset:' + id, file.name);
      }
      this.toast('Model added: ' + file.name);
    } catch (err) { alert('Upload failed: ' + err.message); }
  }

  async openModels(it) {
    this._modelTarget = it;
    const dlg = $('#models-dlg');
    dlg.showModal();
    if (!this.store.server) { $('#models-grid').innerHTML = '<div class="empty">The free model library needs the local server (run start.sh). You can still upload a .glb file.</div>'; return; }
    if (!this._phModels) {
      $('#models-grid').innerHTML = '<div class="empty">Loading the Poly Haven catalogue…</div>';
      try {
        const r = await fetch('api/polyhaven/models');
        if (!r.ok) throw new Error(await r.text());
        this._phModels = (await r.json()).models;
      } catch (err) { $('#models-grid').innerHTML = `<div class="empty">${esc(err.message)} — check your internet connection.</div>`; return; }
      const cats = [...new Set(this._phModels.flatMap((m) => m.categories.filter((c) => !c.startsWith('collection'))))].sort();
      $('#models-cat').innerHTML = '<option value="">All categories</option>' + cats.map((c) => `<option>${esc(c)}</option>`).join('');
      $('#models-search').oninput = () => this.renderModels();
      $('#models-cat').onchange = () => this.renderModels();
    }
    // start with a search that fits the item
    const def = ITEMS[it.catalogId];
    $('#models-search').value = guessModelQuery(def);
    this.renderModels();
  }

  renderModels() {
    const q = $('#models-search').value.trim().toLowerCase(), cat = $('#models-cat').value;
    const words = q.split(/\s+/).filter(Boolean);
    let list = this._phModels.filter((m) => (!cat || m.categories.includes(cat)) && (!words.length || words.some((w) => (m.name + ' ' + m.id + ' ' + m.categories.join(' ') + ' ' + m.tags.join(' ')).toLowerCase().includes(w))));
    if (!list.length && words.length) list = this._phModels.filter((m) => !cat || m.categories.includes(cat));
    const grid = $('#models-grid');
    grid.innerHTML = '';
    for (const m of list.slice(0, 120)) {
      const card = h('div', { class: 'model-card', title: m.categories.join(', ') },
        h('img', { loading: 'lazy', src: `https://cdn.polyhaven.com/asset_img/thumbs/${m.id}.png?width=256&height=256`, alt: m.name }),
        h('div', {}, m.name));
      card.onclick = async () => {
        card.classList.add('busy');
        this.toast('Downloading ' + m.name + '…');
        try {
          const r = await fetch('api/polyhaven/download?id=' + encodeURIComponent(m.id));
          if (!r.ok) throw new Error(await r.text());
          const { url } = await r.json();
          $('#models-dlg').close();
          await this.assignModel(this._modelTarget, url, m.name);
          this.toast(m.name + ' added (CC0, Poly Haven)');
        } catch (err) { alert('Download failed: ' + err.message); }
        card.classList.remove('busy');
      };
      grid.append(card);
    }
    if (!list.length) grid.append(h('div', { class: 'empty' }, 'No models match.'));
  }

  roomProps(root, r) {
    const F = this.F, units = this.design.settings.units;
    this.header(root, r.name, `${fmtArea(polygonArea(r.points), units)} · perimeter ${fmtLen(polygonPerimeter(r.points), units)}`, r.color);
    const items = this.floor.items.filter((it) => pointInRoom(it, r));
    root.append(
      F.row('Name', F.text(r, 'name', { label: 'room name' })),
      F.row('Type', F.select(r, 'type', ROOM_TYPES.map(([n]) => [n, n]), (v) => {
        const t = ROOM_TYPES.find(([n]) => n === v);
        if (t) r.color = t[1];
        if (/^Room \d+$/.test(r.name) || ROOM_TYPES.some(([n]) => n === r.name)) r.name = v;
      })),
      F.row('Floor colour', F.color(r, 'color', '#f5f5f5')),
      F.row('Floor finish', F.text(r, 'floorFinish', { placeholder: 'e.g. vitrified tiles 600×600' })),
      F.row('3D floor', F.select(r, 'material', FLOOR_MATERIALS)),
      F.row('Wall paint', F.text(r, 'wallFinish', { placeholder: 'e.g. Asian Paints “Ivory”' })),
      F.row('Notes', F.area(r, 'notes', { placeholder: 'use, wishes, measurements, to-dos…' })),
      F.row('', F.check(r, 'hideLabel', 'Hide name on plan')),
      F.info('Area', fmtArea(polygonArea(r.points), units)),
      F.info('Corners', String(r.points.length)),
      F.info('Items inside', items.length ? items.map((it) => it.label || ITEMS[it.catalogId]?.name).slice(0, 12).join(', ') + (items.length > 12 ? '…' : '') : '—'),
      h('div', { class: 'btns' },
        F.btn('Build walls around', () => this.wallsAroundRoom(r)),
        F.btn('Delete room', () => this.deleteSelection(), 'danger')),
      h('p', { class: 'tip' }, 'Drag the corner squares to reshape. Rooms are the “places” of your plan — name them and add notes.'),
    );
  }

  noteProps(root, n) {
    const F = this.F;
    if (n.kind === 'pin') {
      const no = this.floor.notes.filter((x) => x.kind === 'pin').indexOf(n) + 1;
      this.header(root, `Note ${no}`, 'Annotation pin', n.color);
      root.append(
        F.row('Title', F.text(n, 'title', { label: 'note title' })),
        F.row('Text', F.area(n, 'text', { placeholder: 'What should happen here?' })),
        F.row('Colour', F.color(n, 'color', '#e53935')),
      );
    } else if (n.kind === 'label') {
      this.header(root, 'Text label', '', n.color);
      root.append(
        F.row('Text', F.text(n, 'text', { label: 'label' })),
        F.row('Font size', F.num(n, 'size', { k: 1, min: 8, max: 72, unit: 'px' })),
        F.row('Colour', F.color(n, 'color', '#263238')),
      );
    } else {
      this.header(root, 'Dimension', fmtLen(Math.hypot(n.x2 - n.x, n.y2 - n.y), this.design.settings.units), n.color);
      const len = { get v() { return Math.hypot(n.x2 - n.x, n.y2 - n.y); }, set v(L) { const a = Math.atan2(n.y2 - n.y, n.x2 - n.x); n.x2 = n.x + Math.cos(a) * L; n.y2 = n.y + Math.sin(a) * L; } };
      root.append(
        F.row('Caption', F.text(n, 'text', { placeholder: 'optional' })),
        F.row('Length', F.num(len, 'v', { min: 1, label: 'length' })),
        F.row('Colour', F.color(n, 'color', '#d81b60')),
      );
    }
    root.append(h('div', { class: 'btns' }, F.btn('Delete', () => this.deleteSelection(), 'danger')));
  }

  // ---------- editing actions ----------
  deleteSelection() {
    const sel = this.selection;
    if (!sel) return;
    const c = this.coll(sel.type);
    const i = c.findIndex((o) => o.id === sel.id);
    if (i < 0) return;
    c.splice(i, 1);
    if (sel.type === 'wall') this.floor.openings = this.floor.openings.filter((o) => o.wallId !== sel.id);
    this.selection = null;
    this.commit('Delete', { props: true });
  }

  duplicate() {
    const sel = this.selection, obj = this.getSelected();
    if (!obj) return;
    const copy = JSON.parse(JSON.stringify(obj));
    copy.id = uid(sel.type[0]);
    const off = 0.3;
    if (sel.type === 'item' || sel.type === 'note') { copy.x += off; copy.y += off; if (copy.x2 !== undefined) { copy.x2 += off; copy.y2 += off; } }
    else if (sel.type === 'wall') { copy.x1 += off; copy.y1 += off; copy.x2 += off; copy.y2 += off; }
    else if (sel.type === 'room') { copy.points = copy.points.map(([x, y]) => [x + off, y + off]); copy.name += ' (copy)'; }
    else if (sel.type === 'opening') { const w = this.floor.walls.find((x) => x.id === copy.wallId); copy.t = clampOpeningT({ ...copy, t: copy.t + (copy.width + 0.1) / wallLength(w) }, w); }
    this.coll(sel.type).push(copy);
    this.selection = { type: sel.type, id: copy.id };
    this.commit('Duplicate', { props: true });
  }

  rotateSelection(deg) {
    const obj = this.getSelected();
    if (this.selection?.type === 'item' && obj) { obj.rot = (((obj.rot + deg) % 360) + 360) % 360; this.commit('Rotate'); }
    else if (this.selection?.type === 'opening' && obj) { obj.swing = -(obj.swing || 1); this.commit('Flip door'); }
  }

  nudge(dx, dy) {
    const sel = this.selection, o = this.getSelected();
    if (!o) return;
    if (sel.type === 'item' || sel.type === 'note') { o.x += dx; o.y += dy; if (o.x2 !== undefined) { o.x2 += dx; o.y2 += dy; } }
    else if (sel.type === 'wall') { o.x1 += dx; o.y1 += dy; o.x2 += dx; o.y2 += dy; }
    else if (sel.type === 'room') o.points = o.points.map(([x, y]) => [x + dx, y + dy]);
    else return;
    this.commit('Nudge');
  }

  splitWall(w) {
    const mx = (w.x1 + w.x2) / 2, my = (w.y1 + w.y2) / 2;
    const L = wallLength(w);
    const w2 = { ...w, id: uid('w'), x1: mx, y1: my };
    for (const o of this.floor.openings.filter((o) => o.wallId === w.id)) {
      const s = o.t * L;
      if (s > L / 2) { o.wallId = w2.id; o.t = (s - L / 2) / (L / 2); } else o.t = s / (L / 2);
    }
    w.x2 = mx; w.y2 = my;
    this.floor.walls.push(w2);
    this.commit('Split wall', { props: true });
  }

  wallsAroundRoom(r) {
    const s = this.design.settings;
    let n = 0;
    for (let i = 0; i < r.points.length; i++) {
      const [ax, ay] = r.points[i], [bx, by] = r.points[(i + 1) % r.points.length];
      const exists = this.floor.walls.some((w) => (Math.hypot(w.x1 - ax, w.y1 - ay) < 1e-3 && Math.hypot(w.x2 - bx, w.y2 - by) < 1e-3) || (Math.hypot(w.x1 - bx, w.y1 - by) < 1e-3 && Math.hypot(w.x2 - ax, w.y2 - ay) < 1e-3));
      if (!exists) { this.floor.walls.push({ id: uid('w'), x1: ax, y1: ay, x2: bx, y2: by, ...this.newWallProps() }); n++; }
    }
    this.commit(`Added ${n} walls`);
    this.toast(n ? `Added ${n} walls` : 'This room already has walls on every side');
  }

  // ---------- files ----------
  async newDesign() {
    if (this.dirty && !confirm('Start a new design? Unsaved changes to the current one will be lost (it stays in the undo draft only until you start).')) return;
    const name = prompt('Name for the new design:', 'My home');
    if (name === null) return;
    this.loadDesign(newDesign(name || 'Untitled design'), null);
    this.saveDraftNow();
  }

  async save(asNew = false) {
    let name = this.fileName;
    if (!name || asNew) {
      const n = prompt(`Save design as (saved to ${this.store.where}):`, safeName(this.design.name));
      if (n === null) return;
      name = safeName(n);
      if (asNew || name !== this.fileName) {
        try {
          const list = await this.store.list();
          if (list.some((d) => d.name === name) && !confirm(`“${name}” already exists. Overwrite it?`)) return;
        } catch { /* ignore */ }
      }
    }
    try {
      this.design.updated = new Date().toISOString();
      await this.store.save(name, this.design);
      this.fileName = name;
      this.dirty = false;
      this.updateTitle();
      this.saveDraftNow();
      this.toast(this.store.server ? `Saved designs/${name}.json` : `Saved “${name}” in this browser`);
    } catch (err) {
      alert('Save failed: ' + err.message);
    }
  }

  async openDialog() {
    const dlg = $('#open-dlg');
    const list = $('#open-list');
    list.innerHTML = '<div class="empty">Loading…</div>';
    $('#open-where').textContent = this.store.where;
    dlg.showModal();
    try {
      const designs = await this.store.list();
      list.innerHTML = '';
      if (!designs.length) list.append(h('div', { class: 'empty' }, 'No saved designs yet.'));
      for (const d of designs) {
        list.append(h('div', { class: 'open-row' },
          h('div', { class: 'grow' }, h('div', { class: 'pt' }, d.name), h('div', { class: 'ps' }, (d.title && d.title !== d.name ? d.title + ' · ' : '') + (d.updated ? new Date(d.updated).toLocaleString() : '') + (d.size ? ` · ${(d.size / 1024).toFixed(1)} KB` : ''))),
          h('button', { class: 'btn primary', onclick: async () => {
            if (this.dirty && !confirm('Open another design? Unsaved changes will be lost.')) return;
            try { const data = await this.store.load(d.name); this.loadDesign(data, d.name); this.saveDraftNow(); dlg.close(); this.toast('Opened ' + d.name); } catch (err) { alert(err.message); }
          } }, 'Open'),
          h('button', { class: 'btn danger', onclick: async () => {
            if (!confirm(`Delete “${d.name}”?${this.store.server ? ' (a backup copy is kept in designs/.history/)' : ''}`)) return;
            await this.store.remove(d.name); this.openDialog();
          } }, 'Delete')));
      }
    } catch (err) {
      list.innerHTML = `<div class="empty">${esc(err.message)}</div>`;
    }
  }

  async importFile(file) {
    if (!file) return;
    try {
      const data = JSON.parse(await file.text());
      if (this.dirty && !confirm('Replace the current design with the imported file?')) return;
      this.loadDesign(data, null, { keepDirty: true });
      this.saveDraftNow();
      this.toast('Imported ' + file.name);
    } catch (err) {
      alert('Could not import: ' + err.message);
    }
    $('#file-import').value = '';
  }

  // ---------- keyboard ----------
  bindKeys() {
    window.addEventListener('keydown', (e) => {
      const tag = (e.target.tagName || '').toLowerCase();
      const typing = tag === 'textarea' || tag === 'select' || (tag === 'input' && !['checkbox', 'radio', 'color', 'range', 'button'].includes(e.target.type));
      const mod = e.ctrlKey || e.metaKey;
      if (mod && e.key.toLowerCase() === 's') { e.preventDefault(); this.save(e.shiftKey); return; }
      if (typing) { if (e.key === 'Escape') e.target.blur(); return; }
      if (mod && e.key.toLowerCase() === 'z') { e.preventDefault(); e.shiftKey ? this.redo() : this.undo(); return; }
      if (mod && e.key.toLowerCase() === 'y') { e.preventDefault(); this.redo(); return; }
      if (mod && e.key.toLowerCase() === 'd') { e.preventDefault(); this.duplicate(); return; }
      if (mod) return;
      // typed exact lengths: wall drawing, room size, or resizing the selected wall / dimension
      const lenCtx = (this.tool === 'wall' && this.plan.wallStart) || this.tool === 'roombox' || this.tool === 'area' ||
        (this.tool === 'select' && this.selection && (this.selection.type === 'wall' || (this.selection.type === 'note' && this.getSelected()?.kind === 'dim')));
      if (lenCtx && (/^[0-9.]$/.test(e.key) || (this.lenBuf && /^[xcmfti'"*]$/i.test(e.key)))) {
        this.lenBuf += e.key; this.plan.requestRender(); this.status('Length: ' + this.lenBuf + ' — Enter to apply, Esc to cancel'); e.preventDefault(); return;
      }
      if (this.lenBuf && e.key === 'Backspace') { this.lenBuf = this.lenBuf.slice(0, -1); this.plan.requestRender(); e.preventDefault(); return; }
      if (this.lenBuf && e.key === 'Escape') { this.lenBuf = ''; this.plan.requestRender(); return; }
      if (this.lenBuf && e.key === 'Enter') {
        const buf = this.lenBuf; this.lenBuf = '';
        if (!this.plan.applyTyped(buf)) this.toast(`Could not use “${buf}” here`);
        this.plan.requestRender(); e.preventDefault(); return;
      }
      if (e.key === 'e' || e.key === 'E') { this.setWallKind(this.design.settings.wallKind === 'interior' ? 'exterior' : 'interior'); return; }
      if (e.key === ' ') { this.plan.spaceDown = true; e.preventDefault(); return; }
      if (e.key === 'Escape') { if (this.tool !== 'select') { this.plan.finishTool(); if (this.tool === 'place') this.setTool('select'); } else this.select(null); return; }
      if (e.key === 'Enter') { this.plan.finishTool(); return; }
      if (e.key === 'Delete' || e.key === 'Backspace') { e.preventDefault(); this.deleteSelection(); return; }
      if (e.key === 'r' || e.key === 'R') {
        if (this.tool === 'place') { this.plan.ghostRot = (this.plan.ghostRot + (e.shiftKey ? 15 : 90)) % 360; this.plan.requestRender(); }
        else this.rotateSelection(e.shiftKey ? 15 : 90);
        return;
      }
      if (e.key === 'f' || e.key === 'F') { this.plan.fit(); return; }
      if (e.key === '+' || e.key === '=') { this.plan.zoomBy(1.25); return; }
      if (e.key === '-') { this.plan.zoomBy(0.8); return; }
      const step = e.shiftKey ? 0.1 : 0.01;
      const arrows = { ArrowLeft: [-step, 0], ArrowRight: [step, 0], ArrowUp: [0, -step], ArrowDown: [0, step] };
      if (arrows[e.key] && this.selection) { e.preventDefault(); this.nudge(...arrows[e.key]); return; }
      const t = TOOLS.find((x) => x[3] === e.key.toLowerCase());
      if (t) this.setTool(t[0]);
      if (e.key === '1') this.setView('2d');
      if (e.key === '2') this.setView('split');
      if (e.key === '3') this.setView('3d');
    });
    window.addEventListener('keyup', (e) => { if (e.key === ' ') this.plan.spaceDown = false; });
  }
}

// ---------- helpers ----------
// a sensible first search in the model library for a catalog item
function guessModelQuery(def) {
  if (!def) return '';
  const n = def.name.toLowerCase();
  const map = [[/sofa|loveseat|sectional|couch/, 'sofa'], [/armchair|recliner|reading chair|rocking/, 'armchair'], [/office chair/, 'office chair'], [/chair|stool/, 'chair'],
    [/bed|crib/, 'bed'], [/coffee table|side table/, 'coffee table'], [/table|desk/, 'table'], [/lamp|light|chandelier|pendant|sconce/, 'lamp'],
    [/plant|tree|shrub|tulsi|hedge|bush/, 'plant'], [/shelf|bookcase|rack/, 'shelf'], [/cabinet|wardrobe|almirah|dresser|chest|drawer|sideboard|crockery/, 'cabinet'],
    [/tv|screen|monitor/, 'tv'], [/vase|urli|sculpture|statue/, 'vase'], [/clock/, 'clock'], [/mirror/, 'mirror'], [/rug|mat/, 'rug'], [/bench/, 'bench'], [/fan/, 'fan']];
  for (const [re, q] of map) if (re.test(n)) return q;
  return n.split(/[ (]/)[0];
}

function stripTime(s) { return s.replace(/"updated":"[^"]*"/, ''); }

function pointInRoom(it, r) {
  if (!r.points || r.points.length < 3) return false;
  let inside = false;
  const pts = r.points, x = it.x, y = it.y;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
    const [xi, yi] = pts[i], [xj, yj] = pts[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

function toHex(c) {
  if (/^#[0-9a-f]{6}$/i.test(c)) return c;
  if (/^#[0-9a-f]{3}$/i.test(c)) return '#' + c.slice(1).split('').map((x) => x + x).join('');
  return '#cccccc';
}

// Draw a catalog thumbnail by reusing the plan renderer on a small canvas.
function drawThumb(cv, def) {
  const fake = Object.create(PlanView.prototype);
  const size = 88, pad = 10;
  const isOpening = def.kind === 'opening';
  const ext = Math.max(def.w + (def.model === 'dining' || def.model === 'rdining' ? 0.8 : 0), isOpening ? def.w : def.d + (def.model === 'dining' ? 0.8 : 0), 0.3);
  fake.ctx = cv.getContext('2d');
  fake.dpr = 1;
  fake.scale = (size - pad * 2) / ext;
  fake.ox = size / 2; fake.oy = size / 2;
  fake.app = { design: { settings: {} } };
  const ctx = fake.ctx;
  ctx.clearRect(0, 0, size, size);
  if (isOpening) {
    ctx.setTransform(fake.scale, 0, 0, fake.scale, size / 2 - (def.w / 2 + 0.2) * fake.scale, size / 2 - 0.15 * fake.scale);
    ctx.fillStyle = '#455a64';
    ctx.fillRect(0, -0.075, def.w + 0.4, 0.15);
    fake.drawOpening({ catalogId: def.id, width: def.w, swing: 1, hinge: 0 }, 0.15, 0.2, 0.2 + def.w);
  } else {
    fake.drawItem({ catalogId: def.id, x: 0, y: 0, rot: 0, w: def.w, d: def.d, h: def.h, elev: 0, color: def.color });
  }
  ctx.setTransform(1, 0, 0, 1, 0, 0);
}

const app = new App();
window.homeDesigner = app;
app.start();
