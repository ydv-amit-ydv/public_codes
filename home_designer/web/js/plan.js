// Home Designer — 2D floor-plan editor (canvas).
import {
  uid, dist, projectOnSegment, polygonArea, polygonCentroid, pointInPolygon, pointInItem,
  toItemLocal, fromItemLocal, itemCorners, openingGeom, clampOpeningT, wallLength, wallAngle,
  wallPieces, fmtLen, fmtArea, parseLen,
} from './model.js';
import { ITEMS } from './catalog.js';

const HANDLE_PX = 7;
const SNAP_PX = 12;
const ROT_HANDLE_PX = 26;
const OVERHEAD = 1.4; // items whose bottom is above this are drawn dashed (overhead)

export class PlanView {
  constructor(canvas, app) {
    this.cv = canvas;
    this.ctx = canvas.getContext('2d');
    this.app = app;
    this.scale = 60; this.ox = 80; this.oy = 80;
    this.dpr = window.devicePixelRatio || 1;
    this.mouse = { sx: 0, sy: 0, x: 0, y: 0, inside: false };
    this.drag = null;
    this.wallStart = null;
    this.roomPts = [];
    this.boxStart = null;
    this.dimStart = null;
    this.ghostRot = 0;
    this.spaceDown = false;
    this.pending = false;
    this.hover = null;

    canvas.addEventListener('pointerdown', (e) => this.onDown(e));
    canvas.addEventListener('pointermove', (e) => this.onMove(e));
    canvas.addEventListener('pointerup', (e) => this.onUp(e));
    canvas.addEventListener('pointerleave', () => { this.mouse.inside = false; this.requestRender(); });
    canvas.addEventListener('wheel', (e) => this.onWheel(e), { passive: false });
    canvas.addEventListener('dblclick', (e) => this.onDbl(e));
    canvas.addEventListener('contextmenu', (e) => { e.preventDefault(); this.finishTool(); });
    new ResizeObserver(() => this.resize()).observe(canvas.parentElement);
    this.resize();
  }

  // ---------- coordinates ----------
  get floor() { return this.app.floor; }
  get px() { return 1 / this.scale; }
  w2s(x, y) { return [x * this.scale + this.ox, y * this.scale + this.oy]; }
  s2w(sx, sy) { return [(sx - this.ox) / this.scale, (sy - this.oy) / this.scale]; }

  resize() {
    const r = this.cv.parentElement.getBoundingClientRect();
    if (!r.width || !r.height) return;
    this.dpr = window.devicePixelRatio || 1;
    this.cv.width = Math.round(r.width * this.dpr);
    this.cv.height = Math.round(r.height * this.dpr);
    this.cv.style.width = r.width + 'px';
    this.cv.style.height = r.height + 'px';
    this.W = r.width; this.H = r.height;
    if (!this.app.design) return;
    if (!this._fitted) { this._fitted = true; this.fit(); }
    this.render();
  }

  bounds(floor = this.floor) {
    let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
    const add = (x, y) => { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); };
    for (const w of floor.walls) { add(w.x1, w.y1); add(w.x2, w.y2); }
    for (const r of floor.rooms) r.points.forEach(([x, y]) => add(x, y));
    for (const it of floor.items) itemCorners(it).forEach(([x, y]) => add(x, y));
    for (const n of floor.notes) { add(n.x, n.y); if (n.kind === 'dim') add(n.x2, n.y2); }
    if (!isFinite(x0)) return null;
    return { x0, y0, x1, y1 };
  }

  fit() {
    if (!this.W) return;
    const b = this.bounds() || { x0: 0, y0: 0, x1: 12, y1: 9 };
    const pad = 50;
    const bw = Math.max(b.x1 - b.x0, 2), bh = Math.max(b.y1 - b.y0, 2);
    this.scale = Math.max(8, Math.min(250, Math.min((this.W - pad * 2) / bw, (this.H - pad * 2) / bh)));
    this.ox = (this.W - bw * this.scale) / 2 - b.x0 * this.scale;
    this.oy = (this.H - bh * this.scale) / 2 - b.y0 * this.scale;
    this.requestRender();
  }

  zoomBy(f, sx = this.W / 2, sy = this.H / 2) {
    const [wx, wy] = this.s2w(sx, sy);
    this.scale = Math.max(5, Math.min(600, this.scale * f));
    this.ox = sx - wx * this.scale;
    this.oy = sy - wy * this.scale;
    this.requestRender();
  }

  focusOn(x, y) {
    this.ox = this.W / 2 - x * this.scale;
    this.oy = this.H / 2 - y * this.scale;
    this.requestRender();
  }

  requestRender() {
    if (this.pending) return;
    this.pending = true;
    requestAnimationFrame(() => { this.pending = false; this.render(); });
  }

  // ---------- snapping ----------
  snap(x, y, e, { endpoints = true } = {}) {
    if (e && e.altKey) return [x, y];
    if (endpoints) {
      let best = null, bd = SNAP_PX / this.scale;
      const test = (px, py) => { const d = dist(x, y, px, py); if (d < bd) { bd = d; best = [px, py]; } };
      for (const w of this.floor.walls) { test(w.x1, w.y1); test(w.x2, w.y2); }
      for (const r of this.floor.rooms) r.points.forEach(([px, py]) => test(px, py));
      if (best) return best;
    }
    const s = this.app.design.settings;
    if (!s.snap) return [x, y];
    const g = s.grid || 0.1;
    return [Math.round(x / g) * g, Math.round(y / g) * g];
  }

  // Shift = constrain to 0/45/90 degrees relative to a start point.
  constrain(sx, sy, x, y, e) {
    if (!e || !e.shiftKey) return [x, y];
    const a = Math.atan2(y - sy, x - sx);
    const step = Math.PI / 4;
    const ra = Math.round(a / step) * step;
    const L = dist(sx, sy, x, y);
    return [sx + Math.cos(ra) * L, sy + Math.sin(ra) * L];
  }

  nearestWall(x, y, maxD = 1.0) {
    let best = null;
    for (const w of this.floor.walls) {
      const p = projectOnSegment(x, y, w.x1, w.y1, w.x2, w.y2);
      if (p.d < maxD && (!best || p.d < best.p.d)) best = { w, p };
    }
    return best;
  }

  // Put a wall-mounted item flat against the nearest wall, facing into the room.
  wallSnapItem(it, x, y) {
    const nw = this.nearestWall(x, y, Math.max(0.8, it.d + 0.3));
    if (!nw) return false;
    const { w, p } = nw;
    const ang = wallAngle(w);
    let nx = -Math.sin(ang), ny = Math.cos(ang);
    if ((x - p.x) * nx + (y - p.y) * ny < 0) { nx = -nx; ny = -ny; }
    const off = (w.thickness || 0.15) / 2 + it.d / 2;
    it.x = p.x + nx * off; it.y = p.y + ny * off;
    it.rot = Math.round(((Math.atan2(-nx, ny) * 180) / Math.PI + 360) % 360);
    return true;
  }

  // ---------- hit testing ----------
  hitTest(x, y) {
    const f = this.floor, tol = 6 / this.scale;
    const sel = this.app.selection;
    // handles of current selection first
    if (sel) {
      const h = this.hitHandle(sel, x, y);
      if (h) return h;
    }
    for (let i = f.notes.length - 1; i >= 0; i--) {
      const n = f.notes[i];
      if (n.kind === 'dim') {
        if (projectOnSegment(x, y, n.x, n.y, n.x2, n.y2).d < tol * 1.5) return { type: 'note', id: n.id, part: 'body' };
      } else if (dist(x, y, n.x, n.y) < 14 / this.scale) return { type: 'note', id: n.id, part: 'body' };
    }
    for (const o of f.openings) {
      const w = f.walls.find((w) => w.id === o.wallId);
      if (!w) continue;
      const g = openingGeom(o, w);
      const [lx, ly] = [(x - g.cx) * g.ux + (y - g.cy) * g.uy, (x - g.cx) * g.nx + (y - g.cy) * g.ny];
      if (Math.abs(lx) <= o.width / 2 && Math.abs(ly) <= (w.thickness || 0.15) / 2 + tol) return { type: 'opening', id: o.id, part: 'body' };
    }
    const items = this.sortedItems().slice().reverse();
    for (const it of items) if (!it.locked && pointInItem(it, x, y, tol * 0.5)) return { type: 'item', id: it.id, part: 'body' };
    for (const w of f.walls) {
      if (projectOnSegment(x, y, w.x1, w.y1, w.x2, w.y2).d <= (w.thickness || 0.15) / 2 + tol) return { type: 'wall', id: w.id, part: 'body' };
    }
    for (const it of items) if (it.locked && pointInItem(it, x, y)) return { type: 'item', id: it.id, part: 'body' };
    for (let i = f.rooms.length - 1; i >= 0; i--) {
      if (pointInPolygon(x, y, f.rooms[i].points)) return { type: 'room', id: f.rooms[i].id, part: 'body' };
    }
    return null;
  }

  hitHandle(sel, x, y) {
    const r = HANDLE_PX / this.scale * 1.4;
    const obj = this.app.getSelected();
    if (!obj) return null;
    if (sel.type === 'wall') {
      if (dist(x, y, obj.x1, obj.y1) < r) return { ...sel, part: 'p1' };
      if (dist(x, y, obj.x2, obj.y2) < r) return { ...sel, part: 'p2' };
    } else if (sel.type === 'room') {
      for (let i = 0; i < obj.points.length; i++) {
        if (dist(x, y, obj.points[i][0], obj.points[i][1]) < r) return { ...sel, part: 'v' + i };
      }
    } else if (sel.type === 'item' && !obj.locked) {
      const [rx, ry] = fromItemLocal(obj, 0, -obj.d / 2 - ROT_HANDLE_PX / this.scale);
      if (dist(x, y, rx, ry) < r) return { ...sel, part: 'rot' };
      const [cx, cy] = fromItemLocal(obj, obj.w / 2, obj.d / 2);
      if (dist(x, y, cx, cy) < r) return { ...sel, part: 'size' };
    } else if (sel.type === 'note' && obj.kind === 'dim') {
      if (dist(x, y, obj.x, obj.y) < r) return { ...sel, part: 'p1' };
      if (dist(x, y, obj.x2, obj.y2) < r) return { ...sel, part: 'p2' };
    }
    return null;
  }

  sortedItems() {
    return this.floor.items.slice().sort((a, b) => {
      const ka = (ITEMS[a.catalogId]?.model === 'rug' || ITEMS[a.catalogId]?.model === 'rrug' ? 0 : 1) + (a.elev >= OVERHEAD ? 2 : 0);
      const kb = (ITEMS[b.catalogId]?.model === 'rug' || ITEMS[b.catalogId]?.model === 'rrug' ? 0 : 1) + (b.elev >= OVERHEAD ? 2 : 0);
      return ka - kb;
    });
  }

  // ---------- events ----------
  evPos(e) {
    const r = this.cv.getBoundingClientRect();
    const sx = e.clientX - r.left, sy = e.clientY - r.top;
    const [x, y] = this.s2w(sx, sy);
    this.mouse = { sx, sy, x, y, inside: true };
    return [x, y];
  }

  onWheel(e) {
    e.preventDefault();
    const r = this.cv.getBoundingClientRect();
    if (e.ctrlKey || !e.shiftKey) this.zoomBy(Math.exp(-e.deltaY * 0.0015), e.clientX - r.left, e.clientY - r.top);
    else { this.ox -= e.deltaY; this.requestRender(); }
  }

  onDown(e) {
    this.cv.setPointerCapture(e.pointerId);
    const [x, y] = this.evPos(e);
    const app = this.app, tool = app.tool;
    if (e.button === 1 || this.spaceDown || tool === 'pan') {
      this.drag = { kind: 'pan', sx: e.clientX, sy: e.clientY, ox: this.ox, oy: this.oy };
      return;
    }
    if (e.button !== 0) return;

    if (tool === 'calib') {
      if (!this.calibStart) { this.calibStart = [x, y]; this.requestRender(); return; }
      const a = this.calibStart; this.calibStart = null;
      this.app.calibrateUnderlay(a, [x, y]);
      return;
    }
    if (tool === 'movebg') {
      const u = this.floor.underlay;
      if (u) this.drag = { kind: 'bg', x0: x, y0: y, ux: u.x, uy: u.y };
      return;
    }
    if (tool === 'select') {
      const hit = this.hitTest(x, y);
      if (hit) {
        if (!app.selection || app.selection.id !== hit.id) app.select(hit.type, hit.id);
        const obj = app.getSelected();
        this.drag = { kind: 'move', hit, x0: x, y0: y, orig: JSON.parse(JSON.stringify(obj)), moved: false };
        if (hit.type === 'wall') this.drag.linked = this.linkedEnds(obj, hit.part);
        if (hit.type === 'room' && hit.part === 'body') this.drag.linkedWalls = null;
      } else {
        app.select(null);
        this.drag = { kind: 'pan', sx: e.clientX, sy: e.clientY, ox: this.ox, oy: this.oy };
      }
      return;
    }
    if (tool === 'wall') {
      let [px, py] = this.snap(x, y, e);
      if (this.wallStart) {
        [px, py] = this.constrain(this.wallStart[0], this.wallStart[1], px, py, e);
        if (dist(px, py, ...this.wallStart) > 0.05) {
          app.floor.walls.push({ id: uid('w'), x1: this.wallStart[0], y1: this.wallStart[1], x2: px, y2: py, ...app.newWallProps() });
          app.commit('Add wall');
        }
      }
      this.wallStart = [px, py];
      this.requestRender();
      return;
    }
    if (tool === 'room') {
      const [px, py] = this.snap(x, y, e);
      if (this.roomPts.length >= 3 && dist(px, py, ...this.roomPts[0]) < SNAP_PX / this.scale) { this.closeRoom(); return; }
      this.roomPts.push([px, py]);
      this.requestRender();
      return;
    }
    if (tool === 'roombox' || tool === 'area') {
      this.boxStart = this.snap(x, y, e);
      this.drag = { kind: 'box' };
      return;
    }
    if (tool === 'place') { this.placeAt(x, y, e); return; }
    if (tool === 'pin' || tool === 'label') {
      const [px, py] = this.snap(x, y, e, { endpoints: false });
      const n = tool === 'pin'
        ? { id: uid('n'), kind: 'pin', x: px, y: py, title: 'Note', text: '', color: '#e53935' }
        : { id: uid('n'), kind: 'label', x: px, y: py, text: 'Label', color: '#263238', size: 14 };
      app.floor.notes.push(n);
      app.commit(tool === 'pin' ? 'Add note' : 'Add label');
      app.setTool('select');
      app.select('note', n.id, { focus: true });
      return;
    }
    if (tool === 'dim') {
      const [px, py] = this.snap(x, y, e);
      if (!this.dimStart) { this.dimStart = [px, py]; this.requestRender(); return; }
      const [qx, qy] = this.constrain(this.dimStart[0], this.dimStart[1], px, py, e);
      if (dist(qx, qy, ...this.dimStart) > 0.02) {
        app.floor.notes.push({ id: uid('n'), kind: 'dim', x: this.dimStart[0], y: this.dimStart[1], x2: qx, y2: qy, text: '', color: '#d81b60', offset: 0.3 });
        app.commit('Add dimension');
      }
      this.dimStart = null;
      this.requestRender();
    }
  }

  onMove(e) {
    const [x, y] = this.evPos(e);
    const d = this.drag;
    if (d && d.kind === 'pan') {
      this.ox = d.ox + (e.clientX - d.sx);
      this.oy = d.oy + (e.clientY - d.sy);
      this.requestRender();
      return;
    }
    if (d && d.kind === 'move') { this.dragMove(d, x, y, e); return; }
    if (d && d.kind === 'bg') { const u = this.floor.underlay; u.x = d.ux + (x - d.x0); u.y = d.uy + (y - d.y0); this.requestRender(); return; }
    this.lastEvent = e;
    if (this.app.tool === 'select') {
      const h = this.hitTest(x, y);
      const key = h ? h.type + h.id + h.part : '';
      if (key !== this._hoverKey) { this._hoverKey = key; this.hover = h; }
      this.cv.style.cursor = h ? (h.part === 'body' ? 'move' : 'crosshair') : 'grab';
    } else {
      this.cv.style.cursor = 'crosshair';
    }
    this.app.status(`x ${x.toFixed(2)}  y ${y.toFixed(2)} m`);
    this.requestRender();
  }

  onUp(e) {
    const d = this.drag;
    this.drag = null;
    if (!d) return;
    if (d.kind === 'move' && d.moved) this.app.commit('Move');
    if (d.kind === 'bg') this.app.commit('Move trace image');
    if (d.kind === 'box' && this.boxStart) {
      const [x, y] = this.snap(this.mouse.x, this.mouse.y, e);
      const [x0, y0] = this.boxStart;
      this.boxStart = null;
      if (Math.abs(x - x0) > 0.2 && Math.abs(y - y0) > 0.2) this.makeBoxRoom(x0, y0, x, y, this.app.tool === 'roombox');
    }
    this.requestRender();
  }

  onDbl(e) {
    const tool = this.app.tool;
    if (tool === 'wall' || tool === 'room') { this.finishTool(); return; }
    if (tool === 'select' && this.app.selection) this.app.focusProps();
  }

  finishTool() {
    if (this.app.tool === 'room' && this.roomPts.length >= 3) { this.closeRoom(); return; }
    const had = this.wallStart || this.roomPts.length || this.dimStart || this.boxStart || this.calibStart;
    this.wallStart = null; this.roomPts = []; this.dimStart = null; this.boxStart = null; this.calibStart = null;
    if (!had && this.app.tool !== 'select') this.app.setTool('select');
    this.requestRender();
  }

  // ---------- creation helpers ----------
  closeRoom() {
    const pts = this.roomPts;
    this.roomPts = [];
    this.addRoom(pts);
  }

  addRoom(points) {
    const n = this.floor.rooms.length + 1;
    const room = { id: uid('r'), name: `Room ${n}`, type: 'Other', points, color: '#f5f5f5', notes: '', floorFinish: '' };
    this.floor.rooms.push(room);
    this.app.commit('Add room');
    this.app.setTool('select');
    this.app.select('room', room.id);
  }

  makeBoxRoom(x0, y0, x1, y1, withWalls) {
    const xa = Math.min(x0, x1), xb = Math.max(x0, x1), ya = Math.min(y0, y1), yb = Math.max(y0, y1);
    if (withWalls) {
      const pts = [[xa, ya], [xb, ya], [xb, yb], [xa, yb]];
      for (let i = 0; i < 4; i++) {
        const [ax, ay] = pts[i], [bx, by] = pts[(i + 1) % 4];
        const exists = this.floor.walls.some((w) =>
          (dist(w.x1, w.y1, ax, ay) < 1e-3 && dist(w.x2, w.y2, bx, by) < 1e-3) ||
          (dist(w.x1, w.y1, bx, by) < 1e-3 && dist(w.x2, w.y2, ax, ay) < 1e-3));
        if (!exists) this.floor.walls.push({ id: uid('w'), x1: ax, y1: ay, x2: bx, y2: by, ...this.app.newWallProps() });
      }
    }
    this.addRoom([[xa, ya], [xb, ya], [xb, yb], [xa, yb]]);
  }

  ghostItem(x, y, e) {
    const def = this.app.placing;
    if (!def) return null;
    let [px, py] = this.snap(x, y, e, { endpoints: false });
    const it = {
      id: uid('i'), catalogId: def.id, x: px, y: py, rot: this.ghostRot,
      w: def.w, d: def.d, h: def.h, elev: def.elev || 0, color: def.color, label: '', notes: '',
    };
    if (def.wall && !(e && e.altKey)) this.wallSnapItem(it, x, y);
    return it;
  }

  ghostOpening(x, y) {
    const def = this.app.placing;
    const nw = this.nearestWall(x, y, 1.0);
    if (!nw) return null;
    const op = {
      id: uid('o'), catalogId: def.id, wallId: nw.w.id, t: nw.p.t, width: def.w, height: def.h,
      sill: def.sill ?? 0, hinge: 0, swing: 1, label: '', notes: '',
    };
    // swing to the side the mouse is on
    const g = openingGeom(op, nw.w);
    op.swing = (x - g.cx) * g.nx + (y - g.cy) * g.ny >= 0 ? 1 : -1;
    const s = this.app.design.settings.snap ? 0.05 : 0;
    if (s) { const L = wallLength(nw.w); op.t = Math.round((op.t * L) / s) * s / L; }
    op.t = clampOpeningT(op, nw.w);
    return op;
  }

  placeAt(x, y, e) {
    const def = this.app.placing;
    if (!def) return;
    if (def.kind === 'opening') {
      const op = this.ghostOpening(x, y);
      if (!op) { this.app.toast('Click on (or near) a wall to place a door or window'); return; }
      this.floor.openings.push(op);
      this.app.commit('Add ' + def.name);
      if (!e.shiftKey) { this.app.setTool('select'); this.app.select('opening', op.id); }
      return;
    }
    const it = this.ghostItem(x, y, e);
    this.floor.items.push(it);
    this.app.commit('Add ' + def.name);
    if (!e.shiftKey) { this.app.setTool('select'); this.app.select('item', it.id); }
  }

  // wall ends that share a point with the dragged end(s)
  linkedEnds(w, part) {
    const pts = part === 'p1' ? [[w.x1, w.y1]] : part === 'p2' ? [[w.x2, w.y2]] : [[w.x1, w.y1], [w.x2, w.y2]];
    const out = [];
    for (const o of this.floor.walls) {
      if (o.id === w.id) continue;
      for (const [px, py] of pts) {
        if (dist(o.x1, o.y1, px, py) < 1e-3) out.push({ w: o, end: 1, ox: o.x1, oy: o.y1 });
        if (dist(o.x2, o.y2, px, py) < 1e-3) out.push({ w: o, end: 2, ox: o.x2, oy: o.y2 });
      }
    }
    return out;
  }

  dragMove(d, x, y, e) {
    const obj = this.app.getSelected();
    if (!obj) return;
    const o = d.orig, part = d.hit.part, type = d.hit.type;
    let dx = x - d.x0, dy = y - d.y0;
    if (!d.moved && Math.hypot(dx, dy) * this.scale < 3) return;
    d.moved = true;
    const g = this.app.design.settings.snap && !e.altKey ? this.app.design.settings.grid || 0.1 : 0;
    const rs = (v) => (g ? Math.round(v / g) * g : v);

    if (type === 'item') {
      if (part === 'rot') {
        let a = (Math.atan2(y - obj.y, x - obj.x) * 180) / Math.PI + 90;
        a = e.shiftKey ? a : Math.round(a / 15) * 15;
        obj.rot = Math.round(((a % 360) + 360) % 360);
      } else if (part === 'size') {
        const [lx, ly] = toItemLocal(o, x, y);
        const nw = Math.max(0.05, rs(lx + o.w / 2)), nd = Math.max(0.05, rs(ly + o.d / 2));
        const [cx, cy] = fromItemLocal(o, -o.w / 2 + nw / 2, -o.d / 2 + nd / 2);
        obj.w = +nw.toFixed(3); obj.d = +nd.toFixed(3); obj.x = cx; obj.y = cy;
      } else {
        obj.x = o.x + dx; obj.y = o.y + dy;
        const def = ITEMS[obj.catalogId];
        if (def?.wall && !e.altKey && this.wallSnapItem(obj, x, y)) { /* snapped */ }
        else { obj.x = rs(obj.x); obj.y = rs(obj.y); }
      }
    } else if (type === 'opening') {
      const nw = this.nearestWall(x, y, 1.5);
      if (nw) {
        obj.wallId = nw.w.id;
        const L = wallLength(nw.w);
        obj.t = g ? Math.round((nw.p.t * L) / 0.05) * 0.05 / L : nw.p.t;
        obj.t = clampOpeningT(obj, nw.w);
      }
    } else if (type === 'wall') {
      let ax = o.x1, ay = o.y1, bx = o.x2, by = o.y2;
      if (part === 'p1') { [ax, ay] = this.snapExcluding(x, y, e, obj); }
      else if (part === 'p2') { [bx, by] = this.snapExcluding(x, y, e, obj); }
      else { dx = rs(dx); dy = rs(dy); ax += dx; ay += dy; bx += dx; by += dy; }
      obj.x1 = ax; obj.y1 = ay; obj.x2 = bx; obj.y2 = by;
      for (const l of d.linked) {
        const nx = l.ox + (part === 'p2' || (part === 'body' && dist(l.ox, l.oy, o.x2, o.y2) < 1e-3) ? bx - o.x2 : ax - o.x1);
        const ny = l.oy + (part === 'p2' || (part === 'body' && dist(l.ox, l.oy, o.x2, o.y2) < 1e-3) ? by - o.y2 : ay - o.y1);
        if (l.end === 1) { l.w.x1 = nx; l.w.y1 = ny; } else { l.w.x2 = nx; l.w.y2 = ny; }
      }
    } else if (type === 'room') {
      if (part.startsWith('v')) {
        const i = +part.slice(1);
        obj.points[i] = this.snap(x, y, e);
      } else {
        dx = rs(dx); dy = rs(dy);
        obj.points = o.points.map(([px, py]) => [px + dx, py + dy]);
      }
    } else if (type === 'note') {
      if (part === 'p1') [obj.x, obj.y] = this.snap(x, y, e);
      else if (part === 'p2') [obj.x2, obj.y2] = this.snap(x, y, e);
      else {
        dx = rs(dx); dy = rs(dy);
        obj.x = o.x + dx; obj.y = o.y + dy;
        if (o.kind === 'dim') { obj.x2 = o.x2 + dx; obj.y2 = o.y2 + dy; }
      }
    }
    this.app.liveUpdate();
  }

  snapExcluding(x, y, e, wall) {
    const saved = [wall.x1, wall.y1, wall.x2, wall.y2];
    wall.x1 = wall.y1 = wall.x2 = wall.y2 = 1e9;
    const r = this.snap(x, y, e);
    [wall.x1, wall.y1, wall.x2, wall.y2] = saved;
    return r;
  }

  // ---------- drawing ----------
  render(opts = {}) {
    const ctx = this.ctx, app = this.app;
    if (!this.W || !app.design) return;
    const exporting = !!opts.exporting;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    ctx.fillStyle = exporting ? '#ffffff' : '#f6f7f9';
    ctx.fillRect(0, 0, this.W, this.H);
    if (!exporting) this.drawGrid();

    const f = this.floor, s = app.design.settings, units = s.units;
    const wt = () => ctx.setTransform(this.dpr * this.scale, 0, 0, this.dpr * this.scale, this.dpr * this.ox, this.dpr * this.oy);
    const st = () => ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    const px = this.px;

    // trace image (underlay)
    const u = f.underlay;
    if (u && u.visible !== false && !exporting) {
      const img = this.getImage(this.app.assetSrc(u.src));
      if (img && img.complete && img.naturalWidth) {
        wt();
        ctx.globalAlpha = u.opacity ?? 0.5;
        ctx.drawImage(img, u.x, u.y, img.naturalWidth * u.mpp, img.naturalHeight * u.mpp);
        ctx.globalAlpha = 1;
      }
    }

    // faint outline of the floor below, to line up walls between levels
    const fi = app.design.floors.indexOf(f);
    if (fi > 0 && s.showBelow !== false) {
      wt();
      ctx.strokeStyle = 'rgba(120,144,156,0.35)';
      ctx.lineWidth = 0.1;
      for (const w of app.design.floors[fi - 1].walls) {
        ctx.beginPath(); ctx.moveTo(w.x1, w.y1); ctx.lineTo(w.x2, w.y2); ctx.stroke();
      }
    }

    // rooms
    wt();
    for (const r of f.rooms) {
      if (r.points.length < 3) continue;
      ctx.beginPath();
      r.points.forEach(([x, y], i) => (i ? ctx.lineTo(x, y) : ctx.moveTo(x, y)));
      ctx.closePath();
      ctx.fillStyle = r.color || '#f5f5f5';
      ctx.fill();
      ctx.strokeStyle = 'rgba(0,0,0,0.12)';
      ctx.lineWidth = px;
      ctx.stroke();
    }

    // items
    for (const it of this.sortedItems()) this.drawItem(it);

    // walls
    wt();
    const opsByWall = f.openings;
    for (const w of f.walls) {
      const th = w.thickness || 0.15;
      const { L, solids, gaps } = wallPieces(w, opsByWall);
      const a = wallAngle(w);
      ctx.save();
      ctx.translate(w.x1, w.y1);
      ctx.rotate(a);
      ctx.fillStyle = w.kind === 'interior' ? '#78909c' : '#37474f';
      for (const [s0, s1] of solids) {
        const e0 = s0 <= 1e-4 ? -th / 2 : s0, e1 = s1 >= L - 1e-4 ? L + th / 2 : s1;
        ctx.fillRect(e0, -th / 2, e1 - e0, th);
      }
      for (const g of gaps) this.drawOpening(g.op, th, g.s0, g.s1);
      ctx.restore();
    }

    // dimension notes (world)
    for (const n of f.notes) if (n.kind === 'dim') this.drawDim(n);

    // ---------- screen-space overlays ----------
    st();
    // room labels
    for (const r of f.rooms) {
      if (r.points.length < 3 || r.hideLabel) continue;
      const [cx, cy] = polygonCentroid(r.points);
      const [sx, sy] = this.w2s(cx, cy);
      const area = polygonArea(r.points);
      ctx.textAlign = 'center';
      ctx.fillStyle = '#263238';
      ctx.font = `600 ${clamp(this.scale * 0.22, 10, 16)}px system-ui, sans-serif`;
      ctx.fillText(r.name, sx, sy);
      ctx.font = `${clamp(this.scale * 0.17, 9, 13)}px system-ui, sans-serif`;
      ctx.fillStyle = '#546e7a';
      ctx.fillText(fmtArea(area, units), sx, sy + clamp(this.scale * 0.22, 11, 17));
    }
    // item labels
    if (s.showLabels !== false) {
      for (const it of f.items) {
        const def = ITEMS[it.catalogId];
        if (it.w * this.scale < 50 || it.d * this.scale < 20) continue;
        if (!it.label && (this.scale < 75 || def?.model === 'rug' || def?.model === 'rrug')) continue;
        const [sx, sy] = this.w2s(it.x, it.y);
        ctx.font = `${clamp(this.scale * 0.13, 9, 12)}px system-ui, sans-serif`;
        ctx.textAlign = 'center';
        ctx.fillStyle = 'rgba(38,50,56,0.85)';
        const txt = it.label || def?.name || it.catalogId;
        ctx.fillText(ellipsize(ctx, txt, it.w * this.scale - 4), sx, sy + 4);
      }
    }
    // wall dimensions
    const sel = app.selection;
    for (const w of f.walls) {
      const selected = sel && sel.type === 'wall' && sel.id === w.id;
      if (!selected && !s.showDims) continue;
      this.drawWallDim(w, units, selected);
    }
    // pins & labels
    let pinNo = 0;
    for (const n of f.notes) {
      if (n.kind === 'pin') { pinNo++; this.drawPin(n, pinNo); }
      else if (n.kind === 'label') this.drawLabel(n);
    }

    if (!exporting) {
      this.drawSelection();
      this.drawToolPreview();
    }
  }

  drawGrid() {
    const ctx = this.ctx;
    const g = this.app.design.settings.grid || 0.1;
    let step = g;
    while (step * this.scale < 8) step *= 5;
    const major = step * 5 * this.scale < 400 ? 1 : 5;
    const [x0, y0] = this.s2w(0, 0), [x1, y1] = this.s2w(this.W, this.H);
    ctx.lineWidth = 1;
    const drawLines = (st, color) => {
      ctx.strokeStyle = color;
      ctx.beginPath();
      for (let x = Math.floor(x0 / st) * st; x <= x1; x += st) { const sx = Math.round(x * this.scale + this.ox) + 0.5; ctx.moveTo(sx, 0); ctx.lineTo(sx, this.H); }
      for (let y = Math.floor(y0 / st) * st; y <= y1; y += st) { const sy = Math.round(y * this.scale + this.oy) + 0.5; ctx.moveTo(0, sy); ctx.lineTo(this.W, sy); }
      ctx.stroke();
    };
    if (step * this.scale >= 8) drawLines(step, '#eceef1');
    drawLines(Math.max(1, step * 5) * major, '#dde1e6');
    // origin
    const [ox, oy] = this.w2s(0, 0);
    ctx.strokeStyle = '#b0bec5';
    ctx.beginPath(); ctx.moveTo(ox - 8, oy); ctx.lineTo(ox + 8, oy); ctx.moveTo(ox, oy - 8); ctx.lineTo(ox, oy + 8); ctx.stroke();
    // scale bar
    let bar = 1; while (bar * this.scale < 60) bar *= 2; while (bar * this.scale > 160) bar /= 2;
    ctx.fillStyle = '#546e7a'; ctx.strokeStyle = '#546e7a';
    const bx = 14, by = this.H - 16;
    ctx.beginPath(); ctx.moveTo(bx, by - 4); ctx.lineTo(bx, by); ctx.lineTo(bx + bar * this.scale, by); ctx.lineTo(bx + bar * this.scale, by - 4); ctx.stroke();
    ctx.font = '11px system-ui, sans-serif'; ctx.textAlign = 'left';
    ctx.fillText(fmtLen(bar, this.app.design.settings.units), bx + 4, by - 6);
  }

  drawOpening(op, th, s0, s1) {
    // drawn in wall-local coords: x along wall, y across (normal)
    const ctx = this.ctx, px = this.px, def = ITEMS[op.catalogId] || {};
    const w = s1 - s0, model = def.model || 'door';
    ctx.fillStyle = '#ffffff';
    ctx.fillRect(s0, -th / 2, w, th);
    ctx.strokeStyle = '#37474f';
    ctx.lineWidth = 1.2 * px;
    // jambs
    ctx.beginPath(); ctx.moveTo(s0, -th / 2); ctx.lineTo(s0, th / 2); ctx.moveTo(s1, -th / 2); ctx.lineTo(s1, th / 2); ctx.stroke();
    const sw = op.swing || 1;
    if (model === 'window') {
      ctx.fillStyle = 'rgba(144,202,249,0.55)';
      ctx.fillRect(s0, -th / 4, w, th / 2);
      ctx.beginPath(); ctx.moveTo(s0, -th / 4); ctx.lineTo(s1, -th / 4); ctx.moveTo(s0, 0); ctx.lineTo(s1, 0); ctx.moveTo(s0, th / 4); ctx.lineTo(s1, th / 4); ctx.stroke();
    } else if (model === 'slider') {
      ctx.lineWidth = 2 * px;
      ctx.strokeStyle = shade(def.color || '#90caf9', -0.3);
      ctx.beginPath(); ctx.moveTo(s0, -th / 6); ctx.lineTo(s0 + w * 0.55, -th / 6); ctx.moveTo(s1, th / 6); ctx.lineTo(s1 - w * 0.55, th / 6); ctx.stroke();
    } else if (model === 'arch') {
      ctx.setLineDash([4 * px, 3 * px]);
      ctx.beginPath(); ctx.moveTo(s0, -th / 2); ctx.lineTo(s1, -th / 2); ctx.moveTo(s0, th / 2); ctx.lineTo(s1, th / 2); ctx.stroke();
      ctx.setLineDash([]);
    } else {
      const leaves = model === 'door2' ? [[s0, w / 2, 1], [s1, w / 2, -1]] : [[op.hinge ? s1 : s0, w, op.hinge ? -1 : 1]];
      ctx.lineWidth = 1.5 * px;
      for (const [hx, lw, dir] of leaves) {
        const hy = (sw * th) / 2;
        ctx.beginPath(); ctx.moveTo(hx, hy); ctx.lineTo(hx, hy + sw * lw); ctx.stroke();
        ctx.save();
        ctx.lineWidth = 0.8 * px;
        ctx.setLineDash([3 * px, 2 * px]);
        ctx.beginPath();
        const a0 = sw > 0 ? Math.PI / 2 : -Math.PI / 2;
        const a1 = dir > 0 ? 0 : Math.PI;
        ctx.arc(hx, hy, lw, a0, a1, (sw > 0) === (dir > 0));
        ctx.stroke();
        ctx.restore();
      }
    }
  }

  drawItem(it) {
    const ctx = this.ctx, px = this.px, def = ITEMS[it.catalogId] || { model: 'box', color: '#bdbdbd' };
    const model = def.model, w = it.w, d = it.d, hw = w / 2, hd = d / 2;
    const color = it.color || def.color;
    const overhead = it.elev >= OVERHEAD && def.kind !== 'opening';
    const edge = shade(color, -0.45);
    ctx.setTransform(this.dpr * this.scale, 0, 0, this.dpr * this.scale, this.dpr * this.ox, this.dpr * this.oy);
    ctx.save();
    ctx.translate(it.x, it.y);
    ctx.rotate((it.rot * Math.PI) / 180);
    ctx.lineWidth = px;
    ctx.strokeStyle = edge;
    ctx.fillStyle = color;
    if (overhead) { ctx.globalAlpha = 0.55; ctx.setLineDash([4 * px, 3 * px]); }
    const rect = (x, y, ww, hh, fill = true) => { if (fill) ctx.fillRect(x, y, ww, hh); ctx.strokeRect(x, y, ww, hh); };
    const ell = (cx, cy, rx, ry, fill = true) => { ctx.beginPath(); ctx.ellipse(cx, cy, Math.abs(rx), Math.abs(ry), 0, 0, Math.PI * 2); if (fill) ctx.fill(); ctx.stroke(); };
    const rrect = (x, y, ww, hh, r, fill = true) => { ctx.beginPath(); ctx.roundRect(x, y, ww, hh, r); if (fill) ctx.fill(); ctx.stroke(); };

    switch (model) {
      case 'rug': case 'pool': case 'shower':
        ctx.globalAlpha = model === 'rug' ? 0.6 : 0.75;
        rect(-hw, -hd, w, d);
        ctx.globalAlpha = 1;
        if (model === 'rug') { ctx.setLineDash([3 * px, 3 * px]); rect(-hw + 0.05, -hd + 0.05, w - 0.1, d - 0.1, false); }
        if (model === 'shower') { ell(0, 0, 0.06, 0.06, false); ctx.beginPath(); ctx.moveTo(-hw, -hd); ctx.lineTo(hw, hd); ctx.moveTo(hw, -hd); ctx.lineTo(-hw, hd); ctx.stroke(); }
        break;
      case 'rrug':
        ctx.globalAlpha = 0.6; ell(0, 0, hw, hd); ctx.globalAlpha = 1;
        break;
      case 'bed': case 'bunk': case 'crib': {
        rrect(-hw, -hd, w, d, 0.04);
        ctx.fillStyle = '#ffffff';
        if (model === 'crib') { rect(-hw + 0.04, -hd + 0.04, w - 0.08, d - 0.08); break; }
        const pw = w > 1.2 ? (w - 0.3) / 2 : w - 0.2;
        const n = w > 1.2 ? 2 : 1;
        for (let i = 0; i < n; i++) rrect(-hw + 0.1 + i * (pw + 0.1), -hd + 0.08, pw, Math.min(0.35, d * 0.18), 0.05);
        ctx.fillStyle = shade(color, 0.25);
        rect(-hw + 0.03, -hd + d * 0.32, w - 0.06, d * 0.66);
        ctx.beginPath(); ctx.moveTo(-hw + 0.03, -hd + d * 0.42); ctx.lineTo(hw - 0.03, -hd + d * 0.42); ctx.stroke();
        break;
      }
      case 'sofa': {
        const bt = Math.min(0.22, d * 0.25), at = Math.min(0.18, w * 0.12);
        rrect(-hw, -hd, w, d, 0.06);
        ctx.fillStyle = shade(color, -0.15);
        rect(-hw, -hd, w, bt);
        if (w > 0.6) { rect(-hw, -hd, at, d); rect(hw - at, -hd, at, d); }
        const seats = Math.max(1, Math.round((w - 2 * at) / 0.65));
        ctx.beginPath();
        for (let i = 1; i < seats; i++) { const x = -hw + at + ((w - 2 * at) * i) / seats; ctx.moveTo(x, -hd + bt); ctx.lineTo(x, hd); }
        ctx.stroke();
        break;
      }
      case 'lsofa': {
        const sd = 0.9, bt = 0.22;
        ctx.beginPath();
        ctx.moveTo(-hw, -hd); ctx.lineTo(hw, -hd); ctx.lineTo(hw, hd); ctx.lineTo(hw - sd, hd); ctx.lineTo(hw - sd, -hd + sd); ctx.lineTo(-hw, -hd + sd); ctx.closePath();
        ctx.fill(); ctx.stroke();
        ctx.fillStyle = shade(color, -0.15);
        rect(-hw, -hd, w, bt); rect(hw - bt, -hd, bt, d);
        break;
      }
      case 'table': case 'dining': case 'ldesk': {
        if (model === 'dining') this.drawChairsRect(def.seats || 4, w, d, color);
        ctx.fillStyle = color;
        if (model === 'ldesk') {
          ctx.beginPath(); ctx.moveTo(-hw, -hd); ctx.lineTo(hw, -hd); ctx.lineTo(hw, hd); ctx.lineTo(hw - 0.6, hd); ctx.lineTo(hw - 0.6, -hd + 0.6); ctx.lineTo(-hw, -hd + 0.6); ctx.closePath(); ctx.fill(); ctx.stroke();
        } else {
          rrect(-hw, -hd, w, d, 0.03);
          rect(-hw + 0.05, -hd + 0.05, w - 0.1, d - 0.1, false);
        }
        break;
      }
      case 'rtable': case 'rdining':
        if (model === 'rdining') this.drawChairsRound(def.seats || 4, hw, hd, color);
        ctx.fillStyle = color;
        ell(0, 0, hw, hd);
        ell(0, 0, hw * 0.85, hd * 0.85, false);
        break;
      case 'chair': case 'stool':
        if (model === 'stool') { ell(0, 0, hw, hd); break; }
        rrect(-hw, -hd, w, d, 0.05);
        ctx.fillStyle = shade(color, -0.2);
        rect(-hw, -hd, w, Math.min(0.1, d * 0.2));
        break;
      case 'counter': case 'vanity': case 'ksink': case 'hob': case 'basin': {
        if (model !== 'basin') rect(-hw, -hd, w, d);
        ctx.fillStyle = '#ffffff';
        if (model === 'ksink') {
          const n = w >= 1.1 ? 2 : 1, bw = (w - 0.15 - 0.05 * (n - 1)) / n;
          for (let i = 0; i < n; i++) rrect(-hw + 0.075 + i * (bw + 0.05), -hd + 0.1, bw, d - 0.2, 0.05);
        } else if (model === 'hob') {
          ctx.fillStyle = '#263238';
          const r = Math.min(w, d) * 0.15;
          const cols = w < 0.72 && d < 0.5 ? [[-w / 4, 0], [w / 4, 0]] : [[-w / 4, -d / 4], [w / 4, -d / 4], [-w / 4, d / 4], [w / 4, d / 4]];
          ctx.strokeStyle = '#90a4ae';
          for (const [cx, cy] of cols) ell(cx, cy, r, r, false);
        } else if (model === 'vanity') {
          const n = w >= 1.3 ? 2 : 1;
          for (let i = 0; i < n; i++) ell(-hw + (w * (i + 0.5)) / n, 0, Math.min(0.22, w / n / 2 - 0.08), d * 0.3);
        } else if (model === 'basin') {
          ctx.fillStyle = color;
          rrect(-hw, -hd, w, d, 0.08);
          ctx.fillStyle = '#ffffff';
          ell(0, d * 0.05, hw * 0.7, hd * 0.6);
        } else {
          ctx.beginPath(); ctx.moveTo(-hw, hd - 0.04); ctx.lineTo(hw, hd - 0.04); ctx.stroke();
        }
        break;
      }
      case 'toilet':
        rrect(-hw, -hd, w, Math.min(0.2, d * 0.3), 0.03);
        ell(0, -hd + d * 0.6, hw * 0.85, d * 0.38);
        break;
      case 'tub': case 'enclosure':
        if (model === 'enclosure') { ctx.fillStyle = 'rgba(179,229,252,0.5)'; rect(-hw, -hd, w, d); ell(0, 0, 0.06, 0.06, false); break; }
        rrect(-hw, -hd, w, d, 0.08);
        ctx.fillStyle = '#e1f5fe';
        rrect(-hw + 0.08, -hd + 0.08, w - 0.16, d - 0.16, Math.min(w, d) * 0.3);
        break;
      case 'stairs': case 'ramp': {
        rect(-hw, -hd, w, d);
        ctx.beginPath();
        if (model === 'stairs') {
          const n = Math.max(2, Math.round(d / 0.28));
          for (let i = 1; i < n; i++) { const y = -hd + (d * i) / n; ctx.moveTo(-hw, y); ctx.lineTo(hw, y); }
        }
        ctx.moveTo(0, hd - 0.1); ctx.lineTo(0, -hd + 0.15);
        ctx.moveTo(-0.1, -hd + 0.3); ctx.lineTo(0, -hd + 0.15); ctx.lineTo(0.1, -hd + 0.3);
        ctx.stroke();
        break;
      }
      case 'spiral': {
        ell(0, 0, hw, hd);
        ctx.beginPath();
        for (let i = 0; i < 12; i++) { const a = (i / 12) * Math.PI * 2; ctx.moveTo(0, 0); ctx.lineTo(Math.cos(a) * hw, Math.sin(a) * hd); }
        ctx.stroke();
        ell(0, 0, 0.08, 0.08);
        break;
      }
      case 'plant': case 'tree': case 'blob': case 'fountain': case 'umbrella': case 'cone': {
        ell(0, 0, hw, hd);
        if (model === 'plant' || model === 'tree') {
          ctx.beginPath();
          for (let i = 0; i < 8; i++) { const a = (i / 8) * Math.PI * 2; ctx.moveTo(0, 0); ctx.lineTo(Math.cos(a) * hw * 0.8, Math.sin(a) * hd * 0.8); }
          ctx.stroke();
        }
        if (model === 'fountain') { ctx.fillStyle = '#e1f5fe'; ell(0, 0, hw * 0.6, hd * 0.6); }
        break;
      }
      case 'cyl': case 'lamp': case 'light': case 'pendant': case 'chandelier':
        if (model === 'light' && (w / d > 2.5 || d / w > 2.5)) { rect(-hw, -hd, w, d); break; }
        ell(0, 0, hw, hd);
        if (model !== 'cyl') { ctx.beginPath(); ctx.moveTo(-hw * 0.6, 0); ctx.lineTo(hw * 0.6, 0); ctx.moveTo(0, -hd * 0.6); ctx.lineTo(0, hd * 0.6); ctx.stroke(); }
        break;
      case 'fan': {
        ell(0, 0, 0.1, 0.1);
        ctx.fillStyle = shade(color, 0.2);
        for (let i = 0; i < 3; i++) {
          ctx.save(); ctx.rotate((i * Math.PI * 2) / 3);
          rrect(0.08, -0.05, hw - 0.08, 0.1, 0.04);
          ctx.restore();
        }
        break;
      }
      case 'car': {
        rrect(-hw, -hd, w, d, 0.35);
        ctx.fillStyle = 'rgba(255,255,255,0.6)';
        rrect(-hw + 0.15, -hd + d * 0.25, w - 0.3, d * 0.12, 0.05);
        rrect(-hw + 0.18, -hd + d * 0.68, w - 0.36, d * 0.08, 0.05);
        break;
      }
      case 'washer':
        rect(-hw, -hd, w, d);
        ctx.fillStyle = '#b0bec5';
        ell(0, hd * 0.15, Math.min(hw, hd) * 0.6, Math.min(hw, hd) * 0.6);
        break;
      case 'shelf': case 'cabinet': case 'mandir':
        rect(-hw, -hd, w, d);
        ctx.beginPath();
        if (model === 'cabinet') { const n = Math.max(1, Math.round(w / 0.5)); for (let i = 1; i < n; i++) { const x = -hw + (w * i) / n; ctx.moveTo(x, -hd); ctx.lineTo(x, hd); } }
        if (model === 'shelf') { ctx.moveTo(-hw, -hd); ctx.lineTo(hw, hd); }
        if (model === 'mandir') { ctx.moveTo(-hw, hd); ctx.lineTo(0, -hd); ctx.lineTo(hw, hd); }
        ctx.stroke();
        break;
      case 'swing': case 'pergola':
        ctx.globalAlpha = 0.35; rect(-hw, -hd, w, d); ctx.globalAlpha = 1;
        rect(-hw, -hd, 0.12, 0.12); rect(hw - 0.12, -hd, 0.12, 0.12); rect(-hw, hd - 0.12, 0.12, 0.12); rect(hw - 0.12, hd - 0.12, 0.12, 0.12);
        if (model === 'swing') { ctx.fillStyle = color; rect(-hw * 0.7, -hd * 0.3, w * 0.7, d * 0.6); }
        break;
      case 'railing': case 'glass': case 'curtain':
        if (model === 'glass') ctx.fillStyle = 'rgba(179,229,252,0.7)';
        rect(-hw, -hd, w, d);
        if (model === 'curtain') { ctx.beginPath(); for (let x = -hw; x < hw; x += 0.1) { ctx.moveTo(x, -hd); ctx.lineTo(x + 0.05, hd); } ctx.stroke(); }
        break;
      case 'screen': case 'art': case 'hood':
        rect(-hw, -hd, w, d);
        break;
      case 'planter':
        rect(-hw, -hd, w, d);
        ctx.fillStyle = '#66bb6a';
        ell(0, 0, hw * 0.8, hd * 0.8);
        break;
      default:
        rect(-hw, -hd, w, d);
    }
    ctx.restore();
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;
  }

  drawChairsRect(seats, w, d, color) {
    const ctx = this.ctx, px = this.px;
    ctx.fillStyle = shade(color, 0.25);
    ctx.strokeStyle = shade(color, -0.4);
    ctx.lineWidth = px;
    const cw = 0.42, cd = 0.42;
    const side = seats >= 6 ? Math.floor((seats - 2) / 2) : seats / 2;
    const ends = seats - side * 2;
    for (let i = 0; i < side; i++) {
      const x = -w / 2 + (w * (i + 0.5)) / side - cw / 2;
      ctx.beginPath(); ctx.roundRect(x, -d / 2 - cd + 0.08, cw, cd, 0.05); ctx.fill(); ctx.stroke();
      ctx.beginPath(); ctx.roundRect(x, d / 2 - 0.08, cw, cd, 0.05); ctx.fill(); ctx.stroke();
    }
    if (ends >= 1) { ctx.beginPath(); ctx.roundRect(-w / 2 - cd + 0.08, -cw / 2, cd, cw, 0.05); ctx.fill(); ctx.stroke(); }
    if (ends >= 2) { ctx.beginPath(); ctx.roundRect(w / 2 - 0.08, -cw / 2, cd, cw, 0.05); ctx.fill(); ctx.stroke(); }
  }

  drawChairsRound(seats, rx, ry, color) {
    const ctx = this.ctx;
    ctx.fillStyle = shade(color, 0.25);
    ctx.strokeStyle = shade(color, -0.4);
    ctx.lineWidth = this.px;
    for (let i = 0; i < seats; i++) {
      const a = (i / seats) * Math.PI * 2 - Math.PI / 2;
      const cx = Math.cos(a) * (rx + 0.15), cy = Math.sin(a) * (ry + 0.15);
      ctx.save(); ctx.translate(cx, cy); ctx.rotate(a + Math.PI / 2);
      ctx.beginPath(); ctx.roundRect(-0.21, -0.21, 0.42, 0.42, 0.05); ctx.fill(); ctx.stroke();
      ctx.restore();
    }
  }

  drawDim(n) {
    const ctx = this.ctx;
    const [ax, ay] = this.w2s(n.x, n.y), [bx, by] = this.w2s(n.x2, n.y2);
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    const L = dist(n.x, n.y, n.x2, n.y2);
    const ang = Math.atan2(by - ay, bx - ax);
    ctx.strokeStyle = n.color || '#d81b60';
    ctx.fillStyle = n.color || '#d81b60';
    ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.moveTo(ax, ay); ctx.lineTo(bx, by); ctx.stroke();
    for (const [x, y, s] of [[ax, ay, 1], [bx, by, -1]]) {
      ctx.save(); ctx.translate(x, y); ctx.rotate(ang);
      ctx.beginPath(); ctx.moveTo(0, -6); ctx.lineTo(0, 6); ctx.moveTo(0, 0); ctx.lineTo(8 * s, -3); ctx.lineTo(8 * s, 3); ctx.closePath(); ctx.fill(); ctx.stroke();
      ctx.restore();
    }
    const txt = (n.text ? n.text + ': ' : '') + fmtLen(L, this.app.design.settings.units);
    ctx.save();
    ctx.translate((ax + bx) / 2, (ay + by) / 2);
    let ra = ang; if (ra > Math.PI / 2 || ra < -Math.PI / 2) ra += Math.PI;
    ctx.rotate(ra);
    ctx.font = '600 12px system-ui, sans-serif';
    ctx.textAlign = 'center';
    const tw = ctx.measureText(txt).width;
    ctx.fillStyle = 'rgba(255,255,255,0.9)';
    ctx.fillRect(-tw / 2 - 3, -17, tw + 6, 15);
    ctx.fillStyle = n.color || '#d81b60';
    ctx.fillText(txt, 0, -5);
    ctx.restore();
  }

  drawWallDim(w, units, strong) {
    const ctx = this.ctx;
    const L = wallLength(w);
    if (L * this.scale < 30) return;
    const a = wallAngle(w);
    const off = (w.thickness || 0.15) / 2 + 14 / this.scale;
    const nx = -Math.sin(a), ny = Math.cos(a);
    const [sx, sy] = this.w2s((w.x1 + w.x2) / 2 + nx * off, (w.y1 + w.y2) / 2 + ny * off);
    ctx.save();
    ctx.translate(sx, sy);
    let ra = a; if (ra > Math.PI / 2 || ra <= -Math.PI / 2) ra += Math.PI;
    ctx.rotate(ra);
    ctx.font = `${strong ? 600 : 400} 11px system-ui, sans-serif`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    const txt = fmtLen(L, units);
    const tw = ctx.measureText(txt).width;
    ctx.fillStyle = 'rgba(255,255,255,0.85)';
    ctx.fillRect(-tw / 2 - 3, -7, tw + 6, 14);
    ctx.fillStyle = strong ? '#1565c0' : '#607d8b';
    ctx.fillText(txt, 0, 0);
    ctx.restore();
    ctx.textBaseline = 'alphabetic';
  }

  drawPin(n, no) {
    const ctx = this.ctx;
    const [sx, sy] = this.w2s(n.x, n.y);
    ctx.fillStyle = n.color || '#e53935';
    ctx.strokeStyle = '#ffffff';
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.arc(sx, sy - 14, 10, Math.PI * 0.8, Math.PI * 0.2);
    ctx.lineTo(sx, sy);
    ctx.closePath();
    ctx.fill(); ctx.stroke();
    ctx.fillStyle = '#ffffff';
    ctx.font = '700 11px system-ui, sans-serif';
    ctx.textAlign = 'center';
    ctx.fillText(String(no), sx, sy - 10);
    if (n.title && this.app.design.settings.showNoteTitles !== false) {
      ctx.font = '600 11px system-ui, sans-serif';
      ctx.textAlign = 'left';
      const tw = ctx.measureText(n.title).width;
      ctx.fillStyle = 'rgba(255,255,255,0.92)';
      ctx.fillRect(sx + 12, sy - 24, tw + 8, 16);
      ctx.strokeStyle = n.color || '#e53935'; ctx.lineWidth = 1;
      ctx.strokeRect(sx + 12, sy - 24, tw + 8, 16);
      ctx.fillStyle = '#263238';
      ctx.fillText(n.title, sx + 16, sy - 12);
    }
  }

  drawLabel(n) {
    const ctx = this.ctx;
    const [sx, sy] = this.w2s(n.x, n.y);
    ctx.font = `600 ${n.size || 14}px system-ui, sans-serif`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.lineWidth = 3;
    ctx.strokeStyle = 'rgba(255,255,255,0.9)';
    ctx.strokeText(n.text || '', sx, sy);
    ctx.fillStyle = n.color || '#263238';
    ctx.fillText(n.text || '', sx, sy);
    ctx.textBaseline = 'alphabetic';
  }

  drawHandle(x, y, shape = 'sq', color = '#1e88e5') {
    const ctx = this.ctx;
    const [sx, sy] = this.w2s(x, y);
    ctx.fillStyle = '#ffffff';
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.beginPath();
    if (shape === 'circle') ctx.arc(sx, sy, HANDLE_PX - 1, 0, Math.PI * 2);
    else ctx.rect(sx - HANDLE_PX + 2, sy - HANDLE_PX + 2, (HANDLE_PX - 2) * 2, (HANDLE_PX - 2) * 2);
    ctx.fill(); ctx.stroke();
  }

  outlinePoly(pts, color, dash = false) {
    const ctx = this.ctx;
    ctx.strokeStyle = color;
    ctx.lineWidth = 2;
    ctx.setLineDash(dash ? [6, 4] : []);
    ctx.beginPath();
    pts.forEach(([x, y], i) => { const [sx, sy] = this.w2s(x, y); i ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy); });
    ctx.closePath();
    ctx.stroke();
    ctx.setLineDash([]);
  }

  drawSelection() {
    const app = this.app, ctx = this.ctx;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    const h = this.hover;
    if (h && app.tool === 'select' && !(app.selection && app.selection.id === h.id)) this.outlineObj(h.type, h.id, 'rgba(30,136,229,0.45)');
    const sel = app.selection, obj = app.getSelected();
    if (!sel || !obj) return;
    this.outlineObj(sel.type, sel.id, '#1e88e5');
    if (sel.type === 'wall') { this.drawHandle(obj.x1, obj.y1); this.drawHandle(obj.x2, obj.y2); }
    if (sel.type === 'room') obj.points.forEach(([x, y]) => this.drawHandle(x, y));
    if (sel.type === 'note' && obj.kind === 'dim') { this.drawHandle(obj.x, obj.y); this.drawHandle(obj.x2, obj.y2); }
    if (sel.type === 'item' && !obj.locked) {
      const [tx, ty] = fromItemLocal(obj, 0, -obj.d / 2);
      const [rx, ry] = fromItemLocal(obj, 0, -obj.d / 2 - ROT_HANDLE_PX / this.scale);
      const [a, b] = this.w2s(tx, ty), [c, d] = this.w2s(rx, ry);
      ctx.strokeStyle = '#1e88e5'; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(a, b); ctx.lineTo(c, d); ctx.stroke();
      this.drawHandle(rx, ry, 'circle');
      const [cx, cy] = fromItemLocal(obj, obj.w / 2, obj.d / 2);
      this.drawHandle(cx, cy);
      // size readout
      const [sx, sy] = this.w2s(...fromItemLocal(obj, 0, obj.d / 2));
      const txt = `${Math.round(obj.w * 100)} × ${Math.round(obj.d * 100)} × ${Math.round(obj.h * 100)} cm`;
      ctx.font = '11px system-ui, sans-serif'; ctx.textAlign = 'center';
      const tw = ctx.measureText(txt).width;
      ctx.fillStyle = 'rgba(21,101,192,0.9)';
      ctx.fillRect(sx - tw / 2 - 4, sy + 8, tw + 8, 16);
      ctx.fillStyle = '#fff';
      ctx.fillText(txt, sx, sy + 20);
    }
  }

  outlineObj(type, id, color) {
    const f = this.floor;
    if (type === 'item') { const it = f.items.find((i) => i.id === id); if (it) this.outlinePoly(itemCorners(it), color); }
    else if (type === 'room') { const r = f.rooms.find((i) => i.id === id); if (r) this.outlinePoly(r.points, color, true); }
    else if (type === 'wall') {
      const w = f.walls.find((i) => i.id === id); if (!w) return;
      const a = wallAngle(w), th = (w.thickness || 0.15) / 2 + 2 / this.scale;
      const nx = -Math.sin(a) * th, ny = Math.cos(a) * th;
      this.outlinePoly([[w.x1 + nx, w.y1 + ny], [w.x2 + nx, w.y2 + ny], [w.x2 - nx, w.y2 - ny], [w.x1 - nx, w.y1 - ny]], color);
    } else if (type === 'opening') {
      const o = f.openings.find((i) => i.id === id); const w = o && f.walls.find((x) => x.id === o.wallId); if (!w) return;
      const g = openingGeom(o, w), th = (w.thickness || 0.15) / 2 + 3 / this.scale, hw = o.width / 2;
      this.outlinePoly([
        [g.cx - g.ux * hw + g.nx * th, g.cy - g.uy * hw + g.ny * th], [g.cx + g.ux * hw + g.nx * th, g.cy + g.uy * hw + g.ny * th],
        [g.cx + g.ux * hw - g.nx * th, g.cy + g.uy * hw - g.ny * th], [g.cx - g.ux * hw - g.nx * th, g.cy - g.uy * hw - g.ny * th]], color);
    } else if (type === 'note') {
      const n = f.notes.find((i) => i.id === id); if (!n) return;
      const ctx = this.ctx; ctx.strokeStyle = color; ctx.lineWidth = 2;
      if (n.kind === 'dim') { const [a, b] = this.w2s(n.x, n.y), [c, d] = this.w2s(n.x2, n.y2); ctx.setLineDash([5, 4]); ctx.beginPath(); ctx.moveTo(a, b); ctx.lineTo(c, d); ctx.stroke(); ctx.setLineDash([]); }
      else { const [sx, sy] = this.w2s(n.x, n.y); ctx.beginPath(); ctx.arc(sx, sy - (n.kind === 'pin' ? 14 : 0), 15, 0, Math.PI * 2); ctx.stroke(); }
    }
  }

  drawToolPreview() {
    const app = this.app, ctx = this.ctx, m = this.mouse, e = this.lastEvent;
    ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
    if (!m.inside && !this.drag) return;
    const tool = app.tool;
    const units = app.design.settings.units;
    const cross = (x, y) => {
      const [sx, sy] = this.w2s(x, y);
      ctx.strokeStyle = '#e91e63'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(sx - 7, sy); ctx.lineTo(sx + 7, sy); ctx.moveTo(sx, sy - 7); ctx.lineTo(sx, sy + 7); ctx.stroke();
    };
    if (tool === 'wall') {
      let [x, y] = this.snap(m.x, m.y, e);
      if (this.wallStart) {
        [x, y] = this.constrain(this.wallStart[0], this.wallStart[1], x, y, e);
        const typed = this.typedWallEnd(x, y, e);
        if (typed) [x, y] = typed;
        const th = app.newWallProps().thickness;
        ctx.strokeStyle = 'rgba(69,90,100,0.6)';
        ctx.lineWidth = th * this.scale;
        ctx.lineCap = 'square';
        const [a, b] = this.w2s(...this.wallStart), [c, d] = this.w2s(x, y);
        ctx.beginPath(); ctx.moveTo(a, b); ctx.lineTo(c, d); ctx.stroke();
        ctx.lineCap = 'butt';
        this.drawWallDim({ x1: this.wallStart[0], y1: this.wallStart[1], x2: x, y2: y, thickness: th }, units, true);
      }
      cross(x, y);
    } else if (tool === 'room') {
      const [x, y] = this.snap(m.x, m.y, e);
      const pts = [...this.roomPts, [x, y]];
      if (pts.length >= 2) {
        ctx.fillStyle = 'rgba(30,136,229,0.12)';
        ctx.beginPath();
        pts.forEach(([px, py], i) => { const [sx, sy] = this.w2s(px, py); i ? ctx.lineTo(sx, sy) : ctx.moveTo(sx, sy); });
        ctx.closePath(); ctx.fill();
        this.outlinePoly(pts, '#1e88e5', true);
        if (pts.length >= 3) {
          const [cx, cy] = this.w2s(...polygonCentroid(pts));
          ctx.fillStyle = '#1565c0'; ctx.font = '600 12px system-ui'; ctx.textAlign = 'center';
          ctx.fillText(fmtArea(polygonArea(pts), units), cx, cy);
        }
      }
      this.roomPts.forEach(([px, py]) => this.drawHandle(px, py));
      cross(x, y);
    } else if ((tool === 'roombox' || tool === 'area')) {
      const [x, y] = this.snap(m.x, m.y, e);
      if (this.boxStart) {
        const [x0, y0] = this.boxStart;
        const pts = [[x0, y0], [x, y0], [x, y], [x0, y]];
        ctx.fillStyle = 'rgba(30,136,229,0.12)';
        const [a, b] = this.w2s(x0, y0), [c, d] = this.w2s(x, y);
        ctx.fillRect(Math.min(a, c), Math.min(b, d), Math.abs(c - a), Math.abs(d - b));
        this.outlinePoly(pts, '#1e88e5', tool === 'area');
        ctx.fillStyle = '#1565c0'; ctx.font = '600 12px system-ui'; ctx.textAlign = 'center';
        ctx.fillText(`${fmtLen(Math.abs(x - x0), units)} × ${fmtLen(Math.abs(y - y0), units)}  (${fmtArea(Math.abs((x - x0) * (y - y0)), units)})`, (a + c) / 2, (b + d) / 2);
      }
      cross(x, y);
    } else if (tool === 'dim') {
      let [x, y] = this.snap(m.x, m.y, e);
      if (this.dimStart) {
        [x, y] = this.constrain(this.dimStart[0], this.dimStart[1], x, y, e);
        this.drawDim({ x: this.dimStart[0], y: this.dimStart[1], x2: x, y2: y, color: '#d81b60' });
      }
      cross(x, y);
    } else if (tool === 'place' && app.placing) {
      if (app.placing.kind === 'opening') {
        const op = this.ghostOpening(m.x, m.y);
        if (op) {
          const w = this.floor.walls.find((x) => x.id === op.wallId);
          const L = wallLength(w), c = op.t * L;
          ctx.setTransform(this.dpr * this.scale, 0, 0, this.dpr * this.scale, this.dpr * this.ox, this.dpr * this.oy);
          ctx.save(); ctx.translate(w.x1, w.y1); ctx.rotate(wallAngle(w));
          ctx.globalAlpha = 0.8;
          this.drawOpening(op, w.thickness || 0.15, c - op.width / 2, c + op.width / 2);
          ctx.restore(); ctx.globalAlpha = 1;
          ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
        } else {
          ctx.fillStyle = '#c62828'; ctx.font = '12px system-ui'; ctx.textAlign = 'left';
          ctx.fillText('move onto a wall', m.sx + 12, m.sy - 8);
        }
      } else {
        const it = this.ghostItem(m.x, m.y, e);
        if (it) {
          ctx.globalAlpha = 0.6;
          this.drawItem(it);
          ctx.globalAlpha = 1;
          ctx.setTransform(this.dpr, 0, 0, this.dpr, 0, 0);
          this.outlinePoly(itemCorners(it), '#1e88e5', true);
        }
      }
    } else if (tool === 'pin' || tool === 'label') {
      cross(m.x, m.y);
    } else if (tool === 'calib') {
      if (this.calibStart) {
        const [a, b] = this.w2s(...this.calibStart);
        ctx.strokeStyle = '#ff6d00'; ctx.lineWidth = 2; ctx.setLineDash([6, 4]);
        ctx.beginPath(); ctx.moveTo(a, b); ctx.lineTo(m.sx, m.sy); ctx.stroke(); ctx.setLineDash([]);
        ctx.fillStyle = '#ff6d00'; ctx.beginPath(); ctx.arc(a, b, 5, 0, Math.PI * 2); ctx.fill();
      }
      cross(m.x, m.y);
    }
    if (app.lenBuf) {
      const txt = `⌨ ${app.lenBuf}  ↵`;
      ctx.font = '600 13px system-ui, sans-serif'; ctx.textAlign = 'left';
      const tw = ctx.measureText(txt).width;
      ctx.fillStyle = '#1565c0'; ctx.beginPath(); ctx.roundRect(m.sx + 14, m.sy + 12, tw + 14, 22, 6); ctx.fill();
      ctx.fillStyle = '#fff'; ctx.fillText(txt, m.sx + 21, m.sy + 27);
    }
  }

  getImage(src) {
    if (!src) return null;
    this._imgs = this._imgs || new Map();
    if (!this._imgs.has(src)) {
      const img = new Image();
      img.onload = () => this.requestRender();
      img.src = src;
      this._imgs.set(src, img);
    }
    return this._imgs.get(src);
  }

  // ---------- typed lengths ----------
  // direction for a typed wall: towards the mouse, snapped to 45° unless Alt is held
  typedWallEnd(mx, my, e, buf = this.app.lenBuf) {
    const L = parseLen(buf, this.app.design.settings.units);
    if (!this.wallStart || !isFinite(L) || L <= 0) return null;
    let a = Math.atan2(my - this.wallStart[1], mx - this.wallStart[0]);
    if (dist(mx, my, ...this.wallStart) < 1e-6) a = 0;
    if (!(e && e.altKey)) a = Math.round(a / (Math.PI / 4)) * (Math.PI / 4);
    return [this.wallStart[0] + Math.cos(a) * L, this.wallStart[1] + Math.sin(a) * L];
  }

  // Enter pressed with a typed value. Returns true if it was used.
  applyTyped(buf) {
    const app = this.app, units = app.design.settings.units, e = this.lastEvent;
    if (app.tool === 'wall' && this.wallStart) {
      const end = this.typedWallEnd(this.mouse.x, this.mouse.y, e, buf);
      if (!end) return false;
      app.floor.walls.push({ id: uid('w'), x1: this.wallStart[0], y1: this.wallStart[1], x2: end[0], y2: end[1], ...app.newWallProps() });
      this.wallStart = end;
      app.commit('Add wall ' + fmtLen(parseLen(buf, units), units));
      return true;
    }
    if (app.tool === 'roombox' || app.tool === 'area') {
      const m = buf.toLowerCase().split(/[x*×]/);
      if (m.length !== 2) return false;
      const W = parseLen(m[0], units), D = parseLen(m[1], units);
      if (!(W > 0 && D > 0)) return false;
      const [x0, y0] = this.boxStart || this.snap(this.mouse.x, this.mouse.y, e);
      this.boxStart = null;
      this.makeBoxRoom(x0, y0, x0 + W, y0 + D, app.tool === 'roombox');
      return true;
    }
    const sel = app.selection, obj = app.getSelected();
    if (sel && sel.type === 'wall' && obj) {
      const L = parseLen(buf, units);
      if (!(L > 0)) return false;
      app.setWallLength(obj, L);
      app.commit('Wall length ' + fmtLen(L, units), { props: true });
      return true;
    }
    if (sel && sel.type === 'note' && obj && obj.kind === 'dim') {
      const L = parseLen(buf, units);
      if (!(L > 0)) return false;
      const a = Math.atan2(obj.y2 - obj.y, obj.x2 - obj.x);
      obj.x2 = obj.x + Math.cos(a) * L; obj.y2 = obj.y + Math.sin(a) * L;
      app.commit('Dimension length', { props: true });
      return true;
    }
    return false;
  }

  // Render the whole floor to a PNG data URL (white background, fitted).
  exportPNG(scale = 2) {
    const saved = { scale: this.scale, ox: this.ox, oy: this.oy, W: this.W, H: this.H, dpr: this.dpr, cvw: this.cv.width, cvh: this.cv.height };
    const b = this.bounds() || { x0: 0, y0: 0, x1: 10, y1: 10 };
    const pxPerM = 80;
    const pad = 60;
    this.W = Math.ceil((b.x1 - b.x0) * pxPerM + pad * 2);
    this.H = Math.ceil((b.y1 - b.y0) * pxPerM + pad * 2 + 40);
    this.dpr = scale;
    this.scale = pxPerM;
    this.ox = pad - b.x0 * pxPerM;
    this.oy = pad + 40 - b.y0 * pxPerM;
    this.cv.width = this.W * scale; this.cv.height = this.H * scale;
    this.render({ exporting: true });
    const ctx = this.ctx;
    ctx.setTransform(scale, 0, 0, scale, 0, 0);
    ctx.fillStyle = '#263238'; ctx.font = '700 18px system-ui, sans-serif'; ctx.textAlign = 'left';
    ctx.fillText(`${this.app.design.name} — ${this.floor.name}`, pad, 32);
    const url = this.cv.toDataURL('image/png');
    Object.assign(this, { scale: saved.scale, ox: saved.ox, oy: saved.oy, W: saved.W, H: saved.H, dpr: saved.dpr });
    this.cv.width = saved.cvw; this.cv.height = saved.cvh;
    this.render();
    return url;
  }
}

// ---------- small helpers ----------
function clamp(v, a, b) { return Math.max(a, Math.min(b, v)); }

function ellipsize(ctx, txt, maxW) {
  if (ctx.measureText(txt).width <= maxW) return txt;
  while (txt.length > 1 && ctx.measureText(txt + '…').width > maxW) txt = txt.slice(0, -1);
  return txt + '…';
}

export function shade(hex, amt) {
  let c = (hex || '#999999').replace('#', '');
  if (c.length === 3) c = c.split('').map((x) => x + x).join('');
  const n = parseInt(c, 16);
  let r = (n >> 16) & 255, g = (n >> 8) & 255, b = n & 255;
  const t = amt < 0 ? 0 : 255, p = Math.abs(amt);
  r = Math.round((t - r) * p + r); g = Math.round((t - g) * p + g); b = Math.round((t - b) * p + b);
  return '#' + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
}
