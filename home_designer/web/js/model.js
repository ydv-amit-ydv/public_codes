// Home Designer — design data model and geometry helpers.
// Coordinates are metres on the plan: x to the right, y downwards (screen-like).
// In 3D the plan y becomes the world z axis and height is world y.

export const FORMAT_VERSION = 1;

export const uid = (p = 'id') => `${p}_${Math.random().toString(36).slice(2, 9)}${Date.now().toString(36).slice(-3)}`;

export function newFloor(name = 'Ground floor', elevation = 0, height = 3.0) {
  return { id: uid('fl'), name, elevation, height, walls: [], openings: [], rooms: [], items: [], notes: [] };
}

export function newDesign(name = 'Untitled design') {
  const now = new Date().toISOString();
  return {
    format: 'home-designer', version: FORMAT_VERSION,
    name, description: '', author: '', created: now, updated: now,
    settings: { units: 'm', grid: 0.1, snap: true, extThickness: 0.23, intThickness: 0.115, wallHeight: 3.0, wallKind: 'exterior' },
    floors: [newFloor()],
  };
}

// Bring older / partial files up to the current shape so the editor never trips on them.
export function normalizeDesign(d) {
  if (!d || typeof d !== 'object') throw new Error('Not a design file');
  const base = newDesign(d.name || 'Imported design');
  const out = { ...base, ...d, settings: { ...base.settings, ...(d.settings || {}) } };
  if (!Array.isArray(out.floors) || !out.floors.length) out.floors = [newFloor()];
  out.floors = out.floors.map((f, i) => ({
    ...newFloor(f.name || `Floor ${i}`, f.elevation ?? i * 3, f.height ?? 3),
    ...f,
    walls: (f.walls || []).map((w) => ({ ...w, kind: w.kind || ((w.thickness || 0.15) >= 0.18 ? 'exterior' : 'interior') })), openings: f.openings || [], rooms: f.rooms || [],
    items: f.items || [], notes: f.notes || [],
  }));
  out.format = 'home-designer';
  out.version = FORMAT_VERSION;
  return out;
}

// ---------- geometry ----------
export const dist = (ax, ay, bx, by) => Math.hypot(bx - ax, by - ay);
export const wallLength = (w) => dist(w.x1, w.y1, w.x2, w.y2);
export const wallAngle = (w) => Math.atan2(w.y2 - w.y1, w.x2 - w.x1);

export function projectOnSegment(px, py, x1, y1, x2, y2) {
  const dx = x2 - x1, dy = y2 - y1;
  const L2 = dx * dx + dy * dy || 1e-9;
  let t = ((px - x1) * dx + (py - y1) * dy) / L2;
  t = Math.max(0, Math.min(1, t));
  const x = x1 + t * dx, y = y1 + t * dy;
  return { t, x, y, d: dist(px, py, x, y) };
}

export function polygonArea(pts) {
  let a = 0;
  for (let i = 0; i < pts.length; i++) {
    const [x1, y1] = pts[i], [x2, y2] = pts[(i + 1) % pts.length];
    a += x1 * y2 - x2 * y1;
  }
  return Math.abs(a) / 2;
}

export function polygonPerimeter(pts) {
  let p = 0;
  for (let i = 0; i < pts.length; i++) {
    const [x1, y1] = pts[i], [x2, y2] = pts[(i + 1) % pts.length];
    p += dist(x1, y1, x2, y2);
  }
  return p;
}

export function polygonCentroid(pts) {
  let a = 0, cx = 0, cy = 0;
  for (let i = 0; i < pts.length; i++) {
    const [x1, y1] = pts[i], [x2, y2] = pts[(i + 1) % pts.length];
    const f = x1 * y2 - x2 * y1;
    a += f; cx += (x1 + x2) * f; cy += (y1 + y2) * f;
  }
  if (Math.abs(a) < 1e-9) {
    const n = pts.length || 1;
    return [pts.reduce((s, p) => s + p[0], 0) / n, pts.reduce((s, p) => s + p[1], 0) / n];
  }
  return [cx / (3 * a), cy / (3 * a)];
}

export function pointInPolygon(x, y, pts) {
  let inside = false;
  for (let i = 0, j = pts.length - 1; i < pts.length; j = i++) {
    const [xi, yi] = pts[i], [xj, yj] = pts[j];
    if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

// Item local frame: centre (x,y), width along local x, depth along local y, rotated by rot degrees.
export function toItemLocal(it, px, py) {
  const a = (-it.rot * Math.PI) / 180;
  const dx = px - it.x, dy = py - it.y;
  return [dx * Math.cos(a) - dy * Math.sin(a), dx * Math.sin(a) + dy * Math.cos(a)];
}
export function fromItemLocal(it, lx, ly) {
  const a = (it.rot * Math.PI) / 180;
  return [it.x + lx * Math.cos(a) - ly * Math.sin(a), it.y + lx * Math.sin(a) + ly * Math.cos(a)];
}
export function pointInItem(it, px, py, pad = 0) {
  const [lx, ly] = toItemLocal(it, px, py);
  return Math.abs(lx) <= it.w / 2 + pad && Math.abs(ly) <= it.d / 2 + pad;
}
export function itemCorners(it) {
  const hw = it.w / 2, hd = it.d / 2;
  return [[-hw, -hd], [hw, -hd], [hw, hd], [-hw, hd]].map(([x, y]) => fromItemLocal(it, x, y));
}

// Where an opening sits: centre point and wall direction.
export function openingGeom(op, wall) {
  const L = wallLength(wall);
  const ang = wallAngle(wall);
  const ux = Math.cos(ang), uy = Math.sin(ang);
  const cx = wall.x1 + ux * op.t * L, cy = wall.y1 + uy * op.t * L;
  return { cx, cy, ux, uy, nx: -uy, ny: ux, L, ang };
}

export function clampOpeningT(op, wall) {
  const L = wallLength(wall);
  if (L <= op.width) return 0.5;
  const m = op.width / 2 / L;
  return Math.max(m, Math.min(1 - m, op.t));
}

// Split a wall into solid pieces and opening gaps (distances in metres from x1,y1).
export function wallPieces(wall, openings) {
  const L = wallLength(wall);
  const gaps = openings
    .filter((o) => o.wallId === wall.id)
    .map((o) => {
      const c = o.t * L;
      return { op: o, s0: Math.max(0, c - o.width / 2), s1: Math.min(L, c + o.width / 2) };
    })
    .sort((a, b) => a.s0 - b.s0);
  const solids = [];
  let cur = 0;
  for (const g of gaps) {
    if (g.s0 > cur + 1e-4) solids.push([cur, g.s0]);
    cur = Math.max(cur, g.s1);
  }
  if (cur < L - 1e-4) solids.push([cur, L]);
  return { L, solids, gaps };
}

// Parse a typed length: "4.5" (current units), "450cm", "4.5m", "12'6", "12ft".
export function parseLen(str, units = 'm') {
  str = String(str).trim().toLowerCase();
  let m;
  if ((m = str.match(/^(\d+(?:\.\d+)?)'(\d+(?:\.\d+)?)?"?$/))) return (+m[1] * 12 + +(m[2] || 0)) * 0.0254;
  if ((m = str.match(/^(\d+(?:\.\d+)?)(cm|mm|m|ft|in)?$/))) {
    const v = +m[1], u = m[2] || (units === 'ft' ? 'ft' : 'm');
    return v * { cm: 0.01, mm: 0.001, m: 1, ft: 0.3048, in: 0.0254 }[u];
  }
  return NaN;
}

// ---------- formatting ----------
export function fmtLen(m, units = 'm') {
  if (units === 'ft') {
    const totalIn = m / 0.0254;
    let ft = Math.floor(totalIn / 12);
    let inch = Math.round(totalIn - ft * 12);
    if (inch === 12) { ft += 1; inch = 0; }
    return `${ft}′${inch}″`;
  }
  return `${m.toFixed(2)} m`;
}
export function fmtArea(m2, units = 'm') {
  if (units === 'ft') return `${(m2 * 10.7639).toFixed(0)} sq ft`;
  return `${m2.toFixed(2)} m²`;
}

// ---------- undo / redo ----------
export class History {
  constructor(limit = 150) { this.limit = limit; this.past = []; this.future = []; }
  reset(snapshot) { this.past = [snapshot]; this.future = []; }
  push(snapshot) {
    if (this.past[this.past.length - 1] === snapshot) return false;
    this.past.push(snapshot);
    if (this.past.length > this.limit) this.past.shift();
    this.future = [];
    return true;
  }
  undo() {
    if (this.past.length < 2) return null;
    this.future.push(this.past.pop());
    return this.past[this.past.length - 1];
  }
  redo() {
    if (!this.future.length) return null;
    const s = this.future.pop();
    this.past.push(s);
    return s;
  }
}
