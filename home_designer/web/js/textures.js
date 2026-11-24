// Home Designer — procedural textures (generated on canvases, no image files needed).
// Each texture knows its real-world size in metres (tex.userData.size) so floors tile at true scale.
import * as THREE from 'three';

const cache = new Map();

// ---------- real photo textures (CC0, Poly Haven) shipped in web/textures ----------
let photoInfo = null;
const photoCache = new Map();
const loader = new THREE.TextureLoader();

export async function loadPhotoManifest(base = 'textures/') {
  try {
    const r = await fetch(base + 'manifest.json', { cache: 'force-cache' });
    if (!r.ok) throw new Error(r.status);
    photoInfo = { base, list: (await r.json()).textures };
  } catch { photoInfo = { base, list: {} }; }
  return photoInfo;
}

// returns { map, normalMap, size } or null when not available (falls back to procedural)
export function photo(key) {
  if (!photoInfo || !photoInfo.list[key]) return null;
  if (photoCache.has(key)) return photoCache.get(key);
  const info = photoInfo.list[key];
  const prep = (file, srgb) => {
    if (!file) return null;
    const t = loader.load(photoInfo.base + file);
    t.wrapS = t.wrapT = THREE.RepeatWrapping;
    t.anisotropy = maxAniso;
    if (srgb) t.colorSpace = THREE.SRGBColorSpace;
    t.repeat.set(1 / info.size_m, 1 / info.size_m);
    return t;
  };
  const out = { map: prep(info.diff, true), normalMap: prep(info.nor, false), size: info.size_m };
  photoCache.set(key, out);
  return out;
}
let maxAniso = 4;
export function setAnisotropy(n) { maxAniso = n; }

// ---------- noise ----------
function rng(seed) { let s = seed >>> 0 || 1; return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 4294967296); }

function valueNoise(size, cells, seed) {
  const r = rng(seed);
  const g = Array.from({ length: (cells + 1) * (cells + 1) }, () => r());
  const at = (x, y) => g[(y % cells) * (cells + 1) + (x % cells)];
  const out = new Float32Array(size * size);
  for (let y = 0; y < size; y++) {
    const fy = (y / size) * cells, y0 = Math.floor(fy), ty = fy - y0, sy = ty * ty * (3 - 2 * ty);
    for (let x = 0; x < size; x++) {
      const fx = (x / size) * cells, x0 = Math.floor(fx), tx = fx - x0, sx = tx * tx * (3 - 2 * tx);
      const a = at(x0, y0), b = at(x0 + 1, y0), c = at(x0, y0 + 1), d = at(x0 + 1, y0 + 1);
      out[y * size + x] = a + (b - a) * sx + (c - a) * sy + (a - b - c + d) * sx * sy;
    }
  }
  return out;
}

function fbm(size, seed, octaves = 4, base = 4) {
  const out = new Float32Array(size * size);
  let amp = 1, tot = 0;
  for (let o = 0; o < octaves; o++) {
    const n = valueNoise(size, base << o, seed + o * 101);
    for (let i = 0; i < out.length; i++) out[i] += n[i] * amp;
    tot += amp; amp *= 0.5;
  }
  for (let i = 0; i < out.length; i++) out[i] /= tot;
  return out;
}

function canvas(size) { const c = document.createElement('canvas'); c.width = c.height = size; return [c, c.getContext('2d')]; }

function finish(cv, sizeM, { color = true } = {}) {
  const t = new THREE.CanvasTexture(cv);
  t.wrapS = t.wrapT = THREE.RepeatWrapping;
  t.anisotropy = maxAniso;
  if (color) t.colorSpace = THREE.SRGBColorSpace;
  t.userData.size = sizeM;
  return t;
}

function hexRGB(h) { const c = new THREE.Color(h); return [c.r * 255, c.g * 255, c.b * 255]; }

// paint per-pixel noise on top of what is already drawn: amount = +/- brightness
function grain(ctx, size, amount, seed, octaves = 4, base = 4, stretchX = 1) {
  const n = fbm(size, seed, octaves, base);
  const img = ctx.getImageData(0, 0, size, size), d = img.data;
  for (let y = 0; y < size; y++) for (let x = 0; x < size; x++) {
    const sx = Math.floor((x * stretchX) % size);
    const v = (n[y * size + sx] - 0.5) * 2 * amount;
    const i = (y * size + x) * 4;
    d[i] = clamp(d[i] + v); d[i + 1] = clamp(d[i + 1] + v); d[i + 2] = clamp(d[i + 2] + v);
  }
  ctx.putImageData(img, 0, 0);
}
const clamp = (v) => (v < 0 ? 0 : v > 255 ? 255 : v);

function memo(key, fn) { if (!cache.has(key)) cache.set(key, fn()); return cache.get(key); }

// ---------- floor finishes ----------
export function woodFloor(tone = '#c19a6b') {
  return memo('wood' + tone, () => {
    const S = 1024, [cv, ctx] = canvas(S);
    const plankW = S / 6, plankL = S / 1.5;
    const r = rng(7);
    const [br, bg, bb] = hexRGB(tone);
    for (let col = 0; col < 6; col++) {
      const off = r() * plankL;
      for (let k = -1; k < 3; k++) {
        const v = (r() - 0.5) * 36;
        ctx.fillStyle = `rgb(${clamp(br + v)},${clamp(bg + v * 0.8)},${clamp(bb + v * 0.6)})`;
        ctx.fillRect(col * plankW, off + k * plankL, plankW, plankL);
      }
    }
    // grain streaks
    const n = fbm(S, 11, 5, 2);
    const img = ctx.getImageData(0, 0, S, S), d = img.data;
    for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
      const i = (y * S + x) * 4;
      const g = Math.sin((x / plankW) * 40 + n[(y >> 3) * S + x] * 18) * 9 + (n[y * S + x] - 0.5) * 30;
      d[i] = clamp(d[i] + g); d[i + 1] = clamp(d[i + 1] + g * 0.85); d[i + 2] = clamp(d[i + 2] + g * 0.7);
    }
    ctx.putImageData(img, 0, 0);
    // joints
    ctx.fillStyle = 'rgba(60,40,25,0.55)';
    for (let col = 0; col <= 6; col++) ctx.fillRect(col * plankW - 1, 0, 2, S);
    const r2 = rng(7);
    for (let col = 0; col < 6; col++) { const off = r2() * plankL; for (let k = -1; k < 3; k++) { r2(); ctx.fillRect(col * plankW, off + k * plankL - 1, plankW, 2); } }
    return finish(cv, 1.2);
  });
}

export function tiles(tileM = 0.6, base = '#e8e2d6', grout = '#b9b2a6', veins = false) {
  return memo(`tile${tileM}${base}${grout}${veins}`, () => {
    const S = 1024, n = Math.max(1, Math.round(1.2 / tileM)), ts = S / n, [cv, ctx] = canvas(S);
    const r = rng(3);
    const [br, bg, bb] = hexRGB(base);
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
      const v = (r() - 0.5) * 14;
      ctx.fillStyle = `rgb(${clamp(br + v)},${clamp(bg + v)},${clamp(bb + v)})`;
      ctx.fillRect(x * ts, y * ts, ts, ts);
    }
    grain(ctx, S, 10, 5, 4, 8);
    if (veins) marbleVeins(ctx, S, 0.35);
    ctx.fillStyle = grout;
    const gw = Math.max(2, S * 0.004);
    for (let i = 0; i <= n; i++) { ctx.fillRect(i * ts - gw / 2, 0, gw, S); ctx.fillRect(0, i * ts - gw / 2, S, gw); }
    return finish(cv, n * tileM);
  });
}

function marbleVeins(ctx, S, alpha) {
  const n = fbm(S, 21, 5, 3);
  const img = ctx.getImageData(0, 0, S, S), d = img.data;
  for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
    const t = Math.abs(Math.sin((x + y) / S * 9 + n[y * S + x] * 14));
    const v = t < 0.06 ? (0.06 - t) / 0.06 : 0;
    const i = (y * S + x) * 4;
    const k = v * 120 * alpha * 2;
    d[i] = clamp(d[i] - k); d[i + 1] = clamp(d[i + 1] - k); d[i + 2] = clamp(d[i + 2] - k * 0.9);
  }
  ctx.putImageData(img, 0, 0);
}

export function marble(base = '#f1efea') {
  return memo('marble' + base, () => {
    const S = 1024, [cv, ctx] = canvas(S);
    ctx.fillStyle = base; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 12, 9, 5, 3);
    marbleVeins(ctx, S, 0.6);
    ctx.fillStyle = 'rgba(150,145,140,0.6)';
    ctx.fillRect(0, 0, S, 2); ctx.fillRect(0, 0, 2, S); ctx.fillRect(0, S / 2, S, 2); ctx.fillRect(S / 2, 0, 2, S);
    return finish(cv, 1.6);
  });
}

export function stone(base = '#8f9a8a') { // kota / granite
  return memo('stone' + base, () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = base; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 26, 13, 6, 4);
    ctx.fillStyle = 'rgba(40,40,40,0.35)';
    ctx.fillRect(0, 0, S, 2); ctx.fillRect(0, 0, 2, S);
    return finish(cv, 0.6);
  });
}

export function grass() {
  return memo('grass', () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#6a9a48'; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 28, 17, 6, 4);
    const r = rng(5);
    for (let i = 0; i < 9000; i++) {
      const g = 90 + r() * 90;
      ctx.fillStyle = `rgba(${g * 0.55},${g},${g * 0.35},0.5)`;
      ctx.fillRect(r() * S, r() * S, 1, 2 + r() * 3);
    }
    return finish(cv, 2);
  });
}

export function pavers(base = '#a8a29a') {
  return memo('pavers' + base, () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#6b665f'; ctx.fillRect(0, 0, S, S);
    const r = rng(8), bw = S / 4, bh = S / 8;
    const [br, bg, bb] = hexRGB(base);
    for (let y = 0; y < 8; y++) for (let x = -1; x < 5; x++) {
      const v = (r() - 0.5) * 30, ox = (y % 2) * bw / 2;
      ctx.fillStyle = `rgb(${clamp(br + v)},${clamp(bg + v)},${clamp(bb + v)})`;
      ctx.fillRect(x * bw + ox + 2, y * bh + 2, bw - 4, bh - 4);
    }
    grain(ctx, S, 14, 4, 4, 8);
    return finish(cv, 1.0);
  });
}

export function concrete(base = '#b5b5b0') {
  return memo('concrete' + base, () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = base; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 18, 31, 6, 3);
    return finish(cv, 2.5);
  });
}

export function carpet(base = '#9c8f84') {
  return memo('carpet' + base, () => {
    const S = 256, [cv, ctx] = canvas(S);
    ctx.fillStyle = base; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 20, 41, 3, 32);
    return finish(cv, 0.5);
  });
}

export function water() {
  return memo('water', () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#3fa9d8'; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 30, 51, 4, 6);
    return finish(cv, 2);
  });
}

// ---------- tintable detail maps (light greyscale, multiplied by the material colour) ----------
export function plaster() {
  return memo('plaster', () => {
    const S = 256, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#f2f2f2'; ctx.fillRect(0, 0, S, S);
    grain(ctx, S, 9, 61, 5, 4);
    return finish(cv, 1.5);
  });
}

export function woodGrain() {
  return memo('grain', () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#e9e9e9'; ctx.fillRect(0, 0, S, S);
    const n = fbm(S, 71, 4, 2);
    const img = ctx.getImageData(0, 0, S, S), d = img.data;
    for (let y = 0; y < S; y++) for (let x = 0; x < S; x++) {
      const i = (y * S + x) * 4;
      const g = Math.sin(y / S * 60 + n[y * S + x] * 12) * 14 + (n[x * S + (y % S)] - 0.5) * 18;
      d[i] = clamp(d[i] + g); d[i + 1] = clamp(d[i + 1] + g); d[i + 2] = clamp(d[i + 2] + g);
    }
    ctx.putImageData(img, 0, 0);
    return finish(cv, 1);
  });
}

export function fabric() {
  return memo('fabric', () => {
    const S = 256, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#eeeeee'; ctx.fillRect(0, 0, S, S);
    for (let y = 0; y < S; y += 2) { ctx.fillStyle = `rgba(0,0,0,${y % 4 ? 0.05 : 0.09})`; ctx.fillRect(0, y, S, 1); }
    for (let x = 0; x < S; x += 2) { ctx.fillStyle = `rgba(0,0,0,${x % 4 ? 0.04 : 0.07})`; ctx.fillRect(x, 0, 1, S); }
    grain(ctx, S, 12, 81, 3, 8);
    return finish(cv, 0.25);
  });
}

export function leaves() {
  return memo('leaves', () => {
    const S = 256, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#dcdcdc'; ctx.fillRect(0, 0, S, S);
    const r = rng(91);
    for (let i = 0; i < 1400; i++) {
      const g = 120 + r() * 135;
      ctx.fillStyle = `rgb(${g},${g},${g})`;
      ctx.beginPath(); ctx.ellipse(r() * S, r() * S, 2 + r() * 5, 1 + r() * 2.5, r() * Math.PI, 0, Math.PI * 2); ctx.fill();
    }
    return finish(cv, 0.5);
  });
}

export function rugPattern(seed = 1) {
  return memo('rug' + seed, () => {
    const S = 512, [cv, ctx] = canvas(S);
    ctx.fillStyle = '#e6e6e6'; ctx.fillRect(0, 0, S, S);
    ctx.strokeStyle = 'rgba(0,0,0,0.28)'; ctx.lineWidth = 14; ctx.strokeRect(24, 24, S - 48, S - 48);
    ctx.lineWidth = 4; ctx.strokeRect(54, 54, S - 108, S - 108);
    ctx.fillStyle = 'rgba(255,255,255,0.35)';
    ctx.save(); ctx.translate(S / 2, S / 2);
    for (let i = 0; i < 8; i++) { ctx.rotate(Math.PI / 4); ctx.beginPath(); ctx.ellipse(70, 0, 60, 20, 0, 0, Math.PI * 2); ctx.fill(); }
    ctx.fillStyle = 'rgba(0,0,0,0.18)'; ctx.beginPath(); ctx.arc(0, 0, 36, 0, Math.PI * 2); ctx.fill();
    ctx.restore();
    grain(ctx, S, 14, 101 + seed, 3, 32);
    const t = finish(cv, 1);
    t.wrapS = t.wrapT = THREE.ClampToEdgeWrapping;
    return t;
  });
}

export function skyTexture() {
  return memo('sky', () => {
    const [cv, ctx] = canvas(512);
    const g = ctx.createLinearGradient(0, 0, 0, 512);
    g.addColorStop(0, '#6fa3d8'); g.addColorStop(0.55, '#b9d5ec'); g.addColorStop(1, '#e9eef0');
    ctx.fillStyle = g; ctx.fillRect(0, 0, 512, 512);
    const t = new THREE.CanvasTexture(cv);
    t.colorSpace = THREE.SRGBColorSpace;
    return t;
  });
}

// ---------- choose a floor finish for a room ----------
export const FLOOR_MATERIALS = [
  ['auto', 'Auto (from finish / room type)'], ['wood', 'Wooden planks'], ['parquet', 'Herringbone parquet'], ['wood_dark', 'Dark wood'],
  ['terrazzo', 'Terrazzo'], ['terracotta', 'Terracotta tiles'],
  ['tiles', 'Tiles 60×60'], ['tiles_large', 'Large tiles 80×80'], ['tiles_small', 'Small tiles 30×30'],
  ['marble', 'Marble'], ['granite', 'Granite / Kota stone'], ['carpet', 'Carpet'], ['grass', 'Grass'],
  ['pavers', 'Paver blocks'], ['concrete', 'Concrete'], ['plain', 'Plain colour'],
];

const OUTDOOR_GRASS = /garden|lawn|courtyard/i;
const OUTDOOR_PAVED = /parking|garage|driveway|utility|wash/i;

export function pickFloor(room) {
  let m = room.material && room.material !== 'auto' ? room.material : null;
  const f = (room.floorFinish || '').toLowerCase(), t = room.type || '';
  if (!m) {
    if (/parquet|herringbone/.test(f)) m = 'parquet';
    else if (/terrazzo/.test(f)) m = 'terrazzo';
    else if (/terracotta|clay tile/.test(f)) m = 'terracotta';
    else if (/wood|laminate|deck|plank/.test(f)) m = 'wood';
    else if (/marble|italian/.test(f)) m = 'marble';
    else if (/granite|kota|stone|slate/.test(f)) m = 'granite';
    else if (/carpet/.test(f)) m = 'carpet';
    else if (/grass|lawn/.test(f)) m = 'grass';
    else if (/paver|brick|cobble/.test(f)) m = 'pavers';
    else if (/cement|concrete|ips/.test(f)) m = 'concrete';
    else if (/300|small|anti.?skid|mosaic/.test(f)) m = 'tiles_small';
    else if (/800|large|vitrified/.test(f)) m = 'tiles_large';
    else if (/tile/.test(f)) m = 'tiles';
    else if (OUTDOOR_GRASS.test(t)) m = 'grass';
    else if (OUTDOOR_PAVED.test(t)) m = 'pavers';
    else if (/bed|study|office|library|theatre|nursery|kids|dressing/i.test(t)) m = 'wood';
    else if (/bath|toilet|wc|powder|laundry/i.test(t)) m = 'tiles_small';
    else if (/terrace|balcony|sit-out|porch|verandah/i.test(t)) m = 'granite';
    else m = 'tiles_large';
  }
  const r = floorProcedural(m);
  r.key = m;
  return r;
}

function floorProcedural(m) {
  switch (m) {
    case 'parquet': return { map: woodFloor('#b8875a'), roughness: 0.4 };
    case 'terrazzo': return { map: tiles(0.6, '#ddd8cf', '#c9c3b8', true), roughness: 0.25 };
    case 'terracotta': return { map: tiles(0.3, '#b86f4b', '#8c5a44'), roughness: 0.7 };
    case 'wood': return { map: woodFloor('#c19a6b'), roughness: 0.45 };
    case 'wood_dark': return { map: woodFloor('#7a5236'), roughness: 0.45 };
    case 'tiles': return { map: tiles(0.6, '#e9e4da', '#bdb6aa'), roughness: 0.25 };
    case 'tiles_large': return { map: tiles(0.8, '#efebe4', '#cdc6ba', true), roughness: 0.18 };
    case 'tiles_small': return { map: tiles(0.3, '#dfe6ea', '#a9b3b8'), roughness: 0.4 };
    case 'marble': return { map: marble(), roughness: 0.12 };
    case 'granite': return { map: stone('#8e9686'), roughness: 0.5 };
    case 'carpet': return { map: carpet(), roughness: 1 };
    case 'grass': return { map: grass(), roughness: 1 };
    case 'pavers': return { map: pavers(), roughness: 0.9 };
    case 'concrete': return { map: concrete(), roughness: 0.9 };
    default: return { map: null, roughness: 0.6, plain: true };
  }
}
