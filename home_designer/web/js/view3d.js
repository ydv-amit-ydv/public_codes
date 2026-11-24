// Home Designer — 3D view.
// three.js with physically based materials, image-based lighting (RoomEnvironment),
// ground-truth ambient occlusion (GTAO), ACES filmic tone mapping, soft shadows and
// procedural textures. Rebuilds the scene from the design data whenever it changes.
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { PointerLockControls } from 'three/addons/controls/PointerLockControls.js';
import { RoomEnvironment } from 'three/addons/environments/RoomEnvironment.js';
import { RoundedBoxGeometry } from 'three/addons/geometries/RoundedBoxGeometry.js';
import { EffectComposer } from 'three/addons/postprocessing/EffectComposer.js';
import { RenderPass } from 'three/addons/postprocessing/RenderPass.js';
import { GTAOPass } from 'three/addons/postprocessing/GTAOPass.js';
import { OutputPass } from 'three/addons/postprocessing/OutputPass.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { ITEMS } from './catalog.js';
import { wallPieces, wallAngle, polygonCentroid, itemCorners } from './model.js';
import * as TX from './textures.js';

const UP = new THREE.Vector3(0, 1, 0);

export class View3D {
  constructor(container, app) {
    this.app = app;
    this.container = container;
    const r = new THREE.WebGLRenderer({ antialias: true, preserveDrawingBuffer: true, powerPreference: 'high-performance' });
    r.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    r.shadowMap.enabled = true;
    r.shadowMap.type = THREE.PCFSoftShadowMap;
    r.toneMapping = THREE.ACESFilmicToneMapping;
    r.toneMappingExposure = 1.15;
    r.outputColorSpace = THREE.SRGBColorSpace;
    container.appendChild(r.domElement);
    this.renderer = r;
    TX.setAnisotropy(Math.min(8, r.capabilities.getMaxAnisotropy()));

    this.scene = new THREE.Scene();
    this.scene.background = TX.skyTexture();
    this.scene.fog = new THREE.Fog('#dfe9ef', 70, 260);
    const pmrem = new THREE.PMREMGenerator(r);
    this.scene.environment = pmrem.fromScene(new RoomEnvironment(r), 0.04).texture;
    this.labelScene = new THREE.Scene();

    this.camera = new THREE.PerspectiveCamera(45, 1, 0.05, 600);
    this.camera.position.set(10, 14, 18);
    this.controls = new OrbitControls(this.camera, r.domElement);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.maxPolarAngle = Math.PI * 0.495;
    this.controls.screenSpacePanning = true;

    this.hemi = new THREE.HemisphereLight('#eef4ff', '#8a7c66', 0.9);
    this.scene.add(this.hemi);
    const sun = new THREE.DirectionalLight('#fff1dc', 2.6);
    sun.castShadow = true;
    sun.shadow.mapSize.set(4096, 4096);
    sun.shadow.bias = -0.0002;
    sun.shadow.normalBias = 0.02;
    sun.shadow.radius = 4;
    this.sun = sun;
    this.scene.add(sun, sun.target);

    const gtex = TX.grass().clone();
    gtex.repeat.set(200, 200);
    gtex.needsUpdate = true;
    const ground = new THREE.Mesh(new THREE.PlaneGeometry(400, 400), new THREE.MeshStandardMaterial({ map: gtex, color: '#cfd8c4', roughness: 1, envMapIntensity: 0.4 }));
    ground.rotation.x = -Math.PI / 2;
    ground.position.y = -0.03;
    ground.receiveShadow = true;
    this.ground = ground;
    this.scene.add(ground);

    this.root = new THREE.Group();
    this.scene.add(this.root);
    this.gltf = new GLTFLoader();
    this.models = new Map();
    TX.loadPhotoManifest().then(() => { this.mats.clear(); this.markDirty(); setTimeout(() => this.markDirty(), 1500); });
    this.mats = new Map();
    this.geos = new Map();
    this.opts = { allFloors: true, cutaway: false, labels: true, hq: true, sunAngle: 35 };
    this.dirty = true;
    this.visible = false;
    this.fitted = false;

    // post-processing: render -> GTAO -> tone map/sRGB output
    this.composer = new EffectComposer(r, new THREE.WebGLRenderTarget(1, 1, { type: THREE.HalfFloatType, samples: 4 }));
    this.composer.addPass(new RenderPass(this.scene, this.camera));
    this.gtao = new GTAOPass(this.scene, this.camera, 1, 1);
    this.gtao.updateGtaoMaterial({ radius: 0.6, distanceExponent: 1.5, thickness: 1.2, scale: 1.0, samples: 16, distanceFallOff: 1.0 });
    this.gtao.updatePdMaterial({ lumaPhi: 10, depthPhi: 2, normalPhi: 3, radius: 6, rings: 2, samples: 16 });
    this.gtao.blendIntensity = 0.9;
    this.composer.addPass(this.gtao);
    this.composer.addPass(new OutputPass());

    // first-person walk
    this.walk = new PointerLockControls(this.camera, r.domElement);
    this.keys = {};
    this.walk.addEventListener('unlock', () => this.exitWalk());
    window.addEventListener('keydown', (e) => { if (this.walk.isLocked) { this.keys[e.code] = true; } });
    window.addEventListener('keyup', (e) => { this.keys[e.code] = false; });
    this.clock = new THREE.Clock();

    new ResizeObserver(() => this.resize()).observe(container);
    this.resize();
    const loop = () => {
      requestAnimationFrame(loop);
      const dt = Math.min(0.05, this.clock.getDelta());
      if (!this.visible) return;
      if (this.dirty) { this.dirty = false; this.rebuild(); }
      if (this.walk.isLocked) this.walkStep(dt); else this.controls.update();
      this.draw();
    };
    loop();
  }

  draw() {
    const r = this.renderer;
    if (this.opts.hq) this.composer.render(); else r.render(this.scene, this.camera);
    if (this.opts.labels && this.labelScene.children.length) {
      r.autoClear = false;
      r.clearDepth();
      r.render(this.labelScene, this.camera);
      r.autoClear = true;
    }
  }

  setVisible(v) { this.visible = v; if (v) { this.resize(); this.dirty = true; } }
  markDirty() { this.dirty = true; }

  resize() {
    const rc = this.container.getBoundingClientRect();
    if (!rc.width || !rc.height) return;
    this.renderer.setSize(rc.width, rc.height);
    this.composer.setPixelRatio(this.renderer.getPixelRatio());
    this.composer.setSize(rc.width, rc.height);
    this.camera.aspect = rc.width / rc.height;
    this.camera.updateProjectionMatrix();
  }

  // ---------- materials ----------
  mat(kind, color = '#cccccc', extra = {}) {
    const key = kind + color + JSON.stringify(extra);
    if (this.mats.has(key)) return this.mats.get(key);
    const c = new THREE.Color(color);
    let m;
    const hsl = {}; c.getHSL(hsl);
    const ph = { wood: hsl.l < 0.22 ? 'f_wood_dark' : 'f_wood', bark: 'f_wood_dark', fabric: 'f_fabric', leather: 'f_leather', paint: 'f_plaster' }[kind];
    const P = ph && TX.photo(ph);
    if (P && P.map) {
      const ns = { fabric: 0.8, leather: 0.6, paint: 0.35, wood: 0.4, bark: 1 }[kind];
      const base = { map: P.map, normalMap: P.normalMap, normalScale: new THREE.Vector2(ns, ns), color: c };
      if (kind === 'fabric') m = new THREE.MeshPhysicalMaterial({ ...base, roughness: 0.95, sheen: 0.5, sheenRoughness: 0.6, sheenColor: new THREE.Color('#ffffff'), envMapIntensity: 0.4 });
      else if (kind === 'leather') m = new THREE.MeshStandardMaterial({ ...base, roughness: 0.5, envMapIntensity: 0.7 });
      else if (kind === 'paint') m = new THREE.MeshStandardMaterial({ ...base, roughness: 0.9, envMapIntensity: 0.75 });
      else m = new THREE.MeshStandardMaterial({ ...base, roughness: kind === 'bark' ? 1 : 0.5, envMapIntensity: 0.65 });
      if (extra.opacity !== undefined) { m.transparent = true; m.opacity = extra.opacity; }
      if (extra.side) m.side = extra.side;
      this.mats.set(key, m);
      return m;
    }
    switch (kind) {
      case 'paint': m = new THREE.MeshStandardMaterial({ color: c, map: TX.plaster(), roughness: 0.92, envMapIntensity: 0.75 }); break;
      case 'wood': m = new THREE.MeshStandardMaterial({ color: c, map: TX.woodGrain(), roughness: 0.55, envMapIntensity: 0.6 }); break;
      case 'fabric': m = new THREE.MeshPhysicalMaterial({ color: c, map: TX.fabric(), roughness: 0.95, sheen: 0.6, sheenRoughness: 0.6, sheenColor: new THREE.Color('#ffffff'), envMapIntensity: 0.35 }); break;
      case 'leather': m = new THREE.MeshStandardMaterial({ color: c, roughness: 0.45, envMapIntensity: 0.6 }); break;
      case 'metal': m = new THREE.MeshStandardMaterial({ color: c, roughness: 0.35, metalness: 0.85, envMapIntensity: 1.0 }); break;
      case 'chrome': m = new THREE.MeshStandardMaterial({ color: '#e8e8e8', roughness: 0.12, metalness: 1, envMapIntensity: 1.2 }); break;
      case 'appliance': m = new THREE.MeshStandardMaterial({ color: c, roughness: 0.28, metalness: 0.35, envMapIntensity: 0.9 }); break;
      case 'ceramic': m = new THREE.MeshPhysicalMaterial({ color: c, roughness: 0.12, clearcoat: 0.6, clearcoatRoughness: 0.1, envMapIntensity: 0.9 }); break;
      case 'glass': m = new THREE.MeshPhysicalMaterial({ color: c, roughness: 0.04, metalness: 0, transparent: true, opacity: 0.28, envMapIntensity: 1.6, depthWrite: false, side: THREE.DoubleSide }); break;
      case 'mirror': m = new THREE.MeshStandardMaterial({ color: '#dfe6ea', roughness: 0.02, metalness: 1, envMapIntensity: 1.5 }); break;
      case 'screen': m = new THREE.MeshStandardMaterial({ color: '#0b0d10', roughness: 0.15, metalness: 0.4, envMapIntensity: 1.2 }); break;
      case 'carpaint': m = new THREE.MeshPhysicalMaterial({ color: c, roughness: 0.35, metalness: 0.6, clearcoat: 1, clearcoatRoughness: 0.05, envMapIntensity: 1.2 }); break;
      case 'rubber': m = new THREE.MeshStandardMaterial({ color: '#1b1b1b', roughness: 0.9 }); break;
      case 'leaf': m = new THREE.MeshStandardMaterial({ color: c, map: TX.leaves(), roughness: 0.85, envMapIntensity: 0.4 }); break;
      case 'bark': m = new THREE.MeshStandardMaterial({ color: c, map: TX.woodGrain(), roughness: 1 }); break;
      case 'stone': m = new THREE.MeshStandardMaterial({ color: c, map: TX.marble(), roughness: 0.2, envMapIntensity: 0.8 }); break;
      case 'granite': m = new THREE.MeshStandardMaterial({ color: c, map: TX.stone('#9a9a9a'), roughness: 0.3, envMapIntensity: 0.8 }); break;
      case 'concrete': m = new THREE.MeshStandardMaterial({ color: c, map: TX.concrete('#d0d0d0'), roughness: 0.9 }); break;
      case 'water': m = new THREE.MeshPhysicalMaterial({ color: c, map: TX.water(), roughness: 0.05, transparent: true, opacity: 0.85, envMapIntensity: 1.5 }); break;
      case 'rug': m = new THREE.MeshStandardMaterial({ color: c, map: TX.rugPattern(1), roughness: 1, envMapIntensity: 0.2 }); break;
      case 'emissive': m = new THREE.MeshStandardMaterial({ color: c, emissive: c, emissiveIntensity: extra.intensity ?? 1.4, roughness: 0.4 }); break;
      case 'plastic': default: m = new THREE.MeshStandardMaterial({ color: c, roughness: 0.5, envMapIntensity: 0.6 });
    }
    if (extra.opacity !== undefined) { m.transparent = true; m.opacity = extra.opacity; }
    if (extra.side) m.side = extra.side;
    this.mats.set(key, m);
    return m;
  }

  // ---------- geometry helpers (all positioned by their base) ----------
  geo(key, make) { if (!this.geos.has(key)) this.geos.set(key, worldUV(make())); return this.geos.get(key); }

  mesh(geo, mat, x, y, z) {
    const m = new THREE.Mesh(geo, mat);
    m.position.set(x, y, z);
    m.castShadow = true; m.receiveShadow = true;
    return m;
  }

  box(w, h, d, mat, x = 0, y = 0, z = 0) {
    w = Math.max(w, 0.004); h = Math.max(h, 0.004); d = Math.max(d, 0.004);
    const g = this.geo(`b${w.toFixed(3)}_${h.toFixed(3)}_${d.toFixed(3)}`, () => new THREE.BoxGeometry(w, h, d));
    return this.mesh(g, mat, x, y + h / 2, z);
  }

  rbox(w, h, d, r, mat, x = 0, y = 0, z = 0) {
    w = Math.max(w, 0.01); h = Math.max(h, 0.01); d = Math.max(d, 0.01);
    r = Math.min(r, w / 2 - 0.001, h / 2 - 0.001, d / 2 - 0.001);
    if (r <= 0.002) return this.box(w, h, d, mat, x, y, z);
    const g = this.geo(`r${w.toFixed(3)}_${h.toFixed(3)}_${d.toFixed(3)}_${r.toFixed(3)}`, () => new RoundedBoxGeometry(w, h, d, 3, r));
    return this.mesh(g, mat, x, y + h / 2, z);
  }

  cyl(rt, rb, h, mat, x = 0, y = 0, z = 0, seg = 32) {
    h = Math.max(h, 0.004);
    const g = this.geo(`c${rt.toFixed(3)}_${rb.toFixed(3)}_${h.toFixed(3)}_${seg}`, () => new THREE.CylinderGeometry(rt, rb, h, seg));
    return this.mesh(g, mat, x, y + h / 2, z);
  }

  sphere(r, mat, x, y, z, sx = 1, sy = 1, sz = 1, detail = 24) {
    const g = this.geo(`s${r.toFixed(3)}_${detail}`, () => new THREE.SphereGeometry(r, detail, Math.round(detail * 0.7)));
    const m = this.mesh(g, mat, x, y, z);
    m.scale.set(sx, sy, sz);
    return m;
  }

  // organic blob (bush / tree crown): a few overlapping spheres
  crown(g, w, h, d, y0, mat, seed = 1) {
    const rnd = mulberry(seed);
    const n = 7;
    for (let i = 0; i < n; i++) {
      const a = (i / n) * Math.PI * 2 + rnd();
      const rr = i === 0 ? 0 : 0.28;
      const s = 0.5 + rnd() * 0.18;
      g.add(this.sphere(0.5, mat, Math.cos(a) * w * rr, y0 + h * (0.45 + rnd() * 0.15), Math.sin(a) * d * rr, w * s, h * s * 1.1, d * s, 16));
    }
  }

  clear() {
    this.root.clear();
    this.labelScene.traverse((o) => { if (o.material) { o.material.map?.dispose(); o.material.dispose(); } });
    this.labelScene.clear();
    // geometries are cached and shared; drop the cache when it gets large
    if (this.geos.size > 3000) { this.geos.forEach((g) => g.dispose()); this.geos.clear(); }
  }

  // ---------- scene ----------
  rebuild() {
    this.clear();
    const d = this.app.design;
    if (!d) return;
    const cur = this.app.floor;
    const floors = this.opts.allFloors ? d.floors.filter((f) => f.elevation <= cur.elevation + 1e-6) : [cur];
    const G = TX.photo('grass');
    if (G && G.map && G.map.image && G.normalMap?.image && !this._groundPhoto) {
      this._groundPhoto = true;
      const t = G.map.clone(); t.repeat.set(200, 200); t.needsUpdate = true;
      const n = G.normalMap.clone(); n.repeat.set(200, 200); n.needsUpdate = true;
      this.ground.material = new THREE.MeshStandardMaterial({ map: t, normalMap: n, color: '#a9c98f', roughness: 1, envMapIntensity: 0.4 });
    }
    for (const f of floors) {
      const g = new THREE.Group();
      g.position.y = f.elevation || 0;
      this.buildFloor(f, g, f === cur);
      this.root.add(g);
    }
    this.fitSun();
    if (!this.fitted) { this.fitted = true; this.fit(); }
  }

  floorMaterial(r) {
    const p = TX.pickFloor(r);
    const P = !p.plain && TX.photo(p.key);
    if (P && P.map) {
      const key = 'pfloor' + p.key;
      if (!this.mats.has(key)) this.mats.set(key, new THREE.MeshStandardMaterial({ map: P.map, normalMap: P.normalMap, normalScale: new THREE.Vector2(0.7, 0.7), roughness: p.roughness, envMapIntensity: 0.8 }));
      return this.mats.get(key);
    }
    const key = 'floor' + (p.map ? p.map.uuid : r.color) + p.roughness;
    if (this.mats.has(key)) return this.mats.get(key);
    let m;
    if (p.map) {
      const t = p.map.clone();
      t.repeat.set(1 / p.map.userData.size, 1 / p.map.userData.size);
      t.needsUpdate = true;
      m = new THREE.MeshStandardMaterial({ map: t, roughness: p.roughness, envMapIntensity: 0.8 });
    } else m = new THREE.MeshStandardMaterial({ color: r.color || '#eeeeee', roughness: p.roughness });
    this.mats.set(key, m);
    return m;
  }

  buildFloor(f, g, isCurrent) {
    const cut = this.opts.cutaway && isCurrent;
    const upper = (f.elevation || 0) > 0.5;
    // room floors (UVs are in metres, so textures tile at real size)
    for (const r of f.rooms) {
      if (r.points.length < 3) continue;
      const shape = new THREE.Shape(r.points.map(([x, y]) => new THREE.Vector2(x, -y)));
      const depth = upper ? 0.15 : 0.05;
      const geo = new THREE.ExtrudeGeometry(shape, { depth, bevelEnabled: false });
      const m = new THREE.Mesh(geo, [this.floorMaterial(r), this.mat('concrete', '#cfcfcf')]);
      m.rotation.x = -Math.PI / 2;
      m.position.y = -depth + 0.001;
      m.receiveShadow = true;
      m.castShadow = upper;
      g.add(m);
      if (this.opts.labels && isCurrent && r.name && !r.hideLabel) {
        const [cx, cy] = polygonCentroid(r.points);
        this.addLabel(r.name, cx, (f.elevation || 0) + 0.3, cy, '#263238', 'rgba(255,255,255,0.88)', 0.8);
      }
    }
    // walls
    for (const w of f.walls) this.buildWall(g, w, f, cut);
    // items
    for (const it of f.items) {
      const def = ITEMS[it.catalogId] || { model: 'box', color: '#bdbdbd' };
      if (cut && (it.elev || 0) > 1.2) continue;
      const ig = new THREE.Group();
      ig.position.set(it.x, it.elev || 0, it.y);
      ig.rotation.y = (-(it.rot || 0) * Math.PI) / 180;
      const murl = it.model3d || this.app.design.modelMap?.[it.catalogId];
      const mdl = murl ? this.getModel(murl) : null;
      try {
        if (mdl && mdl.state === 'ok') this.placeModel(ig, mdl, it);
        else this.buildItem(ig, def, it);
      } catch (e) { console.warn('item', it.catalogId, e); }
      g.add(ig);
    }
    // notes
    if (this.opts.labels && isCurrent) {
      let no = 0;
      const ey = f.elevation || 0;
      for (const n of f.notes) {
        if (n.kind === 'pin') {
          no++;
          const pm = this.mat('emissive', n.color || '#e53935', { intensity: 0.5 });
          g.add(this.cyl(0.015, 0.015, 1.7, pm, n.x, 0, n.y, 8));
          g.add(this.sphere(0.08, pm, n.x, 1.75, n.y));
          this.addLabel(`${no}. ${n.title || 'Note'}`, n.x, ey + 2.1, n.y, '#ffffff', n.color || '#e53935', 1);
        } else if (n.kind === 'label') {
          this.addLabel(n.text || '', n.x, ey + 0.6, n.y, n.color || '#263238', 'rgba(255,255,255,0.85)', 0.9);
        }
      }
    }
  }

  buildWall(g, w, f, cut) {
    const th = w.thickness || 0.15;
    const H = cut ? Math.min(1.1, w.height || f.height || 3) : w.height || f.height || 3;
    const paint = this.mat('paint', w.color || '#f1ede6');
    const capMat = this.mat('plastic', '#4a565e');
    const skirt = this.mat('wood', '#e8e2d8');
    const { L, solids, gaps } = wallPieces(w, f.openings);
    const wg = new THREE.Group();
    wg.position.set(w.x1, 0, w.y1);
    wg.rotation.y = -wallAngle(w);
    const piece = (x0, x1, y0, y1, ends) => {
      const e0 = ends && x0 <= 1e-4 ? -th / 2 : x0, e1 = ends && x1 >= L - 1e-4 ? L + th / 2 : x1;
      wg.add(this.box(e1 - e0, y1 - y0, th, paint, (e0 + e1) / 2, y0, 0));
      return [e0, e1];
    };
    for (const [s0, s1] of solids) {
      const [e0, e1] = piece(s0, s1, 0, H, true);
      if (H > 0.3 && (w.height || 3) > 1.2) wg.add(this.box(e1 - e0, 0.09, th + 0.024, skirt, (e0 + e1) / 2, 0, 0));
      if (cut || H < (w.height || 3)) wg.add(this.box(e1 - e0, 0.012, th + 0.002, capMat, (e0 + e1) / 2, H - 0.006, 0));
    }
    for (const gp of gaps) {
      const o = gp.op, def = ITEMS[o.catalogId] || {};
      const cx = (gp.s0 + gp.s1) / 2, ow = gp.s1 - gp.s0;
      const sill = def.model === 'window' ? o.sill ?? def.sill ?? 0.9 : 0;
      const top = Math.min(H, sill + (o.height || def.h || 2.1));
      if (sill > 0) piece(gp.s0, gp.s1, 0, Math.min(sill, H), false);
      if (top < H) piece(gp.s0, gp.s1, top, H, false);
      if (cut && sill >= H) continue;
      this.buildOpening(wg, o, def, cx, ow, sill, Math.min(top, H) - sill, th);
    }
    g.add(wg);
  }

  buildOpening(wg, o, def, cx, ow, sill, h, th) {
    const model = def.model;
    const frameM = this.mat('plastic', '#f4f4f2');
    if (model === 'window') {
      const fw = 0.05;
      wg.add(this.box(ow + 0.06, 0.03, th + 0.06, this.mat('granite', '#d8d4cc'), cx, sill - 0.03, 0));
      wg.add(this.box(ow, fw, th * 0.6, frameM, cx, sill, 0));
      wg.add(this.box(ow, fw, th * 0.6, frameM, cx, sill + h - fw, 0));
      wg.add(this.box(fw, h, th * 0.6, frameM, cx - ow / 2 + fw / 2, sill, 0));
      wg.add(this.box(fw, h, th * 0.6, frameM, cx + ow / 2 - fw / 2, sill, 0));
      const panes = Math.max(2, Math.round(ow / 0.7));
      for (let i = 1; i < panes; i++) wg.add(this.box(0.035, h, th * 0.5, frameM, cx - ow / 2 + (ow * i) / panes, sill, 0));
      wg.add(this.box(ow - 0.04, h - 0.04, 0.012, this.mat('glass', '#cfe8f7'), cx, sill + 0.02, 0));
    } else if (model === 'slider') {
      const glass = (def.color || '') === '#90caf9';
      const panel = glass ? this.mat('glass', '#cfe8f7') : this.mat('metal', o.color || def.color || '#78909c');
      const fm = this.mat('metal', '#5d6b73');
      wg.add(this.box(ow, 0.04, th, fm, cx, h - 0.04, 0));
      wg.add(this.box(ow, 0.02, th, fm, cx, 0, 0));
      for (const [px, pz] of [[cx - ow * 0.24, -0.025], [cx + ow * 0.24, 0.025]]) {
        wg.add(this.box(ow * 0.52, h - 0.06, 0.02, panel, px, 0.02, pz));
        wg.add(this.box(0.04, h - 0.06, 0.04, fm, px - ow * 0.26, 0.02, pz));
        wg.add(this.box(0.04, h - 0.06, 0.04, fm, px + ow * 0.26, 0.02, pz));
      }
    } else if (model === 'arch') {
      wg.add(this.box(ow, 0.06, th + 0.02, frameM, cx, h - 0.06, 0));
    } else {
      // frame + leaves open at 90° (matches the plan symbol)
      const fm = this.mat('wood', '#efe9e0');
      wg.add(this.box(0.05, h, th + 0.03, fm, cx - ow / 2 + 0.025, 0, 0));
      wg.add(this.box(0.05, h, th + 0.03, fm, cx + ow / 2 - 0.025, 0, 0));
      wg.add(this.box(ow, 0.05, th + 0.03, fm, cx, h - 0.05, 0));
      const sw = o.swing || 1;
      const leafMat = def.id === 'door_safety' ? this.mat('metal', '#37474f') : this.mat('wood', o.color || def.color || '#8d6e63');
      const leaves = model === 'door2' ? [[cx - ow / 2 + 0.05, ow / 2 - 0.05, 1], [cx + ow / 2 - 0.05, ow / 2 - 0.05, -1]] : [[o.hinge ? cx + ow / 2 - 0.05 : cx - ow / 2 + 0.05, ow - 0.1, o.hinge ? -1 : 1]];
      for (const [hx, lw, dir] of leaves) {
        const z = sw * (th / 2 + lw / 2);
        const lx = hx + dir * 0.02;
        wg.add(this.rbox(0.04, h - 0.06, lw, 0.006, leafMat, lx, 0.01, z));
        // panels + handle
        wg.add(this.box(0.046, (h - 0.06) * 0.35, lw * 0.7, this.mat('wood', shade(o.color || def.color || '#8d6e63', -0.1)), lx, 0.25, z));
        wg.add(this.box(0.046, (h - 0.06) * 0.35, lw * 0.7, this.mat('wood', shade(o.color || def.color || '#8d6e63', -0.1)), lx, 1.15, z));
        wg.add(this.box(0.09, 0.025, 0.12, this.mat('chrome'), lx, 1.0, sw * (th / 2 + lw - 0.08)));
      }
    }
  }

  // ---------- items ----------
  kindFor(def, color) {
    const c = new THREE.Color(color), hsl = {};
    c.getHSL(hsl);
    const woody = hsl.h > 0.03 && hsl.h < 0.13 && hsl.s > 0.12 && hsl.l < 0.72;
    return woody ? 'wood' : 'plastic';
  }

  buildItem(g, def, it) {
    const w = it.w, d = it.d, h = it.h, c = it.color || def.color;
    const add = (m) => (g.add(m), m);
    const body = this.mat(this.kindFor(def, c), c);
    const wood = this.mat('wood', c);
    const darkWood = this.mat('wood', shade(c, -0.3));
    const fab = this.mat('fabric', c);
    const fabDark = this.mat('fabric', shade(c, -0.18));
    const fabLight = this.mat('fabric', shade(c, 0.25));
    const legM = this.mat('wood', '#4e342e');
    const legs = (ww, dd, hh, r = 0.02, inset = 0.06, m = legM) => {
      for (const [sx, sz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]]) add(this.cyl(r, r * 0.8, hh, m, sx * (ww / 2 - inset), 0, sz * (dd / 2 - inset), 10));
    };
    const cat = def.cat || '';
    switch (def.model) {
      case 'rug': {
        if (/lawn|grass/i.test(def.name)) { add(this.box(w, 0.02, d, this.floorMaterial({ material: 'grass' }), 0, 0, 0)); break; }
        if (/path|paving|parking/i.test(def.name)) { add(this.box(w, 0.02, d, this.floorMaterial({ material: 'pavers' }), 0, 0, 0)); break; }
        if (/sump|septic|chamber|drain|pit/i.test(def.name)) { add(this.box(w, 0.025, d, this.mat('concrete', '#9e9e9e'), 0, 0, 0)); add(this.box(0.6, 0.03, 0.6, this.mat('metal', '#5d5d5d'), 0, 0, 0)); break; }
        if (/mat/i.test(def.name)) { add(this.rbox(w, 0.012, d, 0.005, this.mat('plastic', c), 0, 0, 0)); break; }
        add(this.box(w, 0.012, d, this.mat('rug', c), 0, 0, 0));
        break;
      }
      case 'rrug': add(this.cyl(w / 2, w / 2, 0.012, this.mat(/pit/i.test(def.name) ? 'concrete' : 'rug', c), 0, 0, 0, 48)).scale.z = d / w; break;
      case 'pool': {
        add(this.box(w + 0.3, 0.12, d + 0.3, this.mat('granite', '#e0ddd5'), 0, 0, 0));
        add(this.box(w, 0.02, d, this.mat('water', '#7fd3f2'), 0, 0.11, 0));
        break;
      }
      case 'shower': add(this.box(w, 0.03, d, this.mat('granite', '#cfd8dc'), 0, 0, 0)); add(this.cyl(0.05, 0.05, 0.035, this.mat('chrome'), 0, 0, 0)); break;
      case 'bed': case 'bunk': case 'crib': {
        if (def.model === 'crib') {
          add(this.box(w, 0.3, d, wood, 0, 0.15, 0));
          add(this.rbox(w - 0.06, 0.12, d - 0.06, 0.04, this.mat('fabric', '#ffffff'), 0, 0.45, 0));
          const n = Math.round(d / 0.08);
          for (let i = 0; i <= n; i++) for (const sx of [-1, 1]) add(this.box(0.02, h - 0.3, 0.02, wood, sx * (w / 2 - 0.02), 0.3, -d / 2 + (d * i) / n));
          add(this.box(w, 0.04, 0.04, wood, 0, h - 0.04, -d / 2)); add(this.box(w, 0.04, 0.04, wood, 0, h - 0.04, d / 2));
          legs(w, d, 0.15, 0.02, 0.04, wood);
          break;
        }
        const frameH = 0.32;
        if (/lounger/i.test(def.name)) { add(this.rbox(w, 0.08, d, 0.03, this.mat('fabric', c), 0, 0.3, 0)); legs(w, d, 0.3, 0.02, 0.05, this.mat('metal', '#9e9e9e')); break; }
        add(this.rbox(w, frameH - 0.08, d, 0.02, darkWood, 0, 0.08, 0));
        legs(w, d, 0.08, 0.03, 0.08);
        const mt = this.mat('fabric', '#f7f5f0');
        add(this.rbox(w - 0.05, 0.22, d - 0.08, 0.06, mt, 0, frameH, 0.02));
        add(this.rbox(w + 0.02, 0.06, d * 0.62, 0.03, fab, 0, frameH + 0.19, d * 0.19 - 0.01));
        add(this.rbox(w + 0.03, 0.02, d * 0.2, 0.01, fabDark, 0, frameH + 0.22, d * 0.32));
        const hbH = Math.max(h, frameH + 0.5);
        if (def.id !== 'divan' && def.id !== 'daybed') add(this.rbox(w + 0.04, hbH - 0.05, 0.1, 0.04, this.mat('fabric', shade(c, -0.1)), 0, 0.05, -d / 2 + 0.05));
        const n = w > 1.2 ? 2 : 1, pw = n === 2 ? (w - 0.32) / 2 : Math.min(0.7, w - 0.2);
        for (let i = 0; i < n; i++) add(this.rbox(pw, 0.14, 0.38, 0.06, this.mat('fabric', '#ffffff'), n === 2 ? -w / 2 + 0.11 + pw / 2 + i * (pw + 0.1) : 0, frameH + 0.2, -d / 2 + 0.33));
        if (def.model === 'bunk') {
          for (const [sx, sz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]]) add(this.box(0.06, h, 0.06, wood, sx * (w / 2 - 0.03), 0, sz * (d / 2 - 0.03)));
          add(this.box(w, 0.1, d, wood, 0, h * 0.58, 0));
          add(this.rbox(w - 0.06, 0.16, d - 0.08, 0.05, mt, 0, h * 0.58 + 0.1, 0));
          add(this.box(w, 0.05, 0.04, wood, 0, h - 0.05, d / 2 - 0.02));
        }
        break;
      }
      case 'sofa': {
        const leather = /recliner|theatre/i.test(def.name);
        const m = leather ? this.mat('leather', c) : fab, md = leather ? this.mat('leather', shade(c, -0.15)) : fabDark;
        const seatH = Math.min(0.42, h * 0.5), at = w > 0.6 ? Math.min(0.18, w * 0.13) : 0, bt = Math.min(0.2, d * 0.24);
        if (/bench/i.test(def.name)) {
          add(this.box(w, 0.05, d * 0.8, wood, 0, 0.42, 0.05));
          add(this.box(w, 0.4, 0.04, wood, 0, 0.45, -d / 2 + 0.05));
          legs(w, d * 0.8, 0.42, 0.025, 0.08, this.mat('metal', '#37474f'));
          break;
        }
        add(this.rbox(w, seatH - 0.1, d, 0.04, md, 0, 0.08, 0));
        legs(w, d, 0.08, 0.025, 0.08, this.mat('wood', '#3e2723'));
        add(this.rbox(w - 2 * at, h - 0.08, bt, 0.06, md, 0, 0.08, -d / 2 + bt / 2));
        if (at) for (const sx of [-1, 1]) add(this.rbox(at, h * 0.72 - 0.08, d, 0.06, md, sx * (w / 2 - at / 2), 0.08, 0));
        const seats = Math.max(1, Math.round((w - 2 * at) / 0.62));
        const sw = (w - 2 * at) / seats;
        for (let i = 0; i < seats; i++) {
          const x = -w / 2 + at + sw * (i + 0.5);
          add(this.rbox(sw - 0.02, 0.14, d - bt - 0.03, 0.06, m, x, seatH - 0.04, bt / 2 + 0.005));
          add(this.rbox(sw - 0.04, (h - seatH) * 0.95, 0.16, 0.07, m, x, seatH + 0.08, -d / 2 + bt + 0.04));
        }
        break;
      }
      case 'lsofa': {
        const sd = 0.9, bt = 0.2, seatH = 0.42;
        add(this.rbox(w, seatH - 0.1, sd, 0.04, fabDark, 0, 0.08, -d / 2 + sd / 2));
        add(this.rbox(sd, seatH - 0.1, d - sd, 0.04, fabDark, w / 2 - sd / 2, 0.08, sd / 2));
        add(this.rbox(w, h - 0.08, bt, 0.06, fabDark, 0, 0.08, -d / 2 + bt / 2));
        add(this.rbox(bt, h - 0.08, d, 0.06, fabDark, w / 2 - bt / 2, 0.08, 0));
        const n = Math.round((w - sd) / 0.62);
        for (let i = 0; i < n; i++) add(this.rbox((w - sd) / n - 0.02, 0.14, sd - bt - 0.02, 0.06, fab, -w / 2 + ((w - sd) / n) * (i + 0.5), seatH - 0.04, -d / 2 + bt + (sd - bt) / 2));
        add(this.rbox(sd - bt - 0.02, 0.14, sd - bt - 0.02, 0.06, fab, w / 2 - bt - (sd - bt) / 2, seatH - 0.04, -d / 2 + bt + (sd - bt) / 2));
        const m2 = Math.round((d - sd) / 0.62);
        for (let i = 0; i < m2; i++) add(this.rbox(sd - bt - 0.02, 0.14, (d - sd) / m2 - 0.02, 0.06, fab, w / 2 - bt - (sd - bt) / 2, seatH - 0.04, -d / 2 + sd + ((d - sd) / m2) * (i + 0.5)));
        break;
      }
      case 'table': case 'dining': case 'ldesk': {
        const t = 0.035;
        const topM = /pool|billiard/i.test(def.name) ? this.mat('fabric', '#2e7d32') : /tennis/i.test(def.name) ? this.mat('plastic', '#1565c0') : wood;
        if (def.model === 'ldesk') {
          add(this.rbox(w, t, 0.6, 0.008, topM, 0, h - t, -d / 2 + 0.3));
          add(this.rbox(0.6, t, d - 0.6, 0.008, topM, w / 2 - 0.3, h - t, 0.3));
          legs(w, d, h - t, 0.025, 0.05, this.mat('metal', '#424242'));
          break;
        }
        add(this.rbox(w, t, d, 0.008, topM, 0, h - t, 0));
        if (/pool|billiard/i.test(def.name)) add(this.box(w + 0.1, 0.12, d + 0.1, wood, 0, h - 0.12, 0));
        if (/console|dressing|desk|study|tv/i.test(def.name) || cat === 'Study & Office') add(this.box(w * 0.42, 0.14, d - 0.05, wood, w * 0.25, h - t - 0.15, 0));
        legs(w, d, h - t, 0.022, 0.05, /standing|gaming|office/i.test(def.name) ? this.mat('metal', '#37474f') : darkWood);
        if (def.model === 'dining') this.chairsRect(g, def.seats || 4, w, d, c);
        if (/dressing/i.test(def.name)) add(this.box(w * 0.7, 0.8, 0.02, this.mat('mirror'), 0, h + 0.1, -d / 2 + 0.03));
        break;
      }
      case 'rtable': case 'rdining': {
        add(this.cyl(w / 2, w / 2, 0.035, wood, 0, h - 0.035, 0, 48)).scale.z = d / w;
        add(this.cyl(0.04, 0.05, h - 0.035, darkWood, 0, 0, 0, 16));
        add(this.cyl(w * 0.22, w * 0.24, 0.025, darkWood, 0, 0, 0, 32));
        if (def.model === 'rdining') this.chairsRound(g, def.seats || 4, w / 2, d / 2, c);
        break;
      }
      case 'chair':
        if (/office/i.test(def.name)) this.officeChair(g, w, d, h, c); else this.chair(g, w, d, h, c);
        break;
      case 'stool':
        add(this.cyl(w / 2, w / 2 - 0.02, 0.06, this.mat('leather', c), 0, h - 0.06, 0));
        add(this.cyl(0.025, 0.025, h - 0.06, this.mat('chrome'), 0, 0, 0, 12));
        add(this.cyl(w * 0.4, w * 0.42, 0.02, this.mat('chrome'), 0, 0, 0));
        break;
      case 'counter': case 'ksink': case 'hob': case 'vanity': {
        const vanity = def.model === 'vanity';
        const top = def.model === 'hob' && h > 0.9 ? h - 0.05 : h - 0.04;
        const cab = vanity ? wood : this.mat(this.kindFor(def, c) === 'wood' ? 'wood' : 'plastic', c);
        add(this.box(w, 0.1, d - 0.08, this.mat('plastic', '#2b2b2b'), 0, 0, -0.04));
        add(this.box(w, top - 0.1, d - 0.03, cab, 0, 0.1, -0.015));
        const doors = Math.max(1, Math.round(w / 0.5));
        for (let i = 0; i < doors; i++) {
          const dx = -w / 2 + (w * (i + 0.5)) / doors;
          add(this.rbox(w / doors - 0.008, top - 0.14, 0.018, 0.004, cab, dx, 0.12, d / 2 - 0.025));
          add(this.box(0.12, 0.012, 0.02, this.mat('chrome'), dx, top - 0.12, d / 2 - 0.005));
        }
        add(this.box(w + 0.01, 0.035, d, vanity ? this.mat('stone', '#f5f3ef') : this.mat('granite', '#3d3a38'), 0, top, 0));
        if (def.model === 'ksink') {
          const n = w >= 1.1 ? 2 : 1;
          for (let i = 0; i < n; i++) add(this.box(Math.min(0.5, w / n - 0.12), 0.008, d * 0.6, this.mat('metal', '#cfd4d8'), -w / 2 + (w * (i + 0.5)) / n, top + 0.033, 0.02));
          add(this.cyl(0.015, 0.015, 0.3, this.mat('chrome'), 0, top + 0.035, -d / 2 + 0.08, 12));
          add(this.box(0.03, 0.03, 0.2, this.mat('chrome'), 0, top + 0.31, -d / 2 + 0.17));
        }
        if (def.model === 'hob') {
          add(this.box(Math.min(w, 0.75) * 0.95, 0.012, d * 0.85, this.mat('screen'), 0, top + 0.035, 0));
          const n = w < 0.72 && d < 0.5 ? [[-w / 4, 0], [w / 4, 0]] : [[-w / 4, -d / 5], [w / 4, -d / 5], [-w / 4, d / 5], [w / 4, d / 5]];
          for (const [x, z] of n) add(this.cyl(0.06, 0.07, 0.025, this.mat('metal', '#1d1d1d'), x, top + 0.045, z, 20));
          if (def.id === 'range') add(this.box(w * 0.8, 0.35, 0.01, this.mat('glass', '#222'), 0, 0.3, d / 2 - 0.01));
        }
        if (vanity) {
          const n = w >= 1.3 ? 2 : 1;
          for (let i = 0; i < n; i++) {
            const x = -w / 2 + (w * (i + 0.5)) / n;
            add(this.rbox(0.42, 0.12, 0.34, 0.05, this.mat('ceramic', '#ffffff'), x, top + 0.035, 0.02));
            add(this.cyl(0.012, 0.012, 0.2, this.mat('chrome'), x, top + 0.035, -d / 2 + 0.07, 10));
          }
          add(this.box(w, 0.75, 0.02, this.mat('mirror'), 0, top + 0.35, -d / 2 + 0.01));
        }
        break;
      }
      case 'basin': {
        const cer = this.mat('ceramic', '#ffffff');
        if (/utility/i.test(def.name)) { add(this.box(w, h, d, this.mat('metal', '#b0bec5'), 0, 0, 0)); break; }
        add(this.cyl(0.06, 0.08, h - 0.16, cer, 0, 0, -d / 6, 16));
        add(this.rbox(w, 0.16, d, 0.06, cer, 0, h - 0.16, 0));
        add(this.rbox(w - 0.1, 0.02, d - 0.14, 0.04, this.mat('ceramic', '#e3eaee'), 0, h - 0.01, 0.02));
        add(this.cyl(0.013, 0.013, 0.18, this.mat('chrome'), 0, h, -d / 2 + 0.06, 10));
        add(this.box(0.025, 0.025, 0.12, this.mat('chrome'), 0, h + 0.16, -d / 2 + 0.11));
        break;
      }
      case 'toilet': {
        const cer = this.mat('ceramic', c);
        if (/bidet/i.test(def.name)) { add(this.rbox(w, h, d, 0.12, cer, 0, 0, 0)); break; }
        const hung = it.catalogId === 'wc_wallhung';
        if (!hung) add(this.rbox(w * 0.95, 0.4, 0.2, 0.04, cer, 0, 0.42, -d / 2 + 0.1));
        else add(this.box(w, 0.3, 0.02, this.mat('chrome'), 0, 1.0, -d / 2 + 0.01));
        if (!hung) add(this.cyl(w * 0.28, w * 0.34, 0.4, cer, 0, 0, 0.0, 24)).scale.z = 1.35;
        add(this.cyl(w / 2, w * 0.42, 0.16, cer, 0, hung ? 0.25 : 0.26, 0.06, 28)).scale.z = 1.4;
        add(this.cyl(w / 2, w / 2, 0.03, this.mat('ceramic', '#fafafa'), 0, hung ? 0.42 : 0.42, 0.06, 28)).scale.z = 1.4;
        break;
      }
      case 'tub': {
        const cer = this.mat('ceramic', c);
        if (/jacuzzi|hot tub/i.test(def.name)) { add(this.cyl(w / 2, w / 2, h, this.mat('wood', '#8d6e63'), 0, 0, 0, 40)); add(this.cyl(w / 2 - 0.1, w / 2 - 0.1, 0.02, this.mat('water', '#6fd1f0'), 0, h - 0.12, 0, 40)); break; }
        add(this.rbox(w, h, d, 0.08, cer, 0, 0, 0));
        add(this.rbox(w - 0.14, 0.02, d - 0.14, 0.06, this.mat('water', '#9fdcf3'), 0, h - 0.12, 0));
        add(this.cyl(0.015, 0.015, 0.12, this.mat('chrome'), -w / 2 + 0.1, h, 0, 10));
        break;
      }
      case 'enclosure': {
        const gl = this.mat('glass', '#d6eef8');
        add(this.box(w, h, 0.01, gl, 0, 0.05, d / 2)); add(this.box(0.01, h, d, gl, w / 2, 0.05, 0));
        add(this.box(w, 0.05, d, this.mat('granite', '#cfd8dc'), 0, 0, 0));
        add(this.box(0.02, 0.02, d, this.mat('chrome'), w / 2, h + 0.05, 0));
        add(this.cyl(0.11, 0.11, 0.015, this.mat('chrome'), 0, 2.05, -d / 2 + 0.25, 24));
        add(this.box(0.02, 0.02, 0.25, this.mat('chrome'), 0, 2.06, -d / 2 + 0.12));
        break;
      }
      case 'stairs': {
        const n = Math.max(2, Math.round(h / 0.17)), rise = h / n, run = d / n;
        const tread = this.mat('stone', '#ece7de'), riser = this.mat('paint', '#f3f0ea');
        for (let i = 0; i < n; i++) {
          add(this.box(w, rise * (i + 1) - 0.03, run, riser, 0, 0, d / 2 - run * (i + 0.5)));
          add(this.box(w + 0.02, 0.03, run + 0.02, tread, 0, rise * (i + 1) - 0.03, d / 2 - run * (i + 0.5)));
        }
        const rail = this.mat('metal', '#5d4037');
        for (let i = 0; i < n; i += 2) add(this.cyl(0.012, 0.012, 0.9, rail, w / 2 - 0.04, rise * (i + 1), d / 2 - run * (i + 0.5), 8));
        const L = Math.hypot(d, h), hr = this.box(0.05, 0.04, L, this.mat('wood', '#5d4037'), w / 2 - 0.04, 0, 0);
        hr.position.y = h / 2 + 0.92; hr.rotation.x = Math.atan2(h, d); add(hr);
        break;
      }
      case 'spiral': {
        const n = 16, rise = h / n;
        add(this.cyl(0.07, 0.07, h + 0.9, this.mat('metal', '#455a64'), 0, 0, 0, 16));
        for (let i = 0; i < n; i++) {
          const a = (i / n) * Math.PI * 2 * 1.1;
          const pg = new THREE.Group();
          const st = this.box(w / 2, 0.04, 0.32, this.mat('wood', c), w / 4, rise * (i + 1) - 0.04, 0);
          pg.add(st);
          pg.add(this.cyl(0.01, 0.01, 0.9, this.mat('metal', '#455a64'), w / 2 - 0.03, rise * (i + 1), 0, 6));
          pg.rotation.y = a; add(pg);
        }
        break;
      }
      case 'ramp': {
        const shape = new THREE.Shape([new THREE.Vector2(-d / 2, 0), new THREE.Vector2(d / 2, 0), new THREE.Vector2(-d / 2, h)]);
        const geo = new THREE.ExtrudeGeometry(shape, { depth: w, bevelEnabled: false });
        const m = new THREE.Mesh(geo, this.mat(/slide/i.test(def.name) ? 'plastic' : 'concrete', c));
        // shape x -> world z, extrusion -> world -x; the high end points "up" the plan (-z)
        m.rotation.y = -Math.PI / 2;
        m.position.set(w / 2, 0, 0);
        m.castShadow = m.receiveShadow = true;
        add(m);
        break;
      }
      case 'plant': {
        const pot = Math.min(0.42, h * 0.3);
        add(this.cyl(w * 0.3, w * 0.22, pot, this.mat('ceramic', '#b0715a'), 0, 0, 0, 20));
        add(this.cyl(w * 0.28, w * 0.28, 0.02, this.mat('concrete', '#5d4037'), 0, pot - 0.03, 0, 20));
        if (h > 1.2) add(this.cyl(0.02, 0.03, h * 0.5, this.mat('bark', '#6d4c41'), 0, pot, 0, 8));
        this.crown(g, w, h - pot, d, pot, this.mat('leaf', c), hashStr(it.id));
        break;
      }
      case 'tree': {
        const bark = this.mat('bark', '#6d4c41');
        if (it.catalogId === 'tree_palm') {
          const segs = 6;
          for (let i = 0; i < segs; i++) add(this.cyl(0.13 - i * 0.008, 0.15 - i * 0.008, h * 0.92 / segs, bark, Math.sin(i * 0.25) * 0.15, (h * 0.92 * i) / segs, 0, 10));
          const leaf = this.mat('leaf', c, { side: THREE.DoubleSide });
          for (let i = 0; i < 9; i++) {
            const pg = new THREE.Group();
            const fr = this.box(0.35, 0.02, w * 0.55, leaf, 0, 0, w * 0.27);
            pg.add(fr); pg.position.set(Math.sin(5 * 0.25) * 0.15, h * 0.92, 0);
            pg.rotation.y = (i / 9) * Math.PI * 2; pg.rotation.x = 0.45 + (i % 2) * 0.2;
            add(pg);
          }
          break;
        }
        add(this.cyl(0.1, 0.17, h * 0.5, bark, 0, 0, 0, 12));
        this.crown(g, w, h * 0.62, d, h * 0.38, this.mat('leaf', c), hashStr(it.id));
        break;
      }
      case 'blob':
        if (cat === 'Outdoor & Garden') this.crown(g, w, h, d, -h * 0.1, this.mat('leaf', c), hashStr(it.id));
        else add(this.sphere(0.5, this.mat(cat === 'Pets' ? 'fabric' : 'leather', c), 0, h / 2, 0, w, h, d));
        break;
      case 'cone': {
        const m = new THREE.Mesh(new THREE.ConeGeometry(w / 2, h, 6, 1, true), this.mat('fabric', c, { side: THREE.DoubleSide }));
        m.position.y = h / 2; m.castShadow = true; add(m);
        break;
      }
      case 'umbrella': {
        add(this.cyl(0.025, 0.025, h, this.mat('metal', '#9e9e9e'), 0, 0, 0, 10));
        add(this.cyl(0.25, 0.3, 0.08, this.mat('concrete', '#777'), 0, 0, 0, 20));
        const m = new THREE.Mesh(new THREE.ConeGeometry(w / 2, 0.45, 16, 1, true), this.mat('fabric', c, { side: THREE.DoubleSide }));
        m.position.y = h - 0.2; m.castShadow = true; add(m);
        break;
      }
      case 'fountain': {
        const st = this.mat('stone', '#d6d2ca');
        add(this.cyl(w / 2, w / 2, 0.4, st, 0, 0, 0, 40));
        add(this.cyl(w / 2 - 0.08, w / 2 - 0.08, 0.02, this.mat('water', '#6fd1f0'), 0, 0.37, 0, 40));
        add(this.cyl(0.08, 0.12, h - 0.4, st, 0, 0.4, 0, 16));
        add(this.cyl(w * 0.25, 0.08, 0.15, st, 0, h - 0.15, 0, 24));
        break;
      }
      case 'cyl': {
        let m = body;
        if (/tank|bin|cylinder|vase|purifier|softener/i.test(def.name)) m = this.mat('plastic', c);
        if (/geyser|clock/i.test(def.name)) m = this.mat('appliance', c);
        if (/column/i.test(def.name)) m = this.mat('paint', c);
        if (/urli|bell/i.test(def.name)) m = this.mat('metal', c);
        add(this.cyl(w / 2, w / 2, h, m, 0, 0, 0, 36)).scale.z = d / w;
        break;
      }
      case 'lamp': {
        const metal = this.mat('metal', '#3a3a3a');
        add(this.cyl(Math.min(0.16, w / 2), Math.min(0.17, w / 2), 0.03, metal, 0, 0, 0, 24));
        add(this.cyl(0.012, 0.012, h * 0.8, metal, 0, 0.03, 0, 8));
        add(this.cyl(w * 0.28, w * 0.45, h * 0.22, this.mat('fabric', c, { side: THREE.DoubleSide }), 0, h * 0.78, 0, 24));
        add(this.sphere(0.04, this.mat('emissive', '#fff2c4', { intensity: 3 }), 0, h * 0.86, 0));
        break;
      }
      case 'light': {
        const long = w / d > 2.5 || d / w > 2.5;
        const em = this.mat('emissive', '#fff6dc', { intensity: 2.5 });
        if (long) add(this.box(w, Math.max(h, 0.02), d, em, 0, 0, 0));
        else { add(this.cyl(w / 2, w / 2, Math.max(h, 0.02), this.mat('plastic', '#f5f5f5'), 0, 0, 0, 24)); add(this.cyl(w / 2 - 0.02, w / 2 - 0.02, 0.005, em, 0, -0.003, 0, 24)); }
        break;
      }
      case 'pendant': {
        add(this.cyl(0.004, 0.004, Math.max(0.05, 2.95 - (it.elev || 0) - h), this.mat('metal', '#222'), 0, h, 0, 6));
        const sh = new THREE.Mesh(new THREE.ConeGeometry(w / 2, h, 24, 1, true), this.mat('metal', c, { side: THREE.DoubleSide }));
        sh.position.y = h / 2; add(sh);
        add(this.sphere(0.05, this.mat('emissive', '#fff2c4', { intensity: 4 }), 0, 0.06, 0));
        break;
      }
      case 'chandelier': {
        add(this.cyl(0.006, 0.006, Math.max(0.05, 2.95 - (it.elev || 0) - h), this.mat('metal', '#b08d57'), 0, h, 0, 6));
        const gold = this.mat('metal', '#c9a45c');
        add(this.cyl(0.04, 0.04, h * 0.7, gold, 0, h * 0.2, 0, 10));
        add(this.cyl(w / 2, w / 2 - 0.02, 0.02, gold, 0, h * 0.3, 0, 32));
        for (let i = 0; i < 8; i++) {
          const a = (i / 8) * Math.PI * 2;
          add(this.sphere(0.045, this.mat('emissive', '#fff0c0', { intensity: 3 }), Math.cos(a) * w * 0.45, h * 0.38, Math.sin(a) * w * 0.45));
        }
        break;
      }
      case 'fan': {
        const m = this.mat('metal', c);
        add(this.cyl(0.012, 0.012, Math.max(0.05, 2.95 - (it.elev || 0)), m, 0, 0.05, 0, 6));
        add(this.cyl(0.11, 0.09, 0.09, m, 0, 0, 0, 24));
        for (let i = 0; i < 3; i++) {
          const pg = new THREE.Group(); pg.rotation.y = (i / 3) * Math.PI * 2;
          const b = this.rbox(w / 2 - 0.1, 0.008, 0.11, 0.004, this.mat('wood', shade(c, 0.15)), w / 4 + 0.05, 0.03, 0);
          b.rotation.x = 0.08; pg.add(b); add(pg);
        }
        break;
      }
      case 'car': {
        const paint = this.mat('carpaint', c);
        const bodyH = h * 0.42;
        add(this.rbox(w, bodyH, d, 0.18, paint, 0, 0.28, 0));
        add(this.rbox(w * 0.88, h * 0.34, d * 0.5, 0.15, paint, 0, 0.28 + bodyH - 0.02, d * 0.02));
        add(this.rbox(w * 0.9, h * 0.28, d * 0.46, 0.12, this.mat('glass', '#20303a'), 0, 0.28 + bodyH, d * 0.02));
        for (const [sx, sz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]]) {
          const wh = this.cyl(0.32, 0.32, 0.22, this.mat('rubber'), 0, 0, 0, 24);
          wh.rotation.z = Math.PI / 2; wh.position.set(sx * (w / 2 - 0.1), 0.32, sz * (d / 2 - 0.8)); add(wh);
          const rim = this.cyl(0.18, 0.18, 0.23, this.mat('chrome'), 0, 0, 0, 16);
          rim.rotation.z = Math.PI / 2; rim.position.copy(wh.position); add(rim);
        }
        for (const sx of [-1, 1]) {
          add(this.box(0.25, 0.08, 0.03, this.mat('emissive', '#fffbe6', { intensity: 1.5 }), sx * (w / 2 - 0.25), 0.55, d / 2));
          add(this.box(0.25, 0.08, 0.03, this.mat('emissive', '#ff3b30', { intensity: 1.2 }), sx * (w / 2 - 0.25), 0.6, -d / 2));
        }
        break;
      }
      case 'washer': {
        add(this.rbox(w, h, d, 0.02, this.mat('appliance', c), 0, 0, 0));
        const ring = this.cyl(Math.min(w, h) * 0.3, Math.min(w, h) * 0.3, 0.03, this.mat('chrome'), 0, 0, 0, 32);
        ring.rotation.x = Math.PI / 2; ring.position.set(0, h * 0.45, d / 2 + 0.005); add(ring);
        const glass = this.cyl(Math.min(w, h) * 0.24, Math.min(w, h) * 0.24, 0.035, this.mat('glass', '#334'), 0, 0, 0, 32);
        glass.rotation.x = Math.PI / 2; glass.position.set(0, h * 0.45, d / 2 + 0.01); add(glass);
        add(this.box(w - 0.04, 0.1, 0.01, this.mat('screen'), 0, h - 0.13, d / 2));
        break;
      }
      case 'shelf': {
        const t = 0.022, m = this.kindFor(def, c) === 'wood' ? wood : this.mat('metal', c);
        add(this.box(t, h, d, m, -w / 2 + t / 2, 0, 0)); add(this.box(t, h, d, m, w / 2 - t / 2, 0, 0));
        add(this.box(w, h, 0.008, m, 0, 0, -d / 2 + 0.004));
        const n = Math.max(1, Math.round(h / 0.36));
        const rnd = mulberry(hashStr(it.id));
        for (let i = 0; i <= n; i++) {
          const y = Math.min(h - t, (h * i) / n);
          add(this.box(w, t, d, m, 0, y, 0));
          // books / objects on shelves
          if (i < n && /book|shelf|rack|unit/i.test(def.name) && def.cat !== 'Kitchen' && h > 0.8) {
            let x = -w / 2 + t + 0.02;
            while (x < w / 2 - 0.12) {
              const bw = 0.025 + rnd() * 0.035, bh = Math.min(0.3, h / n - 0.06) * (0.65 + rnd() * 0.35);
              if (rnd() < 0.15) { x += 0.08; continue; }
              add(this.box(bw, bh, d * 0.7, this.mat('plastic', BOOK_COLORS[Math.floor(rnd() * BOOK_COLORS.length)]), x + bw / 2, y + t, 0.02));
              x += bw + 0.004;
            }
          }
        }
        break;
      }
      case 'cabinet': {
        const m = /almirah|tool|steel/i.test(def.name) ? this.mat('appliance', c) : this.kindFor(def, c) === 'wood' ? wood : this.mat('plastic', c);
        add(this.box(w, h - 0.06, d - 0.02, m, 0, 0.06, -0.01));
        add(this.box(w - 0.02, 0.06, d - 0.06, this.mat('plastic', '#2b2b2b'), 0, 0, -0.03));
        const n = Math.max(1, Math.round(w / 0.5));
        const glassDoors = /display|crockery|showcase|bar/i.test(def.name);
        for (let i = 0; i < n; i++) {
          const x = -w / 2 + (w * (i + 0.5)) / n;
          const door = this.rbox(w / n - 0.006, h - 0.07, 0.02, 0.004, glassDoors ? this.mat('glass', '#dde8ee') : m, x, 0.065, d / 2 - 0.01);
          add(door);
          if (glassDoors) add(this.box(w / n - 0.04, (h - 0.07) * 0.55, 0.012, m, x, 0.065, d / 2 - 0.01));
          add(this.box(0.016, 0.22, 0.025, this.mat('chrome'), x + (i % 2 ? -1 : 1) * (w / n / 2 - 0.05), h * 0.48, d / 2 + 0.008));
        }
        break;
      }
      case 'mandir': {
        const gold = this.mat('metal', '#d4a537');
        add(this.box(w, h * 0.35, d, wood, 0, 0, 0));
        add(this.box(w, h * 0.5, 0.03, this.mat('wood', '#f4e3c1'), 0, h * 0.35, -d / 2 + 0.015));
        for (const sx of [-1, 1]) add(this.cyl(0.035, 0.035, h * 0.5, gold, sx * (w / 2 - 0.05), h * 0.35, d / 2 - 0.06, 12));
        add(this.box(w, 0.05, d, wood, 0, h * 0.85, 0));
        const dome = new THREE.Mesh(new THREE.ConeGeometry(w * 0.55, h * 0.13, 4), gold);
        dome.rotation.y = Math.PI / 4; dome.position.y = h * 0.9 + h * 0.065; dome.castShadow = true; add(dome);
        add(this.sphere(0.03, this.mat('emissive', '#ffb74d', { intensity: 3 }), 0, h * 0.37, 0));
        break;
      }
      case 'hood': {
        const m = this.mat('metal', '#c7ccd0');
        add(this.box(w, h * 0.22, d, m, 0, 0, 0));
        add(this.box(w * 0.3, h * 0.78, d * 0.45, m, 0, h * 0.22, -d / 2 + d * 0.225));
        add(this.box(w * 0.9, 0.005, d * 0.85, this.mat('metal', '#333'), 0, -0.002, 0));
        break;
      }
      case 'screen': {
        add(this.rbox(w, h, d * 0.4, 0.01, this.mat('plastic', '#1a1a1a'), 0, 0, -d * 0.1));
        add(this.box(w - 0.02, h - 0.02, 0.005, this.mat('screen'), 0, 0.01, d * 0.1 + 0.002));
        if (!def.wall) { add(this.box(0.08, 0.1, 0.06, this.mat('plastic', '#222'), 0, -0.1, 0)); add(this.box(0.25, 0.012, 0.18, this.mat('plastic', '#222'), 0, -0.11, 0)); }
        break;
      }
      case 'art': {
        add(this.box(w, h, d, this.mat('wood', '#3e2723'), 0, 0, 0));
        add(this.box(w - 0.08, h - 0.08, 0.005, this.mat('plastic', '#fbf7ef'), 0, 0.04, d / 2));
        const rnd = mulberry(hashStr(it.id));
        for (let i = 0; i < 4; i++) add(this.box((w - 0.2) * (0.3 + rnd() * 0.4), (h - 0.2) * (0.2 + rnd() * 0.4), 0.004, this.mat('plastic', shade(c, (rnd() - 0.5) * 0.7)), (rnd() - 0.5) * (w * 0.4), 0.1 + rnd() * (h * 0.4), d / 2 + 0.004));
        break;
      }
      case 'glass':
        if (/mirror/i.test(def.name)) { add(this.rbox(w, h, 0.02, 0.01, this.mat('metal', '#8d6e63'), 0, 0, -0.005)); add(this.box(w - 0.04, h - 0.04, 0.01, this.mat('mirror'), 0, 0.02, 0.01)); }
        else if (/aquarium/i.test(def.name)) { add(this.box(w, 0.75, d, wood, 0, 0, 0)); add(this.box(w, h - 0.75, d, this.mat('water', '#5ec8ee'), 0, 0.75, 0)); }
        else add(this.box(w, h, d, this.mat('glass', '#d6eef8'), 0, 0, 0));
        break;
      case 'curtain': {
        const n = Math.max(6, Math.round(w / 0.1));
        const cm = this.mat('fabric', c), cm2 = this.mat('fabric', shade(c, -0.12));
        for (let i = 0; i < n; i++) add(this.cyl(w / n * 0.62, w / n * 0.62, h - 0.05, i % 2 ? cm : cm2, -w / 2 + (w * (i + 0.5)) / n, 0, 0, 10));
        add(this.cyl(0.015, 0.015, w + 0.2, this.mat('chrome'), 0, 0, 0, 8)).rotation.z = Math.PI / 2;
        g.children[g.children.length - 1].position.set(0, h - 0.02, 0);
        break;
      }
      case 'railing': {
        const m = /fence|gate/i.test(def.name) && this.kindFor(def, c) === 'wood' ? wood : this.mat('metal', c);
        add(this.box(w, 0.05, Math.max(d, 0.05), m, 0, h - 0.05, 0));
        add(this.box(w, 0.04, Math.max(d, 0.04), m, 0, 0.08, 0));
        const n = Math.max(2, Math.round(w / 0.12));
        for (let i = 0; i <= n; i++) add(this.box(0.02, h - 0.05, 0.02, m, -w / 2 + (w * i) / n, 0, 0));
        break;
      }
      case 'swing': case 'pergola': {
        const p = 0.1, m = this.mat('wood', c), md = this.mat('wood', shade(c, -0.2));
        for (const [sx, sz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]]) add(this.box(p, h, p, m, sx * (w / 2 - p / 2), 0, sz * (d / 2 - p / 2)));
        if (def.model === 'pergola') {
          const n = Math.round(w / 0.35);
          for (let i = 0; i <= n; i++) add(this.box(0.05, 0.14, d + 0.4, md, -w / 2 + (w * i) / n, h, 0));
          for (const sz of [-1, 1]) add(this.box(w + 0.4, 0.16, 0.1, md, 0, h - 0.16, sz * d / 2));
          if (/gazebo/i.test(def.name)) { const roof = new THREE.Mesh(new THREE.ConeGeometry(w * 0.78, 1.0, 4), this.mat('plastic', '#8d4b37')); roof.rotation.y = Math.PI / 4; roof.position.y = h + 0.6; roof.castShadow = true; add(roof); }
        } else {
          add(this.box(w, 0.08, 0.08, md, 0, h - 0.08, 0));
          add(this.box(w * 0.7, 0.06, d * 0.5, m, 0, 0.45, 0));
          add(this.box(w * 0.7, 0.35, 0.04, m, 0, 0.51, -d * 0.24));
          for (const sx of [-1, 1]) add(this.cyl(0.008, 0.008, h - 0.53, this.mat('chrome'), sx * w * 0.33, 0.51, 0, 6));
          add(this.rbox(w * 0.66, 0.08, d * 0.45, 0.03, this.mat('fabric', '#e57373'), 0, 0.51, 0.01));
        }
        break;
      }
      case 'planter': {
        add(this.box(w, h * 0.55, d, /tulsi/i.test(def.name) ? this.mat('paint', '#d7a86e') : this.mat('concrete', c), 0, 0, 0));
        this.crown(g, w * 0.85, h * 0.6, d * 0.85, h * 0.4, this.mat('leaf', '#4caf50'), hashStr(it.id));
        break;
      }
      default: {
        // generic appliance / furniture box with a bit of detail
        let m = body;
        if (/Electrical|Laundry/.test(cat) || /fridge|refrigerator|dishwasher|microwave|ac|cooler|heater|purifier|dehumid|inverter|toaster|coffee|treadmill|bike|trainer/i.test(def.name)) m = this.mat('appliance', c);
        if (/compound|wall|parapet|column|beam|landing|chimney|fireplace|niche|duct|lift/i.test(def.name)) m = this.mat('paint', c);
        if (/hedge|flower/i.test(def.name)) { add(this.rbox(w, h * 0.25, d, 0.04, this.mat('concrete', '#8d6e63'), 0, 0, 0)); this.crown(g, w, h * 0.8, d, h * 0.2, this.mat('leaf', c), hashStr(it.id)); break; }
        add(this.rbox(w, h, d, Math.min(0.02, w / 6, h / 6, d / 6), m, 0, 0, 0));
        if (/fridge|refrigerator/i.test(def.name)) {
          add(this.box(0.004, h - 0.04, 0.01, this.mat('plastic', '#555'), w > 0.8 ? 0 : 0, 0.02, d / 2));
          add(this.box(w - 0.02, 0.004, 0.01, this.mat('plastic', '#555'), 0, h * 0.62, d / 2));
          for (const sx of w > 0.8 ? [-1, 1] : [1]) add(this.box(0.02, 0.35, 0.03, this.mat('chrome'), sx * 0.05 - (w > 0.8 ? 0 : w * 0.38), h * 0.4, d / 2 + 0.015));
        }
        if (/^AC|split ac|ac \(/i.test(def.name) || /ac_split/.test(it.catalogId)) add(this.box(w * 0.9, 0.015, 0.01, this.mat('plastic', '#bbb'), 0, h * 0.2, d / 2));
        if (/fireplace/i.test(def.name)) add(this.box(w * 0.6, h * 0.55, 0.02, this.mat('emissive', '#ff7a1a', { intensity: 1.2 }), 0, 0.1, d / 2));
      }
    }
  }

  // ---------- external glTF / GLB models ----------
  getModel(url, cb) {
    let m = this.models.get(url);
    if (!m) {
      m = { state: 'loading', waiters: [] };
      this.models.set(url, m);
      const src = this.app.assetSrc(url) || url;
      this.gltf.load(src, (g) => {
        const sc = g.scene;
        sc.traverse((o) => { if (o.isMesh) { o.castShadow = true; o.receiveShadow = true; } });
        m.scene = sc;
        m.state = 'ok';
        this.markDirty();
        m.waiters.forEach((f) => f(m));
        m.waiters = [];
      }, undefined, (err) => {
        m.state = 'error';
        console.warn('model failed', url, err);
        this.app.toast('Could not load 3D model: ' + (url.startsWith('asset:') ? 'uploaded file' : url));
      });
    }
    if (cb) { if (m.state === 'ok') cb(m); else if (m.state === 'loading') m.waiters.push(cb); }
    return m;
  }

  // natural size (metres) of a model after the item's extra model rotation
  modelSize(m, rot = 0) {
    const inner = new THREE.Group();
    inner.add(m.scene.clone(true));
    inner.rotation.y = (rot * Math.PI) / 180;
    inner.updateMatrixWorld(true);
    return new THREE.Box3().setFromObject(inner).getSize(new THREE.Vector3());
  }

  placeModel(g, m, it) {
    const inner = new THREE.Group();
    inner.add(m.scene.clone(true));
    inner.rotation.y = ((it.modelRot || 0) * Math.PI) / 180;
    inner.updateMatrixWorld(true);
    const box = new THREE.Box3().setFromObject(inner);
    const size = box.getSize(new THREE.Vector3());
    const holder = new THREE.Group();
    holder.add(inner);
    // stretch to the item's width / height / depth so plan and 3D always agree
    const sx = it.w / Math.max(size.x, 1e-4), sy = it.h / Math.max(size.y, 1e-4), sz = it.d / Math.max(size.z, 1e-4);
    holder.scale.set(sx, sy, sz);
    holder.position.set(-((box.min.x + box.max.x) / 2) * sx, -box.min.y * sy, -((box.min.z + box.max.z) / 2) * sz);
    g.add(holder);
  }

  chair(g, w, d, h, c) {
    const seat = Math.min(0.46, h * 0.5);
    const wood = this.mat('wood', c), fab = this.mat('fabric', shade(c, 0.35));
    g.add(this.rbox(w, 0.06, d, 0.015, wood, 0, seat - 0.06, 0));
    g.add(this.rbox(w - 0.04, 0.04, d - 0.06, 0.02, fab, 0, seat, 0.01));
    g.add(this.rbox(w, h - seat, 0.04, 0.012, wood, 0, seat, -d / 2 + 0.02));
    for (const [sx, sz] of [[-1, -1], [1, -1], [-1, 1], [1, 1]]) g.add(this.box(0.035, seat - 0.06, 0.035, wood, sx * (w / 2 - 0.03), 0, sz * (d / 2 - 0.03)));
  }

  officeChair(g, w, d, h, c) {
    const m = this.mat('fabric', c), metal = this.mat('chrome');
    for (let i = 0; i < 5; i++) {
      const a = (i / 5) * Math.PI * 2;
      const leg = this.box(0.03, 0.03, w * 0.45, metal, 0, 0.06, 0);
      const pg = new THREE.Group(); leg.position.z = w * 0.22; pg.add(leg);
      pg.add(this.sphere(0.03, this.mat('rubber'), 0, 0.03, w * 0.44));
      pg.rotation.y = a; g.add(pg);
    }
    g.add(this.cyl(0.025, 0.025, 0.4, metal, 0, 0.06, 0, 10));
    g.add(this.rbox(w * 0.85, 0.09, d * 0.8, 0.035, m, 0, 0.45, 0.02));
    g.add(this.rbox(w * 0.8, h - 0.58, 0.07, 0.03, m, 0, 0.58, -d / 2 + 0.08));
  }

  chairsRect(g, seats, w, d, c) {
    const side = seats >= 6 ? Math.floor((seats - 2) / 2) : seats / 2;
    const ends = seats - side * 2;
    const put = (x, z, rot) => { const cg = new THREE.Group(); cg.position.set(x, 0, z); cg.rotation.y = rot; this.chair(cg, 0.44, 0.44, 0.92, c); g.add(cg); };
    for (let i = 0; i < side; i++) {
      const x = -w / 2 + (w * (i + 0.5)) / side;
      put(x, -d / 2 - 0.1, Math.PI);
      put(x, d / 2 + 0.1, 0);
    }
    if (ends >= 1) put(-w / 2 - 0.1, 0, -Math.PI / 2);
    if (ends >= 2) put(w / 2 + 0.1, 0, Math.PI / 2);
  }

  chairsRound(g, seats, rx, rz, c) {
    for (let i = 0; i < seats; i++) {
      const a = (i / seats) * Math.PI * 2 - Math.PI / 2;
      const cg = new THREE.Group();
      cg.position.set(Math.cos(a) * (rx + 0.12), 0, Math.sin(a) * (rz + 0.12));
      cg.rotation.y = -a - Math.PI / 2;
      this.chair(cg, 0.44, 0.44, 0.92, c);
      g.add(cg);
    }
  }

  addLabel(text, x, y, z, fg, bg, scale = 1) {
    const cv = document.createElement('canvas');
    const ctx = cv.getContext('2d');
    const fs = 44;
    ctx.font = `600 ${fs}px system-ui, sans-serif`;
    const tw = Math.ceil(ctx.measureText(text).width) + 34;
    cv.width = tw; cv.height = fs + 26;
    ctx.font = `600 ${fs}px system-ui, sans-serif`;
    ctx.fillStyle = bg;
    ctx.beginPath(); ctx.roundRect(0, 0, cv.width, cv.height, 16); ctx.fill();
    ctx.fillStyle = fg; ctx.textBaseline = 'middle'; ctx.textAlign = 'center';
    ctx.fillText(text, cv.width / 2, cv.height / 2 + 2);
    const tex = new THREE.CanvasTexture(cv);
    tex.colorSpace = THREE.SRGBColorSpace;
    // constant on-screen size (not scaled by distance) so labels never cover the view
    const sp = new THREE.Sprite(new THREE.SpriteMaterial({ map: tex, depthTest: false, transparent: true, toneMapped: false, sizeAttenuation: false }));
    const k = 0.00042 * scale;
    sp.scale.set(cv.width * k, cv.height * k, 1);
    sp.position.set(x, y, z);
    this.labelScene.add(sp);
  }

  // ---------- camera / light framing ----------
  bounds() {
    const f = this.app.floor;
    let x0 = Infinity, z0 = Infinity, x1 = -Infinity, z1 = -Infinity;
    const add = (x, z) => { x0 = Math.min(x0, x); z0 = Math.min(z0, z); x1 = Math.max(x1, x); z1 = Math.max(z1, z); };
    for (const fl of this.app.design.floors) {
      fl.walls.forEach((w) => { add(w.x1, w.y1); add(w.x2, w.y2); });
      fl.rooms.forEach((r) => r.points.forEach(([x, y]) => add(x, y)));
      if (fl === f) fl.items.forEach((it) => itemCorners(it).forEach(([x, y]) => add(x, y)));
    }
    if (!isFinite(x0)) { x0 = 0; z0 = 0; x1 = 10; z1 = 8; }
    return { x0, z0, x1, z1, cx: (x0 + x1) / 2, cz: (z0 + z1) / 2, size: Math.max(x1 - x0, z1 - z0, 4) };
  }

  fitSun() {
    const b = this.bounds();
    const a = (this.opts.sunAngle * Math.PI) / 180;
    const R = b.size * 1.5 + 10;
    this.sun.position.set(b.cx + Math.cos(a) * R * 0.6, R, b.cz + Math.sin(a) * R * 0.8);
    this.sun.target.position.set(b.cx, 0, b.cz);
    const sc = this.sun.shadow.camera, s = b.size * 0.85 + 4;
    sc.left = -s; sc.right = s; sc.top = s; sc.bottom = -s; sc.near = 0.5; sc.far = R * 3;
    sc.updateProjectionMatrix();
  }

  fit() {
    const b = this.bounds();
    const y = this.app.floor.elevation || 0;
    this.controls.target.set(b.cx, y, b.cz);
    this.camera.position.set(b.cx + b.size * 0.55, y + b.size * 0.95, b.cz + b.size * 1.05);
    this.controls.update();
  }

  topView() {
    const b = this.bounds();
    this.controls.target.set(b.cx, this.app.floor.elevation || 0, b.cz);
    this.camera.position.set(b.cx, (this.app.floor.elevation || 0) + b.size * 1.45, b.cz + 0.01);
    this.controls.update();
  }

  // ---------- walk-through ----------
  enterWalk() {
    const f = this.app.floor;
    const sel = this.app.getSelected();
    let x, z;
    if (sel && this.app.selection.type === 'room') [x, z] = polygonCentroid(sel.points);
    else { const r = f.rooms.find((r) => /living|hall|lounge/i.test(r.name + r.type)) || f.rooms[0]; if (r) [x, z] = polygonCentroid(r.points); else { const b = this.bounds(); x = b.cx; z = b.cz; } }
    this.saved = { pos: this.camera.position.clone(), target: this.controls.target.clone() };
    this.camera.position.set(x, (f.elevation || 0) + 1.6, z);
    this.camera.lookAt(x + 1, (f.elevation || 0) + 1.5, z);
    this.controls.enabled = false;
    this.walk.lock();
    this.app.toast('Walk mode: W A S D to move, mouse to look, Shift to run, Esc to exit', 4000);
  }

  exitWalk() {
    this.controls.enabled = true;
    if (this.saved) { this.camera.position.copy(this.saved.pos); this.controls.target.copy(this.saved.target); this.saved = null; }
    this.controls.update();
  }

  walkStep(dt) {
    const k = this.keys, sp = (k.ShiftLeft || k.ShiftRight ? 3.2 : 1.6) * dt;
    if (k.KeyW || k.ArrowUp) this.walk.moveForward(sp);
    if (k.KeyS || k.ArrowDown) this.walk.moveForward(-sp);
    if (k.KeyA || k.ArrowLeft) this.walk.moveRight(-sp);
    if (k.KeyD || k.ArrowRight) this.walk.moveRight(sp);
  }

  // High-resolution still image of the current 3D view
  snapshot(scale = 2) {
    const r = this.renderer, pr = r.getPixelRatio();
    r.setPixelRatio(Math.max(pr, scale));
    this.resize();
    this.draw();
    const url = r.domElement.toDataURL('image/png');
    r.setPixelRatio(pr);
    this.resize();
    return url;
  }
}

// ---------- helpers ----------
const BOOK_COLORS = ['#8e2b2b', '#1f4e79', '#2e6b3a', '#c9a227', '#5b3a6e', '#d35400', '#e8e2d0', '#37474f', '#a1887f'];

function shade(hex, amt) {
  const c = new THREE.Color(hex || '#999999');
  if (amt < 0) c.multiplyScalar(1 + amt);
  else c.lerp(new THREE.Color('#ffffff'), amt);
  return '#' + c.getHexString();
}

// planar UVs in metres, picked per face by its normal, so textures keep real-world scale
function worldUV(geo) {
  const pos = geo.attributes.position, nor = geo.attributes.normal;
  if (!pos || !nor) return geo;
  const uv = new Float32Array(pos.count * 2);
  for (let i = 0; i < pos.count; i++) {
    const ax = Math.abs(nor.getX(i)), ay = Math.abs(nor.getY(i)), az = Math.abs(nor.getZ(i));
    const x = pos.getX(i), y = pos.getY(i), z = pos.getZ(i);
    if (ax >= ay && ax >= az) { uv[i * 2] = z; uv[i * 2 + 1] = y; }
    else if (ay >= az) { uv[i * 2] = x; uv[i * 2 + 1] = z; }
    else { uv[i * 2] = x; uv[i * 2 + 1] = y; }
  }
  geo.setAttribute('uv', new THREE.BufferAttribute(uv, 2));
  return geo;
}

function hashStr(s) { let h = 2166136261; for (const ch of String(s)) h = Math.imul(h ^ ch.charCodeAt(0), 16777619); return h >>> 0; }
function mulberry(a) { return () => { a |= 0; a = (a + 0x6d2b79f5) | 0; let t = Math.imul(a ^ (a >>> 15), 1 | a); t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t; return ((t ^ (t >>> 14)) >>> 0) / 4294967296; }; }
