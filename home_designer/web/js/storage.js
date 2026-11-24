// Home Designer — saving and loading.
// With the bundled server (server.py) designs are saved as JSON files in home_designer/designs/
// so they can be committed to git. Without it (static hosting) they go to this browser's storage.

const LIB_PREFIX = 'homeDesigner:lib:';
const DRAFT_KEY = 'homeDesigner:draft';

function safeLS(fn, fallback) { try { return fn(); } catch { return fallback; } }

export class Store {
  constructor() { this.server = false; }

  async init() {
    try {
      const r = await fetch('api/ping', { cache: 'no-store' });
      this.server = r.ok && (await r.json()).ok === true;
    } catch { this.server = false; }
    return this.server;
  }

  get where() { return this.server ? 'designs/ folder (server)' : 'this browser (local storage)'; }

  async list() {
    if (this.server) {
      const r = await fetch('api/designs', { cache: 'no-store' });
      if (!r.ok) throw new Error('Could not list designs: ' + r.status);
      return (await r.json()).designs;
    }
    const out = [];
    safeLS(() => {
      for (let i = 0; i < localStorage.length; i++) {
        const k = localStorage.key(i);
        if (!k.startsWith(LIB_PREFIX)) continue;
        const name = k.slice(LIB_PREFIX.length);
        let updated = '', size = 0;
        try { const v = localStorage.getItem(k); size = v.length; updated = JSON.parse(v).updated || ''; } catch { /* ignore */ }
        out.push({ name, updated, size });
      }
    });
    return out.sort((a, b) => (b.updated || '').localeCompare(a.updated || ''));
  }

  async load(name) {
    if (this.server) {
      const r = await fetch('api/designs/' + encodeURIComponent(name), { cache: 'no-store' });
      if (!r.ok) throw new Error('Could not open "' + name + '": ' + r.status);
      return r.json();
    }
    const v = safeLS(() => localStorage.getItem(LIB_PREFIX + name), null);
    if (!v) throw new Error('Design not found: ' + name);
    return JSON.parse(v);
  }

  async save(name, design) {
    const body = JSON.stringify(design, null, 1);
    if (this.server) {
      const r = await fetch('api/designs/' + encodeURIComponent(name), {
        method: 'PUT', headers: { 'Content-Type': 'application/json' }, body,
      });
      if (!r.ok) throw new Error((await r.text()) || 'Save failed: ' + r.status);
      return r.json();
    }
    safeLS(() => localStorage.setItem(LIB_PREFIX + name, body));
    return { ok: true, name };
  }

  async remove(name) {
    if (this.server) {
      const r = await fetch('api/designs/' + encodeURIComponent(name), { method: 'DELETE' });
      if (!r.ok) throw new Error('Delete failed: ' + r.status);
      return;
    }
    safeLS(() => localStorage.removeItem(LIB_PREFIX + name));
  }

  saveDraft(state) { return safeLS(() => { localStorage.setItem(DRAFT_KEY, JSON.stringify(state)); return true; }, false); }
  loadDraft() { return safeLS(() => JSON.parse(localStorage.getItem(DRAFT_KEY) || 'null'), null); }
  clearDraft() { safeLS(() => localStorage.removeItem(DRAFT_KEY)); }
}

// Turn any design name into a safe file name (same rule as server.py).
export function safeName(name) {
  return (name || 'design').trim().replace(/[^A-Za-z0-9 _.-]+/g, '-').replace(/^[.\s-]+/, '').slice(0, 80) || 'design';
}

export function download(filename, data, type = 'application/json') {
  const blob = data instanceof Blob ? data : new Blob([data], { type });
  const url = typeof data === 'string' && data.startsWith('data:') ? data : URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename;
  document.body.appendChild(a); a.click(); a.remove();
  if (!url.startsWith('data:')) setTimeout(() => URL.revokeObjectURL(url), 2000);
}
