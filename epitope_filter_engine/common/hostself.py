"""Host-proteome identity: does a candidate epitope look like something the vaccinated person already makes?

A T-cell epitope that is identical (or one residue away) to a human protein is a poor vaccine candidate: tolerance, or an
autoimmune risk. This module downloads UniProt's reference proteome of the host once (about 8 MB, ~20,600 proteins for
human), indexes it, and for every candidate peptide reports

    Host_Self            identical (an exact 9-mer shared with the host) | near (<= 1 mismatch in a 9-mer) | distinct
    Host_Exact_Kmers     how many of the peptide's 9-mers occur exactly in the host proteome
    Host_Best_Identity   best fraction of identical residues over any 9-mer window of the peptide against any host 9-mer
    Host_Best_Protein / Host_Best_Window / Host_Best_Start   where the closest host match is

The unit is the 9-mer (the MHC-bound core; the peptide itself for 8-mers). A peptide longer than 9 is judged by its best
9-mer window. Search: a 9-mer with at most 2 mismatches has at least one of its three 3-residue blocks unchanged, so the
host windows are looked up by block and compared in full: exact for matches with <= 2 mismatches; anything worse is
reported as "distinct" (best identity below 7/9). Peptides shorter than 9 residues are compared by brute force.

The index is kept in memory (a few hundred MB for human); it is built on first use.
"""
import gzip, io, os, threading, urllib.request

import numpy as np
import pandas as pd

from . import config

# organism -> (UniProt proteome id, taxon id, name)
ORGANISMS = {"human": ("UP000005640", "9606", "Homo sapiens"), "mouse": ("UP000000589", "10090", "Mus musculus")}
URL = ("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/Eukaryota/"
       "{up}/{up}_{tax}.fasta.gz")
MODES = ("off", "flag", "exclude_identical", "exclude_near")
AA = "ACDEFGHIKLMNPQRSTVWY"
_CODE = np.full(256, 21, np.uint8)               # 21 = unknown residue (X, U, B ...): never matches anything
for _i, _c in enumerate(AA, 1):
    _CODE[ord(_c)] = _i
BASE = 22
K = 9
SEP = 0

_LOCK = threading.RLock()
_STATE = {"downloading": False, "error": None, "bytes": 0, "total": 0, "organism": None}
_INDEX = {}                                       # organism -> Index


def organism():
    o = str((config.load().get("host_proteome") or {}).get("organism", "human")).lower()
    if o not in ORGANISMS:
        raise ValueError(f"host_proteome.organism must be one of {', '.join(ORGANISMS)}")
    return o


def path(org=None):
    """The proteome file in use: host_proteome.fasta if set (any FASTA, plain or .gz), else the downloaded one."""
    org = org or organism()
    own = (config.load().get("host_proteome") or {}).get("fasta") or ""
    if own and own != "auto":
        return os.path.expanduser(own)
    up, tax, _ = ORGANISMS[org]
    d = os.path.join(config.workspace(), "host")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{up}_{tax}.fasta.gz")


def url(org=None):
    up, tax, _ = ORGANISMS[org or organism()]
    return (config.load().get("host_proteome") or {}).get("url") or URL.format(up=up, tax=tax)


def _open(p):
    """Text handle of a FASTA, gzipped or not (decided by the file's first bytes, not its name)."""
    with open(p, "rb") as f:
        gz = f.read(2) == b"\x1f\x8b"
    return gzip.open(p, "rt", errors="replace") if gz else open(p, "rt", errors="replace")


def available(org=None):
    p = path(org)
    return os.path.isfile(p) and os.path.getsize(p) > 0


def status():
    org = organism()
    p = path(org)
    out = {"organism": org, "name": ORGANISMS[org][2], "path": p, "present": available(org), "url": url(org),
           "downloading": _STATE["downloading"], "error": _STATE["error"], "bytes": _STATE["bytes"], "total": _STATE["total"],
           "loaded": org in _INDEX and _INDEX[org].path == p}
    if out["present"]:
        out["size_mb"] = round(os.path.getsize(p) / 1e6, 1)
        if out["loaded"]:
            ix = _INDEX[org]
            out.update(proteins=len(ix.names), residues=int(ix.n_res))
    return out


# ---------------------------------------------------------------------------------------------- download
def download(org=None, progress=None):
    """Fetch the reference proteome to the workspace (written to a temporary file, then renamed; checked to be FASTA)."""
    org = org or organism()
    dest = path(org)
    own = (config.load().get("host_proteome") or {}).get("fasta") or ""
    if own and own != "auto":
        raise ValueError("host_proteome.fasta is set to your own file; nothing to download.")
    tmp = dest + ".part"
    req = urllib.request.Request(url(org), headers={"User-Agent": "epitope-filter-engine"})
    with urllib.request.urlopen(req, timeout=60) as r, open(tmp, "wb") as f:
        total = int(r.headers.get("Content-Length") or 0)
        got = 0
        _STATE.update(bytes=0, total=total)
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            got += len(chunk)
            _STATE["bytes"] = got
            if progress:
                progress(got, total)
    try:
        n = sum(1 for line in _open(tmp) if line.startswith(">"))
    except (OSError, EOFError, UnicodeError) as e:
        os.remove(tmp)
        raise ValueError(f"the download is not a readable FASTA ({e})")
    if n < 1:
        os.remove(tmp)
        raise ValueError("the download contains no FASTA sequences")
    os.replace(tmp, dest)
    with _LOCK:
        _INDEX.pop(org, None)
    return n


def start_download(org=None):
    """Download in the background; status() shows progress. False when one is already running."""
    with _LOCK:
        if _STATE["downloading"]:
            return False
        _STATE.update(downloading=True, error=None, bytes=0, total=0, organism=org or organism())

    def work():
        try:
            download(org)
        except Exception as e:                              # noqa: BLE001 - shown in Settings
            _STATE["error"] = f"{type(e).__name__}: {e}"
        finally:
            _STATE["downloading"] = False
    threading.Thread(target=work, daemon=True).start()
    return True


# ---------------------------------------------------------------------------------------------- index
class Index:
    def __init__(self, p):
        self.path = p
        self.names, parts, starts = [], [], []
        pos = 0
        seq, name = [], None

        def flush():
            nonlocal pos
            if name is not None and seq:
                s = "".join(seq)
                self.names.append(name)
                starts.append(pos)
                parts.append(s)
                pos += len(s) + 1
        with _open(p) as f:
            for line in f:
                if line.startswith(">"):
                    flush()
                    seq, name = [], line[1:].split(None, 1)[0]
                else:
                    seq.append(line.strip().upper())
            flush()
        if not parts:
            raise ValueError(f"no sequences in {p}")
        self.text = "*".join(parts)                                  # for whole-peptide search
        self.prot_start = np.array(starts, np.int64)
        raw = np.frombuffer(self.text.encode("ascii", "replace"), np.uint8)
        self.seq = _CODE[raw].copy()
        self.seq[raw == ord("*")] = SEP
        self.n_res = int((self.seq != SEP).sum())
        n = len(self.seq)
        # a 9-window is usable when none of its residues is a separator or unknown
        ok = (self.seq != SEP) & (self.seq != 21)
        c = np.concatenate([[0], np.cumsum(~ok)])
        self.nwin = max(0, n - K + 1)
        self.valid = (c[K:] - c[:-K]) == 0 if self.nwin else np.zeros(0, bool)
        self.block = []
        for o in (0, 3, 6):
            if self.nwin:
                v = (self.seq[o:o + self.nwin].astype(np.int32) * BASE + self.seq[o + 1:o + 1 + self.nwin]) * BASE + self.seq[o + 2:o + 2 + self.nwin]
                st = np.flatnonzero(self.valid).astype(np.int32)
                vv = v[st].astype(np.int16)
                order = np.argsort(vv, kind="stable")
                self.block.append((vv[order], st[order]))
            else:
                self.block.append((np.zeros(0, np.int16), np.zeros(0, np.int32)))
        self.cache = {}

    def locate(self, i):
        """-> (protein name, 1-based position) of the residue index i."""
        k = int(np.searchsorted(self.prot_start, i, side="right")) - 1
        return self.names[k], int(i - self.prot_start[k] + 1)

    def _encode(self, s):
        return _CODE[np.frombuffer(s.encode("ascii", "replace"), np.uint8)]

    def best9(self, kmer):
        """-> (identity 0-1 or nan, window start or -1). Exact for <= 2 mismatches."""
        hit = self.cache.get(kmer)
        if hit is not None:
            return hit
        q = self._encode(kmer)
        cand = []
        for (vals, starts), o in zip(self.block, (0, 3, 6)):
            v = np.int16((int(q[o]) * BASE + int(q[o + 1])) * BASE + int(q[o + 2]))
            lo, hi = np.searchsorted(vals, v, "left"), np.searchsorted(vals, v, "right")
            if hi > lo:
                cand.append(starts[lo:hi].astype(np.int64))
        res = (float("nan"), -1)
        if cand:
            c = np.unique(np.concatenate(cand))
            c = c[(c >= 0) & (c < self.nwin)]
            c = c[self.valid[c]]
            if len(c):
                same = (self.seq[c[:, None] + np.arange(K)] == q).sum(axis=1)
                j = int(np.argmax(same))
                res = (float(same[j]) / K, int(c[j]))
        self.cache[kmer] = res
        return res

    def best_short(self, kmer):
        """Brute force for peptides shorter than 9: best identity over every window of the same length."""
        hit = self.cache.get(kmer)
        if hit is not None:
            return hit
        L = len(kmer)
        q = self._encode(kmer)
        best, at = -1, -1
        step = 2_000_000
        for a in range(0, len(self.seq) - L + 1, step):
            seg = self.seq[a:a + step + L - 1]
            if len(seg) < L:
                break
            w = np.lib.stride_tricks.sliding_window_view(seg, L)
            ok = (w != SEP).all(axis=1)
            sc = np.where(ok, (w == q).sum(axis=1), -1)
            j = int(np.argmax(sc))
            if sc[j] > best:
                best, at = int(sc[j]), a + j
        res = (best / L, at) if best >= 0 else (float("nan"), -1)
        self.cache[kmer] = res
        return res

    def has_full(self, pep):
        return pep in self.text


def index(org=None):
    org = org or organism()
    p = path(org)
    with _LOCK:
        ix = _INDEX.get(org)
        if ix is None or ix.path != p or ix.mtime != os.path.getmtime(p):
            ix = Index(p)
            ix.mtime = os.path.getmtime(p)
            _INDEX[org] = ix
        return ix


def warm():
    """Build the index in the background so the first run does not wait for it."""
    if available():
        threading.Thread(target=lambda: index(), daemon=True).start()


# ---------------------------------------------------------------------------------------------- per peptide
def check(pep, ix, near_min=8 / 9):
    """-> dict of Host_* values for one peptide."""
    L = len(pep)
    if L < 5 or any(ch not in AA for ch in pep):
        return None
    wins = [pep[i:i + K] for i in range(L - K + 1)] if L >= K else [pep]
    best, at, win, exact = -1.0, -1, "", 0
    for w in wins:
        ident, st = (ix.best9(w) if len(w) == K else ix.best_short(w))
        if ident == ident and ident >= 1.0 - 1e-9:
            exact += 1
        if ident == ident and ident > best:
            best, at, win = ident, st, w
    if exact:
        self_ = "identical"
    elif best >= near_min - 1e-9:
        self_ = "near"
    else:
        self_ = "distinct"
    prot, ppos = ("", None)
    hostwin = ""
    if at >= 0:
        prot, ppos = ix.locate(at)
        hostwin = ix.text[at:at + len(win)]
    return {"Host_Self": self_, "Host_Exact_Kmers": exact, "Host_Best_Identity": round(best, 3) if best >= 0 else np.nan,
            "Host_Best_Protein": prot, "Host_Best_Window": hostwin, "Host_Best_Start": ppos}


COLUMNS = ["Host_Self", "Host_Exact_Kmers", "Host_Best_Identity", "Host_Best_Protein", "Host_Best_Window", "Host_Best_Start"]


def annotate(pep, rules, notes=None):
    """Adds the Host_* columns to the peptide table (candidates only) and, in an exclude mode, removes flagged peptides from
    `Candidate` (column Host_Excluded says which). Nothing happens when the proteome is missing; a note says so."""
    mode = rules.get("host_mode", "flag")
    if mode == "off" or pep is None or pep.empty:
        return pep
    for c in COLUMNS + ["Host_Excluded"]:
        pep[c] = np.nan if c in ("Host_Exact_Kmers", "Host_Best_Identity", "Host_Best_Start") else ("" if c != "Host_Excluded" else False)
    if not available():
        if notes is not None:
            notes.append(f"Host proteome ({organism()}) is not downloaded, so self-similarity was not checked. "
                         "Settings > Host proteome > Download (about 8 MB, once).")
        return pep
    ix = index()
    near_min = float(rules.get("host_near_min", 8 / 9))
    cache = {}
    sel = pep.index[pep["Candidate"].astype(bool)]
    for i in sel:
        p = pep.at[i, "Peptide"]
        r = cache.get(p)
        if r is None:
            r = cache[p] = check(p, ix, near_min) or {}
        for k, v in r.items():
            pep.at[i, k] = v
    pep["Host_Excluded"] = False
    if notes is not None and len(sel):
        vc = pep.loc[sel, "Host_Self"].value_counts()
        notes.append(f"Host check ({ORGANISMS[organism()][2]}): of {len(sel)} candidate peptide(s), {int(vc.get('identical', 0))} identical, "
                     f"{int(vc.get('near', 0))} near (1 mismatch), {int(vc.get('distinct', 0))} distinct. Columns Host_Self, Host_Best_Identity, Host_Best_Protein.")
    if mode in ("exclude_identical", "exclude_near"):
        bad = {"exclude_identical": ("identical",), "exclude_near": ("identical", "near")}[mode]
        flagged = pep["Host_Self"].isin(bad)
        pep.loc[flagged, "Host_Excluded"] = True
        pep["Candidate"] = pep["Candidate"].astype(bool) & ~flagged
        if notes is not None and flagged.any():
            notes.append(f"{int(flagged.sum())} candidate peptide(s) were removed because they match the {ORGANISMS[organism()][2]} proteome "
                         f"({'identical 9-mer' if mode == 'exclude_identical' else 'identical or one mismatch'}); see Host_Self / Host_Excluded.")
    return pep
