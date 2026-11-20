"""A prediction cache: a score computed once for a (tool, settings, allele, peptide) is never computed again.

NetMHCpan took 40 minutes for one 19,000-peptide x 37-allele run; run it again with one more allele, or the same protein under another
strategy, and only what is new is computed. Entries are keyed by everything that changes the numbers (tool path, its options), stored per
allele in the workspace (workspace/cache/<tool>-<hash>/<allele>.pkl.gz), and checked by peptide.
"""
import hashlib
import os
import shutil
import threading

import pandas as pd

from . import config

_LOCK = threading.RLock()


def enabled():
    return bool(config.load().get("performance", {}).get("cache", 1))


def root():
    d = os.path.join(config.workspace(), "cache")
    os.makedirs(d, exist_ok=True)
    return d


def _dir(tool, sig):
    h = hashlib.sha1(repr(sig).encode()).hexdigest()[:10]
    d = os.path.join(root(), f"{tool}-{h}")
    os.makedirs(d, exist_ok=True)
    return d


def _path(d, allele):
    return os.path.join(d, allele.replace("*", "").replace(":", "").replace("/", "_") + ".pkl.gz")


def _load(d, allele):
    try:
        return pd.read_pickle(_path(d, allele))
    except (OSError, ValueError, EOFError):
        return None


def _save(d, allele, df):
    path = _path(d, allele)
    tmp = path + ".tmp"
    df.to_pickle(tmp, compression="gzip")
    os.replace(tmp, path)


def fetch(tool, sig, peptides, alleles, compute, canon, job=None, stage=None):
    """Table tools (several columns per row, one Allele column): compute(missing_peptides, alleles) -> DataFrame[Peptide, Allele, ...].
    `alleles` are the names the tool wants; `canon` maps a name to the key used for storing and matching. Returns the same DataFrame
    compute would, for every peptide and allele asked for, computing only what is not cached."""
    peps = sorted(set(peptides))
    if not enabled() or not peps or not alleles:
        return compute(peps, list(alleles))
    with _LOCK:
        d = _dir(tool, sig)
        store, missing = {}, {}
        for a in alleles:
            df = _load(d, canon(a))
            store[a] = df
            have = set() if df is None else set(df["Peptide"])
            miss = [p for p in peps if p not in have]
            if miss:
                missing[a] = miss
        total, todo = len(peps) * len(alleles), sum(len(m) for m in missing.values())
        if job:
            job.log(f"cache: {total - todo:,} of {total:,} {tool} predictions already known" + (f", computing {todo:,}" if todo else ""))
        groups = {}
        for a, m in missing.items():
            groups.setdefault(hash(tuple(m)), ([], m))[0].append(a)
        for als, m in groups.values():
            res = compute(m, als)
            if res is None or res.empty:
                continue
            ck = res["Allele"].map(canon)
            for a in als:
                new = res[ck == canon(a)]
                if new.empty:
                    continue
                old = store[a]
                merged = new if old is None else pd.concat([old, new], ignore_index=True)
                store[a] = merged.drop_duplicates("Peptide", keep="last").reset_index(drop=True)
                _save(d, canon(a), store[a])
        parts = [store[a][store[a]["Peptide"].isin(peps)] for a in alleles if store[a] is not None]
        parts = [p for p in parts if not p.empty]
        if job and stage and not todo:
            job.update(stage, done=1, total=1, msg="all predictions were already known (cache)")
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def fetch_scores(tool, sig, pairs, compute, canon, job=None):
    """Pair tools (one score per peptide-allele pair): pairs = DataFrame[peptide, allele]; compute(missing_pairs) -> array of scores in
    that order. Returns the scores aligned to `pairs`."""
    if not enabled() or pairs.empty:
        return compute(pairs)
    with _LOCK:
        d = _dir(tool, sig)
        out = pd.Series(float("nan"), index=pairs.index)
        miss_idx, store = [], {}
        for a, g in pairs.groupby("allele"):
            df = _load(d, canon(a))
            store[a] = df if df is not None else pd.DataFrame({"Peptide": [], "score": []})
            known = store[a].drop_duplicates("Peptide").set_index("Peptide")["score"]
            hit = g["peptide"].map(known)
            out.loc[g.index] = hit
            miss_idx += list(g.index[hit.isna()])
        if job:
            job.log(f"cache: {len(pairs) - len(miss_idx):,} of {len(pairs):,} {tool} pairs already known" + (f", computing {len(miss_idx):,}" if miss_idx else ""))
        if miss_idx:
            sub = pairs.loc[miss_idx]
            scores = pd.Series(compute(sub), index=sub.index)
            out.loc[miss_idx] = scores
            for a, g in sub.groupby("allele"):
                add = pd.DataFrame({"Peptide": g["peptide"].values, "score": scores.loc[g.index].values})
                merged = pd.concat([store[a], add], ignore_index=True).drop_duplicates("Peptide", keep="last")
                _save(d, canon(a), merged.reset_index(drop=True))
        return out.values


def size():
    r = root()
    return sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(r) for f in fs)


def clear():
    r = root()
    for name in os.listdir(r):
        shutil.rmtree(os.path.join(r, name), ignore_errors=True)
