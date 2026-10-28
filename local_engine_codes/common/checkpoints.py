"""Finished classes are remembered, so a run that stops later (a tool fails, the browser closes, you press Cancel) does not lose them.

When MHC class I finishes, its result is stored under a fingerprint of everything that decides it (the protein, sizes, alleles, models,
funnel, filters, strategy, engine version, population file). Any later run whose class I has the same fingerprint takes the stored result
instead of running again: Resume is simply "run again". The files are hard-linked into the new run's folder, so deleting the old run
does not break the new one.
"""
import hashlib
import json
import os
import shutil
import time

from . import config

RUN = "{RUN}"
PATH_KEYS = (".csv", ".json", ".html", ".zip", ".gz")


def _dir():
    d = os.path.join(config.workspace(), "checkpoints")
    os.makedirs(d, exist_ok=True)
    return d


def fingerprint(norm, cls, version, population_file=None):
    """norm: the class's normalised parameters (alleles are order-free)."""
    body = {k: v for k, v in norm.items() if k not in ("title",)}
    body["alleles"] = sorted(body.get("alleles", []))
    pop = None
    if population_file and os.path.exists(population_file):
        st = os.stat(population_file)
        pop = [population_file, int(st.st_mtime), st.st_size]
    return hashlib.sha1(json.dumps({"c": cls, "p": body, "v": version, "pop": pop}, sort_keys=True, default=str).encode()).hexdigest()[:20]


def _paths(obj):
    """Every file path in a result (any nesting)."""
    if isinstance(obj, str):
        if os.sep in obj and obj.endswith(PATH_KEYS):
            yield obj
    elif isinstance(obj, dict):
        for k, v in obj.items():
            yield from _paths(v)
    elif isinstance(obj, list):
        for v in obj:
            yield from _paths(v)


def _swap(obj, a, b):
    if isinstance(obj, str):
        return b + obj[len(a):] if obj.startswith(a) else obj
    if isinstance(obj, list):
        return [_swap(x, a, b) for x in obj]
    if isinstance(obj, dict):
        return {k: _swap(v, a, b) for k, v in obj.items()}
    return obj


def save(fp, run_id, cls, result):
    d = config.run_dir(run_id)
    meta = {"fp": fp, "run_id": str(run_id).split("_")[0], "cls": cls, "saved": time.time(), "result": _swap(result, d, RUN)}
    path = os.path.join(_dir(), fp + ".json")
    with open(path + ".tmp", "w") as f:
        json.dump(meta, f, default=str)
    os.replace(path + ".tmp", path)


def load(fp, dest_run_id):
    """-> (result, source run id, age in seconds) with the files linked into dest_run_id's folder, or None."""
    try:
        with open(os.path.join(_dir(), fp + ".json")) as f:
            meta = json.load(f)
    except (OSError, ValueError):
        return None
    src = config.run_dir(meta["run_id"])
    result = _swap(meta["result"], RUN, src)
    if not all(os.path.exists(p) for p in _paths(result)):
        return None
    dest = config.run_dir(dest_run_id)
    if os.path.abspath(dest) != os.path.abspath(src):
        for p in set(_paths(result)):
            q = os.path.join(dest, os.path.basename(p))
            if not os.path.exists(q):
                try:
                    os.link(p, q)
                except OSError:
                    shutil.copy2(p, q)
        result = _swap(result, src, dest)
    return result, meta["run_id"], time.time() - meta["saved"]
