"""The workspace: every run is kept on disk (input, settings, dates, state, tables, reports) and can be listed, reopened,
renamed, annotated or deleted later, also after the server was restarted.

Layout (config.workspace()):   runs/<id>/run.json        what was asked, when, how it ended, where the files are
                               runs/<id>/<id>_*.csv ...  the outputs of that run
run.json holds the full settings (the same json the page saves), so a past run can be re-run or used as a template.
"""
import json
import os
import re
import shutil
import time

from . import config
from .jobs import Job

ENGINE_VERSION = "2026.10"
RUN = "{RUN}"


def _file(rid):
    return os.path.join(config.run_dir(rid), "run.json")


def _swap(obj, a, b):
    """Replace the leading path `a` by `b` in every string of a json-like object (keeps the workspace relocatable)."""
    if isinstance(obj, str):
        return b + obj[len(a):] if obj.startswith(a) else obj
    if isinstance(obj, list):
        return [_swap(x, a, b) for x in obj]
    if isinstance(obj, dict):
        return {k: _swap(v, a, b) for k, v in obj.items()}
    return obj


def _write(rid, meta):
    path = _file(rid)
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        json.dump(meta, f, indent=1, default=str)
    os.replace(tmp, path)


def read(rid):
    try:
        with open(_file(rid)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def describe_input(params):
    text = params.get("text") or ""
    names = re.findall(r"^>\s*(\S+)", text, flags=re.M)
    residues = sum(len(re.sub(r"[^A-Za-z]", "", ln)) for ln in text.splitlines() if not ln.startswith(">"))
    return {"proteins": names[:20], "n_proteins": len(names), "residues": residues, "mode": params.get("mode", "scan"),
            "sha1": __import__("hashlib").sha1(text.encode()).hexdigest()[:10]}


def default_name(kind, params):
    if kind == "single":
        return f"{params.get('model', 'model')} · {(describe_input(params)['proteins'] or ['peptides'])[0]}"
    return params.get("title") or "Epitope prediction report"


def create(job, kind, params, name=None):
    job.id = job.id            # (kept: the id is the folder name)
    _write(job.id, {"id": job.id, "kind": kind, "name": name or default_name(kind, params), "created": time.time(), "state": "queued",
                    "params": params, "input": describe_input(params), "engine_version": ENGINE_VERSION, "notes": "", "starred": False})


def finish(job):
    meta = read(job.id)
    if meta is None:
        return
    d = config.run_dir(job.id)
    snap = job.snapshot(with_log=False)
    meta.update({"state": job.state, "started": job.started, "finished": job.finished, "error": job.error,
                 "summary": (job.result or {}).get("summary"), "stages": snap["stages"], "log": list(job.logs[-120:]),
                 "result": _swap(job.result, d, RUN) if job.result else None, "tasks": job.tasks_snapshot() if job._tasks else None,
                 "tools": {k: {kk: vv for kk, vv in v.items() if kk in ("home", "python", "command", "device")}
                           for k, v in config.effective()["models"].items() if v.get("home") or v.get("python")}})
    _write(job.id, meta)


def listing():
    out = []
    root = os.path.join(config.workspace(), "runs")
    for rid in os.listdir(root):
        m = read(rid)
        if not m:
            continue
        size = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(os.path.join(root, rid)) for f in fs)
        out.append({k: m.get(k) for k in ("id", "kind", "name", "created", "started", "finished", "state", "summary", "input", "notes", "starred", "error", "engine_version")}
                   | {"size": size, "model": (m.get("params") or {}).get("model"),
                      "classes": [c for c in ("c1", "c2") if (m.get("params") or {}).get(c, {}).get("enabled", True) and (m.get("params") or {}).get(c)]})
    return sorted(out, key=lambda r: -(r["created"] or 0))


def update(rid, **fields):
    m = read(rid)
    if m is None:
        return None
    for k in ("name", "notes", "starred"):
        if k in fields:
            m[k] = fields[k]
    _write(rid, m)
    return m


def delete(rid):
    if not re.fullmatch(r"[0-9a-f]{12}", str(rid)):
        raise ValueError("bad run id")
    d = os.path.join(config.workspace(), "runs", rid)
    if os.path.isdir(d):
        shutil.rmtree(d)


def load_job(rid):
    """A finished Job rebuilt from disk (None if there is no such run), so every results endpoint works on old runs."""
    m = read(rid) if re.fullmatch(r"[0-9a-f]{12}", str(rid)) else None
    if not m:
        return None
    j = Job(m["kind"])
    j.id, j.state, j.started, j.finished, j.error = rid, m["state"], m.get("started"), m.get("finished"), m.get("error")
    j.stages, j.logs = m.get("stages") or [], m.get("log") or []
    j.created = m.get("created") or time.time()
    j._saved_tasks = m.get("tasks")
    if m.get("result"):
        j.result = _swap(m["result"], RUN, config.run_dir(rid))
    return j


def mark_interrupted():
    """Runs that were queued or running when the server stopped can never finish: say so."""
    for r in listing():
        if r["state"] in ("queued", "running"):
            m = read(r["id"])
            m.update(state="interrupted", error="The server stopped before this run finished.", finished=time.time())
            _write(r["id"], m)
