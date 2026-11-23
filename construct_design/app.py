#!/usr/bin/env python3
"""
Gene Construct Builder - Flask LAN server.

Everything lives in one directory: this file, the HTML page, and the
constructs themselves as plain .json files you can read, copy or hand-edit.

    pip install flask
    python3 app.py                      # http://<your-ip>:8080/
    python3 app.py --port 5000 --data ./constructs

There is NO login: anyone who can reach the port can read, change and delete
the constructs in the data directory. Run it only on a trusted LAN.
"""
import argparse
import json
import os
import re
import socket
from datetime import datetime, timezone

from flask import Flask, jsonify, request, send_from_directory

HERE = os.path.dirname(os.path.abspath(__file__))
UI_FILE = "gene-construct-builder-server.html"

app = Flask(__name__)
DATA_DIR = HERE  # replaced in main(); the one directory holding the .json files


# ---------- helpers ----------
def now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def slug(name):
    """A safe, readable filename stem built from the construct name."""
    s = re.sub(r"[^a-z0-9]+", "-", str(name or "").strip().lower()).strip("-")[:60]
    return s or "construct"


def path_for(cid):
    """Resolve an id to a path inside DATA_DIR, or None if it tries to escape."""
    if not isinstance(cid, str) or not re.fullmatch(r"[A-Za-z0-9._-]{1,80}", cid) or cid.startswith("."):
        return None
    p = os.path.join(DATA_DIR, cid + ".json")
    return p if os.path.dirname(os.path.abspath(p)) == os.path.abspath(DATA_DIR) else None


def unique_id(stem, except_id=None):
    cid, n = stem, 2
    while cid != except_id and os.path.exists(os.path.join(DATA_DIR, cid + ".json")):
        cid, n = f"{stem}-{n}", n + 1
    return cid


def read_construct(cid):
    p = path_for(cid)
    if not p or not os.path.exists(p):
        return None
    try:
        with open(p, encoding="utf-8") as fh:
            doc = json.load(fh)
    except (OSError, ValueError):
        return None
    if not isinstance(doc.get("blocks"), list):
        return None
    mtime = datetime.fromtimestamp(os.path.getmtime(p), timezone.utc).isoformat(timespec="seconds")
    return {
        "id": cid,
        "name": doc.get("name") or cid,
        "settings": doc.get("settings") or {},
        "blocks": doc["blocks"],
        "created": doc.get("created") or mtime,
        "updated": doc.get("updated") or mtime,
    }


def clean(body, prev=None):
    """Keep only the fields a construct may have; raise ValueError if unusable."""
    if not isinstance(body, dict):
        raise ValueError("body must be a JSON object")
    blocks = body.get("blocks")
    if not isinstance(blocks, list):
        raise ValueError('"blocks" must be an array')
    if len(blocks) > 2000:
        raise ValueError("too many elements")
    settings = body.get("settings")
    name = body.get("name") or (settings or {}).get("title") or "Untitled construct"
    return {
        "name": str(name)[:200],
        "settings": settings if isinstance(settings, dict) else {},
        "blocks": blocks,
        "created": (prev or {}).get("created") or body.get("created") or now(),
        "updated": now(),
    }


def write_construct(cid, doc):
    p = path_for(cid)
    if not p:
        raise ValueError("bad id")
    tmp = p + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(doc, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    os.replace(tmp, p)  # readers never see a half-written file


# ---------- pages ----------
@app.get("/")
def index():
    return send_from_directory(HERE, UI_FILE)


@app.get("/files/<path:name>")
def raw_file(name):
    """Download one construct file as-is."""
    cid = name[:-5] if name.endswith(".json") else name
    if not path_for(cid):
        return jsonify(error="bad name"), 400
    return send_from_directory(DATA_DIR, cid + ".json", mimetype="application/json")


# ---------- API (full CRUD) ----------
@app.get("/api/constructs")
def list_constructs():
    items = []
    for fn in sorted(os.listdir(DATA_DIR)):
        if not fn.endswith(".json"):
            continue
        cid = fn[:-5]
        c = read_construct(cid)
        if c:
            items.append({"id": cid, "name": c["name"], "created": c["created"],
                          "updated": c["updated"], "elements": len(c["blocks"]), "file": fn})
        else:
            items.append({"id": cid, "name": fn, "invalid": True, "elements": 0,
                          "updated": "", "file": fn})
    items.sort(key=lambda c: c["updated"], reverse=True)
    return jsonify(dataDir=os.path.abspath(DATA_DIR), constructs=items)


@app.get("/api/constructs/<cid>")
def get_construct(cid):
    c = read_construct(cid)
    return (jsonify(c), 200) if c else (jsonify(error="no such construct"), 404)


@app.post("/api/constructs")
def create_construct():
    doc = clean(request.get_json(force=True, silent=False))
    cid = unique_id(slug(doc["name"]))
    write_construct(cid, doc)
    return jsonify(id=cid, **doc), 201


@app.route("/api/constructs/<cid>", methods=["PUT", "PATCH"])
def update_construct(cid):
    prev = read_construct(cid)
    if not prev:
        return jsonify(error="no such construct"), 404
    doc = clean(request.get_json(force=True, silent=False), prev)
    new_id = cid
    if slug(doc["name"]) != slug(prev["name"]):          # name changed -> rename the file too
        new_id = unique_id(slug(doc["name"]), cid)
        if new_id != cid:
            write_construct(new_id, doc)
            try:
                os.remove(path_for(cid))
            except OSError:
                pass
    if new_id == cid:
        write_construct(cid, doc)
    return jsonify(id=new_id, **doc)


@app.delete("/api/constructs/<cid>")
def delete_construct(cid):
    p = path_for(cid)
    if not p:
        return jsonify(error="bad id"), 400
    try:
        os.remove(p)
    except OSError:
        return jsonify(error="no such construct"), 404
    return jsonify(deleted=cid)


@app.after_request
def no_cache(resp):
    resp.headers["Cache-Control"] = "no-store"
    resp.headers["Access-Control-Allow-Origin"] = "*"
    return resp


@app.errorhandler(ValueError)
def bad_value(e):
    return jsonify(error=str(e)), 400


@app.errorhandler(400)
def bad_request(e):
    return jsonify(error="invalid JSON in request body"), 400


# ---------- start ----------
def main():
    global DATA_DIR
    ap = argparse.ArgumentParser(description="Gene Construct Builder LAN server")
    ap.add_argument("--data", default=HERE, help="directory holding the construct .json files (default: this folder)")
    ap.add_argument("--port", type=int, default=8080)
    ap.add_argument("--host", default="0.0.0.0", help="0.0.0.0 serves the whole LAN")
    args = ap.parse_args()

    DATA_DIR = os.path.abspath(args.data)
    os.makedirs(DATA_DIR, exist_ok=True)

    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
    except OSError:
        ip = "your-ip"
    print("Gene Construct Builder")
    print("  data directory :", DATA_DIR)
    print("  local          : http://localhost:%d/" % args.port)
    print("  on your LAN    : http://%s:%d/" % (ip, args.port))
    print("  WARNING: no login - anyone who can reach this port can edit or delete these constructs.")
    app.run(host=args.host, port=args.port, threaded=True)


if __name__ == "__main__":
    main()
