#!/usr/bin/env python3
"""Home Designer — tiny local server (Python standard library only).

Serves the web app from ./web and stores designs as JSON files in ./designs
so they can be committed to git. Every save keeps the previous version in
./designs/.history/ (git-ignored).

    python3 server.py                 # http://127.0.0.1:8765
    python3 server.py --port 9000 --open
    python3 server.py --host 0.0.0.0  # reachable from other devices on your network
"""
from __future__ import annotations

import argparse
import base64
import json
import mimetypes
import os
import re
import shutil
import sys
import time
import webbrowser
from datetime import datetime, timezone
from http import HTTPStatus
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, unquote, urlparse
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
WEB = ROOT / "web"
DESIGNS = Path(os.environ.get("HOME_DESIGNER_DATA", ROOT / "designs")).resolve()
HISTORY = DESIGNS / ".history"
MODELS = DESIGNS / "models"          # downloaded / uploaded 3D models (glTF / GLB)
PH_API = "https://api.polyhaven.com"
_ph_cache: dict = {}
MAX_BYTES = 25 * 1024 * 1024
MAX_IMAGE = 15 * 1024 * 1024
KEEP_VERSIONS = 20
NAME_RE = re.compile(r"^[A-Za-z0-9 _.\-]{1,80}$")

mimetypes.add_type("text/javascript", ".js")
mimetypes.add_type("text/css", ".css")


def safe_name(name: str) -> str | None:
    name = name.strip()
    if not NAME_RE.match(name) or name.startswith(".") or ".." in name:
        return None
    return name


def design_path(name: str) -> Path:
    return DESIGNS / f"{name}.json"


def backup(path: Path, name: str) -> None:
    if not path.exists():
        return
    HISTORY.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    shutil.copy2(path, HISTORY / f"{name}__{stamp}.json")
    old = sorted(HISTORY.glob(f"{glob_escape(name)}__*.json"))
    for p in old[:-KEEP_VERSIONS]:
        p.unlink(missing_ok=True)


def http_get(url: str, limit: int) -> bytes:
    req = Request(url, headers={"User-Agent": "HomeDesigner (local design tool)"})
    with urlopen(req, timeout=60) as r:
        data = r.read(limit + 1)
    if len(data) > limit:
        raise ValueError("file too large")
    return data


def glob_escape(s: str) -> str:
    return re.sub(r"([\[\]\*\?])", r"[\1]", s)


class Handler(SimpleHTTPRequestHandler):
    server_version = "HomeDesigner/1.0"

    def __init__(self, *a, **kw):
        super().__init__(*a, directory=str(WEB), **kw)

    def translate_path(self, path):
        # /models/... is served from designs/models so models live (and get committed) with designs
        p = urlparse(path).path
        if p.startswith("/models/"):
            rel = Path(unquote(p[len("/models/"):]))
            target = (MODELS / rel).resolve()
            if MODELS.resolve() in target.parents:
                return str(target)
            return str(MODELS / "__invalid__")
        return super().translate_path(path)

    def log_message(self, fmt, *args):  # quieter log
        if "/api/" in (self.path or ""):
            sys.stderr.write("%s %s\n" % (self.command, self.path))

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    # ---- helpers ----
    def send_json(self, obj, status=HTTPStatus.OK):
        body = json.dumps(obj).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def send_err(self, status, msg):
        body = msg.encode()
        self.send_response(status)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def api_name(self) -> str | None:
        path = urlparse(self.path).path
        prefix = "/api/designs/"
        if not path.startswith(prefix):
            return None
        return safe_name(unquote(path[len(prefix):]))

    # ---- routes ----
    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/ping":
            return self.send_json({"ok": True, "designs": str(DESIGNS)})
        if path == "/api/designs":
            items = []
            for p in sorted(DESIGNS.glob("*.json")):
                st = p.stat()
                title = p.stem
                try:
                    with p.open(encoding="utf-8") as f:
                        title = json.load(f).get("name") or p.stem
                except Exception:
                    pass
                items.append({
                    "name": p.stem,
                    "title": title,
                    "size": st.st_size,
                    "updated": datetime.fromtimestamp(st.st_mtime, timezone.utc).isoformat(),
                })
            items.sort(key=lambda d: d["updated"], reverse=True)
            return self.send_json({"designs": items})
        if path == "/api/fetch-image":
            return self.fetch_image()
        if path == "/api/polyhaven/models":
            return self.ph_models()
        if path == "/api/polyhaven/download":
            return self.ph_download()
        if path == "/api/models":
            items = [str(p.relative_to(MODELS)).replace(os.sep, "/") for p in MODELS.rglob("*") if p.suffix.lower() in (".glb", ".gltf")] if MODELS.exists() else []
            return self.send_json({"models": sorted(items)})
        if path.startswith("/api/designs/"):
            name = self.api_name()
            if not name:
                return self.send_err(HTTPStatus.BAD_REQUEST, "Bad design name")
            p = design_path(name)
            if not p.exists():
                return self.send_err(HTTPStatus.NOT_FOUND, "No such design")
            data = p.read_bytes()
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        if path.startswith("/api/"):
            return self.send_err(HTTPStatus.NOT_FOUND, "Unknown API route")
        return super().do_GET()

    def fetch_image(self):
        """Download an image from the web for tracing (avoids browser CORS limits)."""
        url = (parse_qs(urlparse(self.path).query).get("url") or [""])[0].strip()
        if not re.match(r"^https?://", url, re.I):
            return self.send_err(HTTPStatus.BAD_REQUEST, "Only http(s) image URLs are supported")
        try:
            req = Request(url, headers={"User-Agent": "Mozilla/5.0 (HomeDesigner image import)"})
            with urlopen(req, timeout=20) as r:
                ctype = (r.headers.get("Content-Type") or "").split(";")[0].strip().lower()
                data = r.read(MAX_IMAGE + 1)
        except Exception as e:  # network errors, 404s, timeouts
            return self.send_err(HTTPStatus.BAD_GATEWAY, f"Download failed: {e}")
        if len(data) > MAX_IMAGE:
            return self.send_err(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, "Image is larger than 15 MB")
        if not ctype.startswith("image/"):
            guess = mimetypes.guess_type(urlparse(url).path)[0] or ""
            if not guess.startswith("image/"):
                return self.send_err(HTTPStatus.UNSUPPORTED_MEDIA_TYPE, f"That address is not an image ({ctype or 'unknown type'}). Open the image itself and copy its address.")
            ctype = guess
        return self.send_json({"dataUrl": f"data:{ctype};base64,{base64.b64encode(data).decode()}"})

    # ---- 3D models ----
    def ph_models(self):
        """List Poly Haven's free (CC0) 3D models, cached for the session."""
        if "models" not in _ph_cache:
            try:
                _ph_cache["models"] = json.loads(http_get(f"{PH_API}/assets?t=models", 8 * 1024 * 1024))
            except Exception as e:
                return self.send_err(HTTPStatus.BAD_GATEWAY, f"Could not reach Poly Haven: {e}")
        out = []
        for mid, info in _ph_cache["models"].items():
            out.append({"id": mid, "name": info.get("name", mid), "categories": info.get("categories", []),
                        "tags": info.get("tags", [])[:12], "downloads": info.get("download_count", 0)})
        out.sort(key=lambda d: -d["downloads"])
        return self.send_json({"models": out})

    def ph_download(self):
        """Download one Poly Haven model (glTF, 1k textures) into designs/models/polyhaven/<id>/."""
        mid = (parse_qs(urlparse(self.path).query).get("id") or [""])[0]
        res = (parse_qs(urlparse(self.path).query).get("res") or ["1k"])[0]
        if not re.match(r"^[A-Za-z0-9_\-]{1,80}$", mid) or res not in ("1k", "2k"):
            return self.send_err(HTTPStatus.BAD_REQUEST, "Bad model id")
        folder = MODELS / "polyhaven" / mid
        main = folder / f"{mid}_{res}.gltf"
        if not main.exists():
            try:
                files = json.loads(http_get(f"{PH_API}/files/{mid}", 4 * 1024 * 1024))
                g = files["gltf"][res]["gltf"]
                folder.mkdir(parents=True, exist_ok=True)
                for rel, inc in (g.get("include") or {}).items():
                    target = (folder / rel).resolve()
                    if folder.resolve() not in target.parents:
                        continue
                    target.parent.mkdir(parents=True, exist_ok=True)
                    target.write_bytes(http_get(inc["url"], 80 * 1024 * 1024))
                main.write_bytes(http_get(g["url"], 20 * 1024 * 1024))
            except Exception as e:
                shutil.rmtree(folder, ignore_errors=True)
                return self.send_err(HTTPStatus.BAD_GATEWAY, f"Download failed: {e}")
        return self.send_json({"ok": True, "url": f"models/polyhaven/{mid}/{main.name}"})

    def upload_model(self):
        name = unquote(urlparse(self.path).path[len("/api/models/upload/"):])
        if not re.match(r"^[A-Za-z0-9 _.\-]{1,80}\.glb$", name, re.I) or name.startswith("."):
            return self.send_err(HTTPStatus.BAD_REQUEST, "Only .glb files with a simple name can be uploaded")
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0 or length > 60 * 1024 * 1024:
            return self.send_err(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, "Model is empty or larger than 60 MB")
        data = self.rfile.read(length)
        if data[:4] != b"glTF":
            return self.send_err(HTTPStatus.BAD_REQUEST, "Not a binary glTF (.glb) file")
        (MODELS / "uploads").mkdir(parents=True, exist_ok=True)
        (MODELS / "uploads" / name).write_bytes(data)
        return self.send_json({"ok": True, "url": f"models/uploads/{name}"})

    def do_PUT(self):
        if urlparse(self.path).path.startswith("/api/models/upload/"):
            return self.upload_model()
        name = self.api_name()
        if not name:
            return self.send_err(HTTPStatus.BAD_REQUEST, "Bad design name (letters, digits, space, _ . - only)")
        length = int(self.headers.get("Content-Length") or 0)
        if length <= 0 or length > MAX_BYTES:
            return self.send_err(HTTPStatus.REQUEST_ENTITY_TOO_LARGE, "Design is empty or too large")
        raw = self.rfile.read(length)
        try:
            data = json.loads(raw)
        except json.JSONDecodeError as e:
            return self.send_err(HTTPStatus.BAD_REQUEST, f"Invalid JSON: {e}")
        if not isinstance(data, dict) or data.get("format") != "home-designer" or not isinstance(data.get("floors"), list):
            return self.send_err(HTTPStatus.BAD_REQUEST, "Not a Home Designer design")
        DESIGNS.mkdir(parents=True, exist_ok=True)
        p = design_path(name)
        backup(p, name)
        tmp = p.with_suffix(".json.tmp")
        with tmp.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=1, ensure_ascii=False)
            f.write("\n")
        os.replace(tmp, p)
        return self.send_json({"ok": True, "name": name, "path": str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)})

    do_POST = do_PUT

    def do_DELETE(self):
        name = self.api_name()
        if not name:
            return self.send_err(HTTPStatus.BAD_REQUEST, "Bad design name")
        p = design_path(name)
        if not p.exists():
            return self.send_err(HTTPStatus.NOT_FOUND, "No such design")
        backup(p, name)
        p.unlink()
        return self.send_json({"ok": True})


def open_browser(url: str) -> None:
    """Open the app in a browser. Under WSL use the Windows browser; never fail or spam errors."""
    import subprocess
    try:
        if "microsoft" in Path("/proc/version").read_text().lower():
            for cmd in (["wslview", url], ["cmd.exe", "/c", "start", "", url], ["explorer.exe", url]):
                try:
                    subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
                    return
                except FileNotFoundError:
                    continue
    except OSError:
        pass
    # quiet the "x-www-browser: not found" noise that xdg-open prints on systems without a browser
    devnull = os.open(os.devnull, os.O_WRONLY)
    saved = os.dup(2)
    try:
        os.dup2(devnull, 2)
        ok = webbrowser.open(url)
    except Exception:
        ok = False
    finally:
        os.dup2(saved, 2)
        os.close(devnull)
        os.close(saved)
    if not ok:
        print(f"(Could not open a browser automatically - open {url} yourself.)")


def main():
    ap = argparse.ArgumentParser(description="Home Designer local server")
    ap.add_argument("--host", default="127.0.0.1", help="interface to bind (default 127.0.0.1; use 0.0.0.0 for LAN)")
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--open", action="store_true", help="open the browser")
    args = ap.parse_args()

    DESIGNS.mkdir(parents=True, exist_ok=True)
    httpd = ThreadingHTTPServer((args.host, args.port), Handler)
    url = f"http://{'127.0.0.1' if args.host in ('0.0.0.0', '') else args.host}:{args.port}/"
    print(f"Home Designer running at {url}   <- open this in your browser")
    print(f"Designs are saved in {DESIGNS}")
    print("Press Ctrl+C to stop.")
    if args.open:
        open_browser(url)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nbye")


if __name__ == "__main__":
    main()
