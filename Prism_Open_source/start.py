#!/usr/bin/env python3
"""OpenPrism server launcher — needs only Python 3.10+.

    python3 start.py
        Start the server (the first run installs the packages into .venv/ next to this file).

    python3 start.py --reset --reset --reset
        Delete all current data (accounts, projects, shares, settings) and start fresh.

    python3 start.py --reset --reset --reset --loadolddata BACKUP
        Delete all current data, load an old backup instead, and start. Use it for the first start after
        installing a new version. BACKUP is a backup .zip file or the folder that contains it (the newest
        backup in the folder is used). The backup is checked before anything is deleted.

    python3 start.py --selfpurge --selfpurge --selfpurge --output BACKUP_FOLDER
        Save all data as a backup .zip in BACKUP_FOLDER, check it, then delete this OpenPrism installation
        completely (code, Python environment, data). Afterwards clone a new copy and start it with --loadolddata.

Repeating a flag three times is the confirmation. Stop a running server (Ctrl+C) before these commands.
Other options are passed to the server, e.g.  python3 start.py --port 8443  (see python3 run.py --help).
"""
from __future__ import annotations

import datetime
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
VENV = HERE / ".venv"
REQ = HERE / "requirements.txt"
STAMP = VENV / ".requirements.sha256"
DATA = Path(os.environ.get("OPENPRISM_DATA") or HERE / "projects").resolve()
DB = DATA / "openprism.db"
PIDFILE = DATA / "server.pid"
WIN = os.name == "nt"
EXPORT_FORMAT = "openprism-server-export-1"


def say(msg=""):
    print(msg, flush=True)


def fail(msg):
    sys.exit(f"\n  {msg}\n")


# ---------------------------------------------------------------- environment
def venv_python() -> Path:
    return VENV / ("Scripts/python.exe" if WIN else "bin/python")


def create_venv():
    say("  Creating a private Python environment in .venv/ …")
    try:
        subprocess.run([sys.executable, "-m", "venv", str(VENV)], check=True)
    except subprocess.CalledProcessError:
        shutil.rmtree(VENV, ignore_errors=True)
        fail("Could not create the environment. On Debian/Ubuntu install the venv module first:\n"
             "      sudo apt install python3-venv\n  then run  python3 start.py  again.")


def install_requirements():
    want = hashlib.sha256(REQ.read_bytes()).hexdigest()
    if STAMP.exists() and STAMP.read_text().strip() == want:
        return
    say("  Installing the required packages (first run only, about a minute) …")
    py = str(venv_python())
    subprocess.run([py, "-m", "pip", "install", "--quiet", "--upgrade", "pip"], check=False)
    r = subprocess.run([py, "-m", "pip", "install", "--quiet", "-r", str(REQ)])
    if r.returncode != 0:
        fail("Installing the packages failed (see the messages above). Check the internet connection and run\n"
             "      python3 start.py\n  again.")
    STAMP.write_text(want)
    say("  Packages installed.")


# ---------------------------------------------------------------- safety checks
def running_server():
    """(pid, port) of a server using this data folder, or None. run.py writes the marker file."""
    try:
        info = json.loads(PIDFILE.read_text())
        pid = int(info["pid"])
    except (OSError, ValueError, KeyError, TypeError):
        return None
    try:
        if WIN:
            out = subprocess.run(["tasklist", "/FI", f"PID eq {pid}"], capture_output=True, text=True).stdout
            alive = str(pid) in out
        else:
            os.kill(pid, 0)
            alive = True
    except (OSError, ProcessLookupError):
        alive = False
    if not alive:
        PIDFILE.unlink(missing_ok=True)
        return None
    return pid, info.get("port")


def require_stopped(action):
    r = running_server()
    if r:
        fail(f"The OpenPrism server is running (process {r[0]}, port {r[1]}). Stop it first (Ctrl+C in its window), then {action}.")


def data_files():
    return [p for p in DATA.iterdir() if p.name != ".gitkeep"] if DATA.exists() else []


def this_schema() -> int:
    """Database schema version this copy of the code understands (read from the source, no imports needed)."""
    import re
    try:
        m = re.search(r"^SCHEMA_VERSION\s*=\s*(\d+)", (HERE / "openprism" / "store.py").read_text(), re.M)
        return int(m.group(1)) if m else 0
    except OSError:
        return 0


def db_schema(db_path: Path) -> int | None:
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        row = con.execute("SELECT value FROM settings WHERE key = 'schema'").fetchone()
    except sqlite3.DatabaseError:
        return None
    finally:
        con.close()
    try:
        return int(row[0]) if row else None
    except (TypeError, ValueError):
        return None


def counts(db_path: Path):
    """(users, projects) in a database file, for messages and checks."""
    con = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    try:
        users = con.execute("SELECT COUNT(*) FROM users").fetchone()[0]
        projects = con.execute("SELECT COUNT(*) FROM projects").fetchone()[0]
        ok = con.execute("PRAGMA integrity_check").fetchone()[0]
    finally:
        con.close()
    if ok != "ok":
        raise ValueError(f"database integrity check failed: {ok}")
    return users, projects


def copy_db(src: Path, dest: Path):
    """Consistent copy of a (possibly WAL-mode) SQLite database."""
    s = sqlite3.connect(str(src))
    try:
        d = sqlite3.connect(str(dest))
        with d:
            s.backup(d)
        d.close()
    finally:
        s.close()


# ---------------------------------------------------------------- export / import
def export(output: str | None) -> Path:
    if not output:
        fail("Say where to save the backup.")
    out = Path(output).expanduser().resolve()
    if out.is_dir():
        out = out / f"openprism-backup-{datetime.datetime.now():%Y%m%d-%H%M%S}.zip"
    if out.suffix.lower() != ".zip":
        out = out.with_suffix(out.suffix + ".zip") if out.suffix else out.with_suffix(".zip")
    if not DB.exists():
        fail(f"There is no server data to export in {DATA}.")
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        snap = Path(tmp) / "openprism.db"
        copy_db(DB, snap)
        users, projects = counts(snap)
        manifest = {"format": EXPORT_FORMAT, "created": datetime.datetime.now().isoformat(timespec="seconds"),
                    "openprism": _version(), "schema": db_schema(snap), "users": users, "projects": projects}
        part = out.with_name(out.name + ".part")
        with zipfile.ZipFile(part, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("manifest.json", json.dumps(manifest, indent=2))
            z.write(snap, "openprism.db")
            tls = DATA / "tls"
            if tls.is_dir():
                for f in tls.iterdir():
                    if f.is_file():
                        z.write(f, f"tls/{f.name}")
        part.replace(out)
    try:
        os.chmod(out, 0o600)        # contains password hashes and the TLS key
    except OSError:
        pass
    verify_export(out)
    say(f"  Backed up {users} account(s) and {projects} project(s) to\n      {out}")
    return out


def verify_export(path: Path):
    with tempfile.TemporaryDirectory() as tmp:
        db = _extract_db(path, Path(tmp))
        counts(db)


def _extract_db(src: Path, tmp: Path) -> Path:
    """Database file from an export: .zip (OpenPrism export or any zip holding openprism.db), a .db file, or a folder."""
    if src.is_dir():
        cand = src / "openprism.db"
        if not cand.exists():
            fail(f"{src} does not contain openprism.db.")
        for sub in ("tls",):
            if (src / sub).is_dir():
                shutil.copytree(src / sub, tmp / sub)
        shutil.copy2(cand, tmp / "openprism.db")
        return tmp / "openprism.db"
    if zipfile.is_zipfile(src):
        with zipfile.ZipFile(src) as z:
            names = z.namelist()
            mf = next((n for n in names if n.rsplit("/", 1)[-1] == "manifest.json"), None)
            if mf:
                (tmp / "manifest.json").write_bytes(z.read(mf)[:100000])
            dbn = next((n for n in names if n.rsplit("/", 1)[-1] == "openprism.db"), None)
            if not dbn:
                fail(f"{src} is not an OpenPrism export (no openprism.db inside).")
            root = dbn[: -len("openprism.db")]
            for n in names:
                if not n.startswith(root) or n.endswith("/"):
                    continue
                rel = n[len(root):]
                if rel == "openprism.db" or rel.startswith("tls/"):
                    target = (tmp / rel).resolve()
                    if tmp.resolve() not in target.parents:
                        fail(f"{src} contains an unsafe path ({n}).")
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with z.open(n) as fsrc, open(target, "wb") as fdst:
                        shutil.copyfileobj(fsrc, fdst)
        return tmp / "openprism.db"
    with open(src, "rb") as f:
        if f.read(16) == b"SQLite format 3\x00":
            shutil.copy2(src, tmp / "openprism.db")
            return tmp / "openprism.db"
    fail(f"{src} is not a .zip export, an openprism.db file or a folder containing one.")


def check_compatible(tmp: Path, db: Path, src: Path):
    """The import only works when this code understands the data format; refuse clearly otherwise."""
    manifest = {}
    if (tmp / "manifest.json").exists():
        try:
            manifest = json.loads((tmp / "manifest.json").read_text())
        except ValueError:
            fail(f"{src}: the manifest inside the export is damaged.")
        if manifest.get("format") != EXPORT_FORMAT:
            fail(f"{src} was made by a different, incompatible OpenPrism version (format “{manifest.get('format')}”;\n"
                 f"  this version reads “{EXPORT_FORMAT}”). Import it with the OpenPrism version that created it.")
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
    except sqlite3.DatabaseError as e:
        fail(f"{src} cannot be read: {e}.")
    finally:
        con.close()
    missing = {"settings", "users", "projects"} - tables
    if missing:
        fail(f"{src} is not an OpenPrism server database (missing: {', '.join(sorted(missing))}).")
    theirs, ours = db_schema(db), this_schema()
    if theirs is None:
        fail(f"{src} has no schema version, so it is not known whether this version can read it. Nothing was changed.")
    if ours and theirs > ours:
        fail(f"{src} was made by a newer OpenPrism (data schema {theirs}; this version understands up to {ours}).\n"
             "  Install the newer version first, then import. Nothing was changed.")
    return manifest, theirs, ours


def find_backup(path: str) -> Path:
    """A backup .zip, a database file, or a folder: a folder holding openprism.db, or the newest backup .zip in it."""
    src = Path(path).expanduser().resolve()
    if not src.exists():
        fail(f"{src} does not exist.")
    if src.is_dir() and not (src / "openprism.db").exists():
        zips = sorted(src.glob("openprism-backup-*.zip")) or sorted(src.glob("*.zip"), key=lambda p: p.stat().st_mtime)
        if not zips:
            fail(f"No backup (.zip) was found in {src}.")
        src = zips[-1]
    return src


def load_backup(src: Path, tmp: Path):
    """Open and check a backup without touching the current data. Returns (db, users, projects, manifest)."""
    db = _extract_db(src, tmp)
    manifest, theirs, ours = check_compatible(tmp, db, src)
    try:
        users, projects = counts(db)
    except (sqlite3.DatabaseError, ValueError) as e:
        fail(f"{src} cannot be loaded: {e}. Nothing was deleted.")
    say(f"  Backup {src}\n      {users} account(s), {projects} project(s)"
        + (f", saved {manifest.get('created')} by OpenPrism {manifest.get('openprism')}" if manifest else "")
        + f" — data format {theirs} (this version reads up to {ours}): OK")
    return db, users, projects


def install_backup(db: Path, tmp: Path):
    DATA.mkdir(parents=True, exist_ok=True)
    shutil.copy2(db, DB)
    try:
        os.chmod(DB, 0o600)
    except OSError:
        pass
    if (tmp / "tls").is_dir():
        shutil.copytree(tmp / "tls", DATA / "tls", dirs_exist_ok=True)


def reset():
    files = data_files()
    for p in files:
        shutil.rmtree(p) if p.is_dir() and not p.is_symlink() else p.unlink()
    say("  All current data were deleted." if files else "  There were no data to delete.")


def selfpurge(output: str | None):
    if not output:
        fail("Say where to save the backup first, e.g.\n"
             "      python3 start.py --selfpurge --selfpurge --selfpurge --output ~/openprism-backups")
    folder = Path(output).expanduser().resolve()
    if folder == HERE or HERE in folder.parents:
        fail(f"The backup folder must be outside {HERE}, because that folder is about to be deleted.")
    if folder.exists() and not folder.is_dir():
        fail(f"{folder} is a file; give a folder for the backup.")
    folder.mkdir(parents=True, exist_ok=True)
    saved = None
    if DB.exists():
        saved = export(str(folder))
    elif data_files():
        fail(f"{DATA} contains files but no database; nothing was deleted. Check the folder by hand.")
    else:
        say("  There are no data to back up.")
    if DATA.exists() and HERE not in DATA.parents and DATA != HERE:
        shutil.rmtree(DATA)
    say(f"  Deleting {HERE} (code, Python environment, data, certificate, logs) …")
    os.chdir(HERE.parent)
    shutil.rmtree(HERE)
    say("  Done — this OpenPrism server was removed completely.\n")
    say("  To install the new version: git clone the repository, go to its Prism_Open_source folder and run")
    say(f"      python3 start.py --reset --reset --reset --loadolddata {saved or folder}" if saved else "      python3 start.py")
    say("")


# ---------------------------------------------------------------- main
def _version():
    try:
        txt = (HERE / "openprism" / "__init__.py").read_text()
        return txt.split('__version__ = "', 1)[1].split('"', 1)[0]
    except (OSError, IndexError):
        return "?"


def _value(args, flag):
    """Value after a flag (removed from args)."""
    if flag in args:
        i = args.index(flag)
        if i + 1 >= len(args) or args[i + 1].startswith("--"):
            fail(f"{flag} needs a path.")
        v = args[i + 1]
        del args[i:i + 2]
        return v
    return None


COMMANDS = """
      python3 start.py                                                   start the server
      python3 start.py --reset --reset --reset                           delete all data, start fresh
      python3 start.py --reset --reset --reset --loadolddata BACKUP      delete all data, load a backup, start
      python3 start.py --selfpurge --selfpurge --selfpurge --output DIR  back up to DIR, then remove this installation"""


def main():
    if sys.version_info < (3, 10):
        fail(f"OpenPrism needs Python 3.10 or newer (this is {sys.version.split()[0]}).")
    args = sys.argv[1:]
    if args[:1] in (["-h"], ["--help"]):
        say(__doc__)
        return
    old = [a for a in ("--update", "--preserve", "--export", "--import", "--input", "--data", "--yes") if a in args]
    if old:
        fail(f"{' '.join(old)} no longer exist{'s' if len(old) == 1 else ''}. The commands are:" + COMMANDS)
    n_purge, n_reset = args.count("--selfpurge"), args.count("--reset")
    output, olddata = _value(args, "--output"), _value(args, "--loadolddata")
    args = [a for a in args if a not in ("--selfpurge", "--reset")]

    if n_purge:
        if n_purge < 3:
            fail("--selfpurge removes this whole installation. Give it three times to confirm:\n"
                 "      python3 start.py --selfpurge --selfpurge --selfpurge --output ~/openprism-backups")
        if n_reset or olddata or args:
            fail("--selfpurge cannot be combined with other options." + COMMANDS)
        require_stopped("purge")
        selfpurge(output)
        return
    if output:
        fail("--output is only used with --selfpurge." + COMMANDS)
    if olddata and n_reset < 3:
        fail("--loadolddata replaces all current data, so it needs  --reset --reset --reset :\n"
             f"      python3 start.py --reset --reset --reset --loadolddata {olddata}")
    if n_reset:
        if n_reset < 3:
            fail("--reset deletes all accounts and projects. Give it three times to confirm:\n"
                 "      python3 start.py --reset --reset --reset")
        require_stopped("reset")
        if olddata:
            src = find_backup(olddata)
            with tempfile.TemporaryDirectory() as tmp:
                tmp = Path(tmp)
                db, users, projects = load_backup(src, tmp)      # checked before anything is deleted
                reset()
                install_backup(db, tmp)
            say(f"  Loaded {users} account(s) and {projects} project(s) from the backup.")
        else:
            reset()
            say("  Starting fresh: the first start runs the setup in the administration console.")
    if os.environ.get("OPENPRISM_NO_START") == "1":      # used by the automated tests
        return
    if not venv_python().exists():
        create_venv()
    install_requirements()
    os.chdir(HERE)
    cmd = [str(venv_python()), str(HERE / "run.py"), *args]
    if WIN:
        sys.exit(subprocess.call(cmd))
    os.execv(cmd[0], cmd)


if __name__ == "__main__":
    main()
