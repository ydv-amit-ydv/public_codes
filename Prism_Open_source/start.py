#!/usr/bin/env python3
"""OpenPrism server launcher and manager — needs only Python 3.10+.

Start
    python3 start.py                                start the server (installs the packages on the first run)
    python3 start.py --port 8443                    other options are passed to the server (python3 run.py --help)

Manage the data (stop the server with Ctrl+C first)
    python3 start.py --export --output FILE.zip     save all accounts, projects, shares and the certificate to FILE.zip
    python3 start.py --import --input FILE.zip      replace ALL current data with an export (asks you to type IMPORT;
                                                    only exports this version can read are accepted)
    python3 start.py --reset --reset --reset        delete ALL data and start fresh (the flag three times = confirmation)
    python3 start.py --selfpurge --export FILE.zip  export everything, check the export, then delete everything that
                                                    belongs to this server (code, environment, data) — re-clone later

Moving to a new version:  export → delete or purge the old copy → git clone the new one → --import → start.
The first run creates a private Python environment in .venv/ next to this file (nothing is installed system-wide).
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
MANAGE = ("--export", "--import", "--reset", "--selfpurge", "--output", "--input")


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
        fail("Say where to save the export, e.g.\n      python3 start.py --export --output ~/openprism-backup.zip")
    out = Path(output).expanduser().resolve()
    if out.is_dir():
        out = out / f"openprism-export-{datetime.datetime.now():%Y%m%d-%H%M%S}.zip"
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
    say(f"  Exported {users} account(s) and {projects} project(s) to\n      {out}")
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


def do_import(inp: str | None, assume_yes: bool = False):
    if not inp:
        fail("Say which export to restore, e.g.\n      python3 start.py --import --input ~/openprism-backup.zip")
    src = Path(inp).expanduser().resolve()
    if not src.exists():
        fail(f"{src} does not exist.")
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        db = _extract_db(src, tmp)
        manifest, theirs, ours = check_compatible(tmp, db, src)
        try:
            users, projects = counts(db)
        except (sqlite3.DatabaseError, ValueError) as e:
            fail(f"{src} cannot be restored: {e}.")
        say(f"\n  Import from {src}")
        say(f"      {users} account(s), {projects} project(s)"
            + (f", exported {manifest.get('created')} by OpenPrism {manifest.get('openprism')}" if manifest else "")
            + f" (data schema {theirs}; this version: {ours})")
        if DB.exists():
            try:
                cu, cp = counts(DB)
                say(f"  This REPLACES everything currently on this server: {cu} account(s), {cp} project(s),\n"
                    "  shares, settings and the certificate (a safety copy of the current database is kept in projects/backups/).")
            except (sqlite3.DatabaseError, ValueError):
                say("  This REPLACES the current server data (a safety copy is kept in projects/backups/).")
        elif data_files():
            say("  This replaces the files currently in the data folder.")
        if not assume_yes:
            try:
                ans = input("  Type IMPORT to continue: ").strip()
            except EOFError:
                ans = ""
            if ans != "IMPORT":
                fail("Cancelled — nothing was changed.")
        if data_files():
            if DB.exists():
                safety = export_safety_copy()
                say(f"  The current data were saved first to\n      {safety}")
            for p in data_files():
                if p.name != "backups":
                    shutil.rmtree(p) if p.is_dir() and not p.is_symlink() else p.unlink()
        DATA.mkdir(parents=True, exist_ok=True)
        shutil.copy2(db, DB)
        try:
            os.chmod(DB, 0o600)
        except OSError:
            pass
        if (tmp / "tls").is_dir():
            shutil.copytree(tmp / "tls", DATA / "tls", dirs_exist_ok=True)
    say(f"  Restored {users} account(s) and {projects} project(s).")


def export_safety_copy() -> Path:
    bdir = DATA / "backups"
    bdir.mkdir(parents=True, exist_ok=True)
    out = bdir / f"before-import-{datetime.datetime.now():%Y%m%d-%H%M%S}.db"
    copy_db(DB, out)
    return out


# ---------------------------------------------------------------- reset / purge
def reset():
    files = data_files()
    if not files:
        say(f"  Nothing to reset: {DATA} is already empty.")
        return
    for p in files:
        shutil.rmtree(p) if p.is_dir() and not p.is_symlink() else p.unlink()
    say(f"  All server data in {DATA} were deleted (accounts, projects, shares, certificate, logs, backups).\n"
        "  The next start runs the first-time setup again.")


def selfpurge(output: str | None):
    if not output:
        fail("--selfpurge needs an export first, so nothing is lost:\n"
             "      python3 start.py --selfpurge --export --output ~/openprism-backup.zip")
    out = Path(output).expanduser().resolve()
    if HERE in out.parents or out == HERE:
        fail(f"The export must be saved outside {HERE}, because that folder is about to be deleted.")
    exported = None
    if DB.exists():
        exported = export(output)
    elif data_files():
        fail(f"{DATA} contains files but no database; nothing was deleted. Check the folder by hand.")
    else:
        say("  There are no server data to export.")
    repo = next((p for p in [HERE, *HERE.parents] if (p / ".git").exists()), None)
    if DATA.exists() and HERE not in DATA.parents and DATA != HERE:
        say(f"  Deleting the data folder {DATA} …")
        shutil.rmtree(DATA)
    say(f"  Deleting {HERE} (code, Python environment, data, certificate, logs, backups) …")
    os.chdir(HERE.parent)
    shutil.rmtree(HERE)
    say("  Done. Everything belonging to this OpenPrism server was removed.")
    say("\n  To set it up again:")
    if repo and repo != HERE:
        say(f"      cd {repo} && git checkout -- {HERE.relative_to(repo)}      (restore this folder from git), or delete")
        say(f"      {repo}  and clone the repository again;")
    else:
        say("      git clone <repository> and cd into its Prism_Open_source folder;")
    if exported:
        say(f"      python3 start.py --import --input {exported}")
    say("      python3 start.py\n")


# ---------------------------------------------------------------- main
def _version():
    try:
        txt = (HERE / "openprism" / "__init__.py").read_text()
        return txt.split('__version__ = "', 1)[1].split('"', 1)[0]
    except (OSError, IndexError):
        return "?"


def _value(args, flag):
    """Value after --output / --input, removed from args."""
    if flag in args:
        i = args.index(flag)
        if i + 1 >= len(args) or args[i + 1].startswith("--"):
            fail(f"{flag} needs a path, e.g. {flag} ~/openprism-backup.zip")
        v = args[i + 1]
        del args[i:i + 2]
        return v
    return None


def main():
    if sys.version_info < (3, 10):
        fail(f"OpenPrism needs Python 3.10 or newer (this is {sys.version.split()[0]}).")
    args = sys.argv[1:]
    if args[:1] in (["-h"], ["--help"]):
        say(__doc__)
        return
    if "--update" in args or "--preserve" in args:
        fail("--update and --preserve were replaced by export / import:\n"
             "      python3 start.py --export --output ~/openprism-backup.zip      (save everything)\n"
             "      … get the new version (git clone / git pull / ZIP) …\n"
             "      python3 start.py --import --input ~/openprism-backup.zip       (restore it)\n"
             "  Data in the projects/ folder are never touched by git, so after  git pull  a plain  python3 start.py  keeps them.")
    # the path may also follow the command directly:  --export FILE.zip  /  --import FILE.zip
    for flag, opt in (("--export", "--output"), ("--import", "--input")):
        if flag in args:
            i = args.index(flag)
            if i + 1 < len(args) and not args[i + 1].startswith("--") and opt not in args:
                args[i + 1:i + 1] = [opt]
    output, inp = _value(args, "--output"), _value(args, "--input")
    n_reset = args.count("--reset")
    actions = [a for a in ("--export", "--import", "--selfpurge") if a in args] + (["--reset"] if n_reset else [])
    if len(actions) > 1 and not (set(actions) == {"--selfpurge", "--export"}):
        fail(f"Use one management command at a time (got {' '.join(actions)}).")
    if actions:
        extra = [a for a in args if a not in ("--export", "--import", "--selfpurge", "--reset", "--yes")]
        if extra:
            fail(f"Management commands do not start the server; remove {' '.join(extra)}.")
    if "--selfpurge" in args:
        if "--export" not in args:
            fail("--selfpurge always exports first:\n      python3 start.py --selfpurge --export --output ~/openprism-backup.zip")
        require_stopped("purge")
        selfpurge(output)
        return
    if "--export" in args:
        require_stopped("export")       # a running server keeps writing; export after stopping it
        export(output)
        return
    if "--import" in args:
        require_stopped("import")
        do_import(inp, "--yes" in args)
        say("  Start the server with   python3 start.py")
        return
    if n_reset:
        if n_reset < 3:
            fail(f"Resetting deletes every account and project. To confirm, give the flag three times:\n"
                 f"      python3 start.py --reset --reset --reset\n"
                 f"  (you gave it {n_reset} time{'s' if n_reset > 1 else ''}). Tip: export first with  python3 start.py --export --output FILE.zip")
        require_stopped("reset")
        reset()
        return
    if output or inp:
        fail("--output belongs to --export and --input to --import.")
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
