#!/usr/bin/env python3
"""OpenPrism launcher: only Python 3.10+ is needed.

    python3 start.py              install what is needed (first run only) and start the server
    python3 start.py --reset      delete ALL server data (accounts, projects, certificate, logs) and start fresh
    python3 start.py --update --preserve   download the latest version, keep all accounts and projects (backup first)
    python3 start.py --update --reset      download the latest version and start with no accounts or projects
    python3 start.py --preserve            (copies without git, e.g. from a ZIP) back up the data, keep it, and start
    python3 start.py --port 8443  any other option is passed to the server (see python3 run.py --help)

The first run creates a private Python environment in .venv/ next to this file and installs the
packages from requirements.txt into it (nothing is installed system-wide). Later runs start at once;
the packages are reinstalled automatically only when requirements.txt changes.
"""
from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
VENV = HERE / ".venv"
REQ = HERE / "requirements.txt"
STAMP = VENV / ".requirements.sha256"
DATA = Path(os.environ.get("OPENPRISM_DATA") or HERE / "projects")
WIN = os.name == "nt"


def say(msg=""):
    print(msg, flush=True)


def venv_python() -> Path:
    return VENV / ("Scripts/python.exe" if WIN else "bin/python")


def create_venv():
    say("  Creating a private Python environment in .venv/ …")
    try:
        subprocess.run([sys.executable, "-m", "venv", str(VENV)], check=True)
    except subprocess.CalledProcessError:
        shutil.rmtree(VENV, ignore_errors=True)
        sys.exit("\n  Could not create the environment. On Debian/Ubuntu install the venv module first:\n"
                 "      sudo apt install python3-venv\n  then run  python3 start.py  again.\n")


def install_requirements():
    want = hashlib.sha256(REQ.read_bytes()).hexdigest()
    if STAMP.exists() and STAMP.read_text().strip() == want:
        return
    say("  Installing the required packages (first run only, about a minute) …")
    py = str(venv_python())
    subprocess.run([py, "-m", "pip", "install", "--quiet", "--upgrade", "pip"], check=False)
    r = subprocess.run([py, "-m", "pip", "install", "--quiet", "-r", str(REQ)])
    if r.returncode != 0:
        sys.exit("\n  Installing the packages failed (see the messages above). Check the internet connection and run\n"
                 "      python3 start.py\n  again.\n")
    STAMP.write_text(want)
    say("  Packages installed.")


def server_running() -> bool:
    pid_ports = []
    try:
        import socket
        for port in (5000, 47219):
            with socket.socket() as s:
                s.settimeout(0.3)
                if s.connect_ex(("127.0.0.1", port)) == 0:
                    pid_ports.append(port)
    except OSError:
        pass
    return bool(pid_ports)


def reset(assume_yes: bool):
    files = [p for p in DATA.iterdir() if p.name != ".gitkeep"] if DATA.exists() else []
    if not files:
        say(f"  Nothing to reset: {DATA} is already empty.")
        return
    say(f"\n  This permanently deletes ALL OpenPrism server data in\n      {DATA}\n"
        "  (every account, every project, shares, the HTTPS certificate and the logs).")
    if server_running():
        say("  ⚠ A server seems to be running on this computer — stop it first (Ctrl+C in its terminal).")
    if not assume_yes:
        try:
            ans = input("  Type RESET to confirm: ").strip()
        except EOFError:
            ans = ""
        if ans != "RESET":
            sys.exit("  Cancelled — nothing was deleted.")
    for p in files:
        shutil.rmtree(p) if p.is_dir() and not p.is_symlink() else p.unlink()
    say("  All server data deleted. The next start runs the first-time setup again.\n")


def backup():
    """Copy the database (accounts, encrypted projects, shares) to projects/backups/ before updating."""
    db = DATA / "openprism.db"
    if not db.exists():
        say("  No existing server data to back up.")
        return
    import datetime
    import sqlite3
    dest_dir = DATA / "backups"
    dest_dir.mkdir(exist_ok=True)
    dest = dest_dir / f"openprism-{datetime.datetime.now():%Y%m%d-%H%M%S}.db"
    src = sqlite3.connect(str(db))
    try:
        out = sqlite3.connect(str(dest))
        with out:
            src.backup(out)          # consistent copy even if the database is in WAL mode
        out.close()
    finally:
        src.close()
    old = sorted(dest_dir.glob("openprism-*.db"))
    for f in old[:-10]:              # keep the 10 most recent backups
        f.unlink()
    say(f"  Backup of all accounts and projects: {dest}")


def git_pull():
    if not (HERE.parent / ".git").exists() and not (HERE / ".git").exists():
        sys.exit("\n  This copy was not downloaded with git (for example from a ZIP file), so it cannot update itself.\n"
                 "  To update it by hand:\n"
                 "    1. Stop the server (Ctrl+C).\n"
                 "    2. Download the new version and copy its files over this folder. Keep the folder\n"
                 f"         {DATA}\n"
                 "       — it holds all accounts and projects (do not replace or delete it).\n"
                 "    3. Start with   python3 start.py --preserve   (keeps all data, makes a backup first)\n"
                 "       or with      python3 start.py --reset      (starts with no data).\n"
                 "  Tip: installing with  git clone  lets  python3 start.py --update --preserve  do all of this for you.")
    say("  Downloading the latest version (git pull) …")
    r = subprocess.run(["git", "-C", str(HERE), "pull", "--ff-only"])
    if r.returncode != 0:
        sys.exit("\n  git pull failed (see above). Your data were not changed. If you edited files in this folder, run\n"
                 "      git -C . status\n  to see them, or clone a fresh copy and copy the projects/ folder into it.")


UPDATE_HELP = """
  Choose what happens to the existing accounts and projects when updating:

      python3 start.py --update --preserve    keep all accounts and projects (a backup is made first)
      python3 start.py --update --reset       delete all accounts and projects and start fresh
"""


def main():
    if sys.version_info < (3, 10):
        sys.exit(f"  OpenPrism needs Python 3.10 or newer (this is {sys.version.split()[0]}).")
    args = sys.argv[1:]
    if args[:1] in (["-h"], ["--help"]):
        say(__doc__)
        return
    preserve, wipe, updating, updated = "--preserve" in args, "--reset" in args, "--update" in args, "--updated" in args
    if preserve and wipe:
        sys.exit("  Use either --preserve or --reset, not both." + UPDATE_HELP)
    if updating and not (preserve or wipe):
        sys.exit("  --update needs --preserve or --reset." + UPDATE_HELP)
    yes = "--yes" in args
    if updating and not updated:
        if not (HERE.parent / ".git").exists() and not (HERE / ".git").exists():
            git_pull()                    # explains how to update a copy without git, then stops
        if server_running():
            sys.exit("  OpenPrism seems to be running on this computer. Stop it first (Ctrl+C in its window), then update.")
        if wipe and not yes:
            say(f"\n  After updating, ALL accounts and projects in {DATA} will be deleted.")
            try:
                ans = input("  Type RESET to confirm: ").strip()
            except EOFError:
                ans = ""
            if ans != "RESET":
                sys.exit("  Cancelled — nothing was downloaded or deleted.")
            yes = True
        if preserve:
            backup()
        git_pull()
        # continue with the freshly downloaded launcher
        rest = [a for a in args if a not in ("--update", "--yes")] + ["--updated"] + (["--yes"] if yes else [])
        cmd = [sys.executable, str(HERE / "start.py"), *rest]
        if WIN:
            sys.exit(subprocess.call(cmd))
        os.execv(cmd[0], cmd)
    if wipe:
        reset(yes)
        if "--no-start" in args:
            return
    elif preserve:
        if not updated:
            backup()
        say(f"  Keeping all accounts and projects in {DATA}.")
    args = [a for a in args if a not in ("--reset", "--preserve", "--update", "--updated", "--yes", "--no-start")]
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
