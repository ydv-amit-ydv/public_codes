#!/usr/bin/env python3
"""OpenPrism launcher: only Python 3.10+ is needed.

    python3 start.py              install what is needed (first run only) and start the server
    python3 start.py --reset      delete ALL server data (accounts, projects, certificate, logs) and start fresh
    python3 start.py --update     pull the latest version from git and update the requirements
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


def update():
    if not (HERE.parent / ".git").exists() and not (HERE / ".git").exists():
        sys.exit("  This copy was not cloned with git; download the new version instead.")
    say("  Updating from git …")
    subprocess.run(["git", "-C", str(HERE), "pull", "--ff-only"], check=True)


def main():
    if sys.version_info < (3, 10):
        sys.exit(f"  OpenPrism needs Python 3.10 or newer (this is {sys.version.split()[0]}).")
    args = sys.argv[1:]
    if args[:1] in (["-h"], ["--help"]):
        say(__doc__)
        return
    yes = "--yes" in args
    if "--reset" in args:
        reset(yes)
        args = [a for a in args if a not in ("--reset", "--yes")]
        if "--no-start" in args:
            return
    if "--update" in args:
        update()
        args = [a for a in args if a != "--update"]
    args = [a for a in args if a != "--no-start"]
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
