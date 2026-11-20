"""Is a newer engine available, and install it.

The engine is normally a git clone of its public repository, so an update is a fast-forward
pull. Only the code folder changes: runs, settings, caches and the learned timings live in
the workspace, which git ignores, so **updating never changes a finished run**. A run also
records the version and commit it was produced with (see common/runs.py), so an old result
stays labelled with the code that made it.

When the engine was installed by copying the files instead of cloning, there is nothing to
pull; `status()` says so and gives the command to do it by hand.
"""
import os
import subprocess
import threading
import time
import re

_ANSI = re.compile("\x1b\\[[0-9;]*m")       # the installer colours its headings

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO = "https://github.com/ydv-amit-ydv/public_codes"
TIMEOUT = 120                       # a fetch over a slow link, not a long build


def _git(*args, timeout=TIMEOUT):
    """Run git in the engine folder. Never raises: returns (ok, stdout, stderr)."""
    try:
        r = subprocess.run(["git", "-C", ROOT, *args], capture_output=True, text=True,
                           timeout=timeout, stdin=subprocess.DEVNULL)
        return r.returncode == 0, r.stdout.strip(), r.stderr.strip()
    except FileNotFoundError:
        return False, "", "git is not installed"
    except subprocess.TimeoutExpired:
        return False, "", f"git took longer than {timeout}s"
    except OSError as e:
        return False, "", str(e)


def version():
    try:
        with open(os.path.join(ROOT, "VERSION")) as f:
            return f.read().strip()
    except OSError:
        return "unknown"


def is_git():
    ok, out, _ = _git("rev-parse", "--is-inside-work-tree", timeout=20)
    return ok and out == "true"


def commit(short=True):
    """The installed commit, or None when this is not a clone."""
    args = ("rev-parse", "--short=12", "HEAD") if short else ("rev-parse", "HEAD")
    ok, out, _ = _git(*args, timeout=20)
    return out if (ok and out) else None


def _branch():
    ok, out, _ = _git("rev-parse", "--abbrev-ref", "HEAD", timeout=20)
    return out if ok and out != "HEAD" else "main"


def _count(rev):
    ok, out, _ = _git("rev-list", "--count", rev, timeout=30)
    try:
        return int(out) if ok else None
    except ValueError:
        return None


def status(check_remote=True):
    """What is installed, and whether something newer is on the remote.

    `check_remote=False` skips the network, for the page's first paint.
    """
    d = {"version": version(), "root": ROOT, "is_git": False, "available": False,
         "behind": None, "ahead": None, "branch": None, "commit": None, "remote": None,
         "dirty": [], "checked": None, "error": None,
         "manual": f"git -C {ROOT} pull --ff-only"}
    if not is_git():
        d["error"] = ("This copy of the engine is not a git clone, so it cannot update itself. "
                      f"Clone it once ({REPO}) and run it from the clone, or copy the new files over this folder.")
        return d
    d["is_git"] = True
    d["branch"] = _branch()
    ok, out, _ = _git("rev-parse", "HEAD", timeout=20)
    d["commit"] = out if ok else None
    ok, out, _ = _git("config", "--get", "remote.origin.url", timeout=20)
    d["remote"] = out if ok else None
    # Only *tracked* edits can be overwritten by a fast-forward, so untracked files are not
    # dirt: __pycache__, a .venv or a stray CSV in the folder must not block an update. If an
    # untracked file really would be overwritten, git itself refuses the merge and says which.
    ok, out, _ = _git("status", "--porcelain", "--untracked-files=no", timeout=30)
    if ok and out:
        d["dirty"] = [l.strip() for l in out.splitlines()[:20]]
        # porcelain is "XY path", but the output has been stripped, so the first line has lost
        # its leading space and a fixed slice would eat the first character of the path
        d["dirty_paths"] = [re.sub(r"^\s*\S{1,2}\s+", "", l).strip('"').split(" -> ")[-1]
                            for l in out.splitlines()]
    if not check_remote:
        return d

    # fetch writes FETCH_HEAD and the remote-tracking ref only: the working tree is untouched
    ok, _, err = _git("fetch", "--quiet", "origin", d["branch"])
    d["checked"] = time.time()
    if not ok:
        d["error"] = f"Could not reach the remote: {err or 'unknown error'}"
        return d
    ok, head, _ = _git("rev-parse", "FETCH_HEAD", timeout=20)
    if not ok:
        d["error"] = "Could not read the fetched revision."
        return d
    d["latest"] = head
    if head == d["commit"]:
        return d                                       # already current

    # Only the local changes the update would actually touch can block it. This repository
    # holds several projects, so an edit under another one is none of the engine's business:
    # git merges around it and keeps it. Blocking on every change in the tree meant a deleted
    # file in a sibling project stopped the engine updating, which is what happened.
    ok, names, _ = _git("diff", "--name-only", "HEAD..FETCH_HEAD", timeout=30)
    incoming = set(names.splitlines()) if ok and names else set()
    d["incoming"] = len(incoming)
    d["blocking"] = sorted(set(d.get("dirty_paths") or []) & incoming)
    behind, ahead = _count(f"HEAD..FETCH_HEAD"), _count(f"FETCH_HEAD..HEAD")
    d["behind"], d["ahead"] = behind, ahead
    if behind is None:                                 # shallow clone: counting can fail
        d["available"] = True
        d["note"] = "The remote has a different revision (this is a shallow clone, so the number of commits is unknown)."
    elif behind > 0 and (ahead or 0) == 0:
        d["available"] = True
    elif behind > 0 and ahead:
        d["error"] = (f"This folder has {ahead} commit(s) of its own and is {behind} behind the remote, so it cannot "
                      f"fast-forward. Sort it out with git, or re-clone.")
    return d


# ---------------------------------------------------------------- repairing the predictors
#
# Updating the code never touches the predictors in the models folder: that is the point of
# keeping them apart. But a fix to the installer is only worth anything once the installer
# runs again, so the page can run it. It is the same install_tools.sh the command line uses,
# it is idempotent, and it repairs in place rather than downloading again.
_repair = {"running": False, "log": [], "started": None, "finished": None, "ok": None}
_repair_lock = threading.Lock()
MAX_LOG = 4000


def repair_status():
    with _repair_lock:
        return {"running": _repair["running"], "started": _repair["started"],
                "finished": _repair["finished"], "ok": _repair["ok"],
                "log": list(_repair["log"])[-400:], "lines": len(_repair["log"])}


def repair_start(models_dir=None, light=False):
    """Run install_tools.sh in the background. Returns at once; follow it with repair_status()."""
    script = os.path.join(ROOT, "install_tools.sh")
    if not os.path.exists(script):
        return {"ok": False, "error": "install_tools.sh is not in this folder."}
    with _repair_lock:
        if _repair["running"]:
            return {"ok": False, "error": "A repair is already running."}
        _repair.update({"running": True, "log": [], "started": time.time(), "finished": None, "ok": None})

    cmd = ["bash", script]
    if models_dir:
        cmd += ["--dir", models_dir]
    if not light:
        cmd += ["--iedb"]

    def work():
        rc = 1
        try:
            # ~/.local/bin holds uv, which the IEDB bundles need; a server started
            # outside start.sh may not have it on PATH
            local_bin = os.path.join(os.path.expanduser("~"), ".local", "bin")
            path = os.environ.get("PATH", "")
            if os.path.isdir(local_bin) and local_bin not in path.split(os.pathsep):
                path = local_bin + os.pathsep + path
            env = dict(os.environ, ENGINE_FROM_START="1", PYTHONUNBUFFERED="1", PATH=path)
            p = subprocess.Popen(cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                 stdin=subprocess.DEVNULL, text=True, bufsize=1, env=env)
            for line in p.stdout:
                line = _ANSI.sub("", line.rstrip("\n"))      # the installer colours its headings
                with _repair_lock:
                    if len(_repair["log"]) < MAX_LOG:
                        _repair["log"].append(line)
                    elif len(_repair["log"]) == MAX_LOG:
                        _repair["log"].append("... (log truncated; run it in a terminal to see everything)")
            rc = p.wait()
        except Exception as e:                           # noqa: BLE001 - always report, never crash the server
            with _repair_lock:
                _repair["log"].append(f"could not run the installer: {e}")
        finally:
            with _repair_lock:
                _repair.update({"running": False, "finished": time.time(), "ok": rc == 0})

    threading.Thread(target=work, daemon=True).start()
    return {"ok": True, "message": "Repair started."}


def apply():
    """Fast-forward to the remote. Returns a dict with ok, message and what changed."""
    st = status()
    if st.get("error"):
        return {"ok": False, "error": st["error"], "status": st}
    if not st["available"]:
        return {"ok": False, "error": "Already up to date.", "status": st}
    blocking = st.get("blocking") or []
    if blocking:
        return {"ok": False, "status": st,
                "error": ("The update changes files you have edited locally, so it would "
                          "overwrite them: " + ", ".join(blocking[:5])
                          + (f" and {len(blocking) - 5} more" if len(blocking) > 5 else "")
                          + ". Commit, stash or undo those, then update.")}
    before = st["commit"]
    ok, out, err = _git("merge", "--ff-only", "FETCH_HEAD", timeout=TIMEOUT)
    if not ok:
        return {"ok": False, "error": f"Update failed: {err or out}", "status": status(check_remote=False)}
    after = status(check_remote=False)
    files = []
    ok2, names, _ = _git("diff", "--name-only", f"{before}..HEAD", timeout=30)
    if ok2 and names:
        files = names.splitlines()
    return {"ok": True, "from": before, "to": after["commit"], "version": after["version"],
            "files": files[:200], "n_files": len(files), "message": out,
            "restart_needed": True, "status": after}
