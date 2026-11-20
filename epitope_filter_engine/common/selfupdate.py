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
import time

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


def apply():
    """Fast-forward to the remote. Returns a dict with ok, message and what changed."""
    st = status()
    if st.get("error"):
        return {"ok": False, "error": st["error"], "status": st}
    if not st["available"]:
        return {"ok": False, "error": "Already up to date.", "status": st}
    if st["dirty"]:
        return {"ok": False, "status": st,
                "error": ("This folder has local changes, so an update would overwrite them: "
                          + ", ".join(st["dirty"][:5]) + ". Commit, stash or undo them first.")}
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
