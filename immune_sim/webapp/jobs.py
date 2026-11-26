#!/usr/bin/env python3
"""Background jobs for the immune_sim GUI.

Every long operation in this project is an existing, validated CLI script.
The GUI runs those scripts as subprocesses rather than reimplementing them,
so the browser cannot drift from the pipeline the README documents and
VALIDATION.md measures. A job is therefore: an allowlisted argv, a worker
thread, a live log, a parsed progress line, and a cancel that kills the whole
process group (so calibrate.py's worker pool dies with it).

Nothing here accepts a command from the browser. The browser sends a task
NAME; tasks.py turns it into argv. That is the only reason this is safe to
bind to a port at all, and it is why `args` are validated per task instead of
being passed through.
"""
import os
import re
import signal
import subprocess
import threading
import time
import uuid
from collections import deque


class Job:
    """One run of one task. Thread-safe; the browser polls snapshot()."""

    def __init__(self, task, label, argv, cwd, produces=()):
        self.id = uuid.uuid4().hex[:12]
        self.task = task
        self.label = label
        self.argv = list(argv)
        self.cwd = str(cwd)
        self.produces = list(produces)
        self.state = "queued"          # queued|running|done|failed|cancelled
        self.created = time.time()
        self.started = None
        self.finished = None
        self.returncode = None
        self.error = None
        self.progress = None           # 0..1 when the task reports it
        self.note = ""                 # one line of parsed status
        self._log = deque(maxlen=600)
        self._lock = threading.RLock()
        self._proc = None
        self._cancelled = False

    # -- log -----------------------------------------------------------------
    def log(self, line):
        with self._lock:
            self._log.append(line.rstrip("\n"))

    def tail(self, n=200):
        with self._lock:
            return list(self._log)[-n:]

    # -- progress ------------------------------------------------------------
    def set_progress(self, frac, note=None):
        with self._lock:
            if frac is not None:
                self.progress = max(0.0, min(1.0, float(frac)))
            if note is not None:
                self.note = note[:300]

    # -- lifecycle -----------------------------------------------------------
    def run(self, parse=None):
        """Run argv to completion, streaming stdout into the log.

        `parse(job, line)` may set progress/note from each output line.
        """
        with self._lock:
            if self._cancelled:
                self.state = "cancelled"
                self.finished = time.time()
                return
            self.state = "running"
            self.started = time.time()
        self.log(f"$ {' '.join(self.argv)}")
        env = dict(os.environ, PYTHONUNBUFFERED="1")
        try:
            # start_new_session so cancel can signal the whole group
            proc = subprocess.Popen(
                self.argv, cwd=self.cwd, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1,
                start_new_session=True)
        except OSError as e:
            with self._lock:
                self.state = "failed"
                self.error = str(e)
                self.finished = time.time()
            self.log(f"could not start: {e}")
            return
        with self._lock:
            self._proc = proc
        try:
            for line in proc.stdout:
                self.log(line)
                if parse is not None:
                    try:
                        parse(self, line)
                    except Exception:
                        pass            # a progress parser must never kill a job
        finally:
            proc.stdout.close()
            rc = proc.wait()
            with self._lock:
                self.returncode = rc
                self.finished = time.time()
                if self._cancelled:
                    self.state = "cancelled"
                elif rc == 0:
                    self.state = "done"
                    self.progress = 1.0
                else:
                    self.state = "failed"
                    self.error = f"exited with code {rc}"
        self.log(f"[{self.state}] exit code {rc} after {self.elapsed:.0f} s")

    def cancel(self):
        with self._lock:
            self._cancelled = True
            proc = self._proc
            if self.state in ("done", "failed", "cancelled"):
                return False
        if proc is not None and proc.poll() is None:
            for sig in (signal.SIGTERM, signal.SIGKILL):
                try:
                    os.killpg(os.getpgid(proc.pid), sig)
                except (ProcessLookupError, PermissionError):
                    break
                for _ in range(20):
                    if proc.poll() is not None:
                        return True
                    time.sleep(0.1)
        return True

    @property
    def elapsed(self):
        end = self.finished or time.time()
        return end - (self.started or self.created)

    def snapshot(self, log_lines=200):
        with self._lock:
            return {
                "id": self.id, "task": self.task, "label": self.label,
                "state": self.state, "progress": self.progress,
                "note": self.note, "elapsed": round(self.elapsed, 1),
                "returncode": self.returncode, "error": self.error,
                "argv": self.argv, "produces": self.produces,
                "created": self.created,
                "log": self.tail(log_lines),
            }


class JobManager:
    """Runs one job per task at a time; keeps the recent history in memory.

    Serialising per task matters: two calibrate.py runs writing the same
    output file would interleave checkpoints, and two harvests would fight
    over the same CSV. Different tasks may run concurrently.
    """

    def __init__(self, keep=40):
        self._jobs = {}
        self._order = deque(maxlen=keep)
        self._lock = threading.RLock()

    def running(self, task=None):
        with self._lock:
            return [j for j in self._jobs.values()
                    if j.state in ("queued", "running")
                    and (task is None or j.task == task)]

    def submit(self, task, label, argv, cwd, produces=(), parse=None):
        with self._lock:
            busy = self.running(task)
            if busy:
                raise RuntimeError(
                    f"{label} is already running (job {busy[0].id}). "
                    "Wait for it or cancel it first.")
            job = Job(task, label, argv, cwd, produces)
            self._jobs[job.id] = job
            self._order.append(job.id)
            for stale in list(self._jobs):
                if stale not in self._order:
                    self._jobs.pop(stale, None)
        threading.Thread(target=job.run, args=(parse,),
                         name=f"job-{job.task}", daemon=True).start()
        return job

    def get(self, job_id):
        with self._lock:
            return self._jobs.get(job_id)

    def recent(self, n=20):
        with self._lock:
            ids = list(self._order)[-n:][::-1]
            return [self._jobs[i] for i in ids if i in self._jobs]


# ---------------------------------------------------------------------------
# progress parsers for the specific scripts this project ships
# ---------------------------------------------------------------------------
_CAL_EVAL = re.compile(r"^\s*eval\s+(\d+)")
_CAL_BEST = re.compile(r"best\s+([0-9.eE+-]+)")
_CAL_SCORE = re.compile(r"score\s+([0-9.eE+-]+)")


def calibrate_parser(total_iters):
    """calibrate.py prints `eval N ... best S`; turn that into a bar."""
    def parse(job, line):
        m = _CAL_EVAL.search(line)
        if m:
            n = int(m.group(1))
            best = _CAL_BEST.search(line) or _CAL_SCORE.search(line)
            note = f"evaluation {n:,} of {total_iters:,}"
            if best:
                note += f"  ·  best score {best.group(1)}"
            job.set_progress(n / max(total_iters, 1), note)
        elif "baseline" in line.lower():
            job.set_progress(None, line.strip()[:200])
    return parse


def counting_parser(markers):
    """Generic: advance one step each time an expected marker appears.

    `markers` is a list of (substring, label). Progress is the fraction of
    markers seen, which is honest about ordering but not about duration: a
    bar that moves in equal steps for unequal work is marked approximate in
    the GUI rather than dressed up as a time estimate.
    """
    def parse(job, line):
        low = line.lower()
        seen = getattr(job, "_seen", None)
        if seen is None:
            seen = job._seen = set()
        for i, (sub, label) in enumerate(markers):
            if sub.lower() in low and i not in seen:
                seen.add(i)
                job.set_progress(len(seen) / len(markers), label)
                break
        else:
            s = line.strip()
            if s and not s.startswith("$"):
                job.set_progress(None, s[:200])
    return parse
