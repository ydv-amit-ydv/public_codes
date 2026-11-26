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
        self.stages = []               # [{name,label,state,done,total,msg}]
        self._stage = None             # name of the running stage
        self._plan_mode = False        # True while run_plan drives the job
        self._progress_at = None       # when `progress` last actually moved
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

    # -- stages ---------------------------------------------------------------
    def add_stages(self, specs):
        """Declare the steps up front so the browser can show the whole plan.

        `specs` is [(name, label)]. Everything is 'pending' until it starts,
        which is the point: a run that is going to take an hour should say
        what it is going to do before it does it.
        """
        with self._lock:
            self.stages = [{"name": n, "label": l, "state": "pending",
                            "done": 0, "total": None, "msg": ""}
                           for n, l in specs]

    def stage_begin(self, name, msg=""):
        with self._lock:
            self._stage = name
            for st in self.stages:
                if st["name"] == name:
                    st["state"] = "running"
                    st["msg"] = msg

    def stage_update(self, done=None, total=None, msg=None, name=None):
        with self._lock:
            for st in self.stages:
                if st["name"] == (name or self._stage):
                    if done is not None:
                        st["done"] = done
                    if total is not None:
                        st["total"] = total
                    if msg is not None:
                        st["msg"] = msg[:300]

    def stage_end(self, name=None, state="done", msg=None):
        with self._lock:
            for st in self.stages:
                if st["name"] == (name or self._stage):
                    st["state"] = state
                    if msg is not None:
                        st["msg"] = msg[:300]
                    if st["total"]:
                        st["done"] = st["total"]

    def _set_overall(self, f):
        """Set overall progress, remembering when it last actually moved."""
        with self._lock:
            if f is None:
                return
            if self.progress is None or abs(f - self.progress) > 1e-9:
                self._progress_at = time.time()
            self.progress = f

    def _stage_fraction(self):
        """Overall progress from the stages, when there is more than one."""
        with self._lock:
            if len(self.stages) < 2:
                return None
            total = len(self.stages)
            whole = sum(1 for st in self.stages if st["state"] == "done")
            part = 0.0
            for st in self.stages:
                if st["state"] == "running" and st["total"]:
                    part = min(1.0, st["done"] / st["total"])
            return min(1.0, (whole + part) / total)

    # -- progress ------------------------------------------------------------
    def set_progress(self, frac, note=None):
        """Record progress reported by the running script.

        Inside a plan this feeds the CURRENT STEP, not the overall bar: a
        step's own parser knows how far that step has got and nothing about
        the plan. Letting it write the overall bar made a five-step plan read
        100% while its second step was still running.
        """
        with self._lock:
            plan = self._plan_mode and self._stage is not None
            if note is not None:
                self.note = note[:300]
        if frac is not None:
            f = max(0.0, min(1.0, float(frac)))
            if plan:
                self.stage_update(done=int(round(f * 1000)), total=1000)
                self._set_overall(self._stage_fraction())
            else:
                self._set_overall(f)

    # -- lifecycle -----------------------------------------------------------
    def run_plan(self, steps):
        """Run several commands in order, as one job with one stage each.

        `steps` is [(name, label, argv, parser)]. A failing step stops the
        plan: the rest stay 'pending' rather than being marked skipped,
        because what the user needs to see is where it stopped. Cancelling
        stops the plan too, not only the step that is running.
        """
        self.add_stages([(n, l) for n, l, _a, _p in steps])
        with self._lock:
            self._plan_mode = True
            self.state = "running"
            self.started = time.time()
        ok = True
        for name, label, argv, parser in steps:
            with self._lock:
                if self._cancelled:
                    ok = False
                    break
            self.stage_begin(name)
            self.log(f"\n=== {label} ===")
            rc = self._spawn(argv, parser, stage=name)
            if rc != 0:
                self.stage_end(name, "failed",
                               f"exit code {rc} - see the log")
                ok = False
                break
            self.stage_end(name, "done")
            with self._lock:
                self._stage = None
                self.note = ""
            self._set_overall(self._stage_fraction())
        with self._lock:
            self.finished = time.time()
            if self._cancelled:
                self.state = "cancelled"
            elif ok:
                self.state = "done"
                self.progress = 1.0
            else:
                self.state = "failed"
                self.error = self.error or "a step failed; the log says which"
        self.log(f"[{self.state}] after {self.elapsed:.0f} s")

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
        rc = self._spawn(self.argv, parse)
        with self._lock:
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

    def _spawn(self, argv, parse=None, stage=None):
        """Start one command, stream its output, return its exit code."""
        self.log(f"$ {' '.join(argv)}")
        env = dict(os.environ, PYTHONUNBUFFERED="1")
        try:
            # start_new_session so cancel can signal the whole group
            proc = subprocess.Popen(
                argv, cwd=self.cwd, env=env, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, text=True, bufsize=1,
                start_new_session=True)
        except OSError as e:
            with self._lock:
                self.error = str(e)
            self.log(f"could not start: {e}")
            return 127
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
                if stage is not None and self.note:
                    self.stage_update(msg=self.note, name=stage)
        finally:
            proc.stdout.close()
            rc = proc.wait()
            with self._lock:
                self.returncode = rc
        return rc

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

    def _eta(self):
        """Seconds left, extrapolated from progress. None when meaningless.

        Two guards, both learned from watching it be wrong:

        - Nothing before 5% in. The first few percent of an annealed search or
          a per-fold retrain are not representative.
        - Nothing while progress is STALE. A step that reports no progress of
          its own freezes the fraction, and extrapolating a frozen fraction
          against a growing elapsed time produces an estimate that climbs for
          as long as you watch it - it read "about 60 s left" and then "about
          25 min left" on the same unchanged step. No estimate is better than
          one that grows.
        """
        p = self.progress
        if (self.state != "running" or p is None or p < 0.05 or p >= 1.0
                or not self.started):
            return None
        with self._lock:
            moved = self._progress_at
        if moved is None or time.time() - moved > 90:
            return None
        return round(self.elapsed * (1.0 - p) / p)

    @property
    def elapsed(self):
        end = self.finished or time.time()
        return end - (self.started or self.created)

    def snapshot(self, log_lines=200):
        with self._lock:
            return {
                "id": self.id, "task": self.task, "label": self.label,
                "state": self.state, "progress": self.progress,
                "stages": list(self.stages), "eta": self._eta(),
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

    def submit_plan(self, task, label, steps, cwd, produces=()):
        """Queue a multi-step plan as one job (see Job.run_plan)."""
        with self._lock:
            busy = self.running(task)
            if busy:
                raise RuntimeError(
                    f"{label} is already running (job {busy[0].id}). "
                    "Wait for it or cancel it first.")
            job = Job(task, label, [], cwd, produces)
            self._jobs[job.id] = job
            self._order.append(job.id)
        threading.Thread(target=job.run_plan, args=(steps,),
                         name=f"plan-{job.task}", daemon=True).start()
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
