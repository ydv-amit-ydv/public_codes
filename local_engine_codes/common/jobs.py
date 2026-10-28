"""Background jobs with per-stage progress, live log and cancel.

A job runs on a worker thread; the browser polls Job.snapshot(). External
programs are started through run_cmd(), which registers the process with the job
so Cancel can stop it (its whole process group, so BigMHC's data-loader workers
die too) and can stream output lines for progress parsing.
"""
import contextvars, os, shutil, signal, subprocess, threading, time, uuid
from concurrent.futures import ThreadPoolExecutor


class Cancelled(Exception):
    pass


class Job:
    def __init__(self, kind="pipeline"):
        self.id = uuid.uuid4().hex[:12]
        self.kind = kind
        self.state = "queued"           # queued | running | done | failed | cancelled
        self.created = time.time()
        self.started = self.finished = None
        self.stages = []                # [{name, label, weight, done, total, state, msg}]
        self.logs = []
        self.error = None
        self.result = None              # filled by the pipeline (summary dict + file paths)
        self.cancelled = False
        self._lock = threading.RLock()
        self._procs = set()
        self._peak = 0.0

    # --- progress ---------------------------------------------------------------
    def stage(self, name, label, weight=1.0, total=None):
        with self._lock:
            self.stages = [s for s in self.stages if s["name"] != name]
            self.stages.append({"name": name, "label": label, "weight": weight, "done": 0,
                                "total": total, "state": "pending", "msg": ""})

    def update(self, name, done=None, total=None, msg=None, state=None):
        with self._lock:
            for s in self.stages:
                if s["name"] == name:
                    if done is not None: s["done"] = done
                    if total is not None: s["total"] = total
                    if msg is not None: s["msg"] = msg
                    if state: s["state"] = state
                    elif s["state"] == "pending": s["state"] = "running"

    def finish_stage(self, name, msg=None):
        with self._lock:
            for s in self.stages:
                if s["name"] == name:
                    s["state"] = "done"
                    if s["total"]: s["done"] = s["total"]
                    if msg is not None: s["msg"] = msg

    def log(self, line):
        with self._lock:
            self.logs.append(line)
            del self.logs[:-400]

    def check(self):
        if self.cancelled:
            raise Cancelled()

    # --- processes ----------------------------------------------------------------
    def _register(self, p):
        with self._lock:
            self._procs.add(p)

    def _unregister(self, p):
        with self._lock:
            self._procs.discard(p)

    def cancel(self):
        with self._lock:
            self.cancelled = True
            procs = list(self._procs)
        for p in procs:
            try:
                os.killpg(p.pid, signal.SIGTERM)
            except (ProcessLookupError, PermissionError):
                pass

    # --- view ----------------------------------------------------------------------
    def fraction(self):
        with self._lock:
            tot = sum(s["weight"] for s in self.stages) or 1.0
            acc = 0.0
            for s in self.stages:
                if s["state"] == "done":
                    f = 1.0
                elif s["total"]:
                    f = min(1.0, s["done"] / s["total"])
                else:
                    f = 0.0
                acc += s["weight"] * f
            self._peak = max(self._peak, acc / tot)             # a progress bar must never go backwards
            return 1.0 if self.state == "done" else self._peak

    def snapshot(self, with_log=True):
        with self._lock:
            now = time.time()
            elapsed = (self.finished or now) - self.started if self.started else 0.0
            frac = self.fraction()
            eta = elapsed * (1 - frac) / frac if self.state == "running" and frac > 0.05 else None
            return {"id": self.id, "kind": self.kind, "state": self.state, "progress": round(frac, 4),
                    "elapsed": round(elapsed, 1), "eta": None if eta is None else round(eta, 1),
                    "stages": [dict(s) for s in self.stages], "error": self.error,
                    "log": list(self.logs[-80:]) if with_log else [],
                    "has_result": self.result is not None,
                    "summary": (self.result or {}).get("summary")}


class ToolHung(RuntimeError):
    """A program stopped making progress (or ran past its time limit) and was killed."""


THREADS = contextvars.ContextVar("engine_threads", default=None)      # per-thread CPU budget set by a parallel stage
_THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                "TF_NUM_INTRAOP_THREADS", "TF_NUM_INTEROP_THREADS", "RAYON_NUM_THREADS")


def _liveness(pgid):
    """A number that changes whenever the process group computes, reads or writes (Linux /proc); None where unavailable."""
    total = 0
    try:
        for name in os.listdir("/proc"):
            if not name.isdigit():
                continue
            try:
                with open(f"/proc/{name}/stat") as f:
                    st = f.read()
                fields = st[st.rindex(")") + 2:].split()
                if int(fields[2]) != pgid:                       # pgrp
                    continue
                total += int(fields[11]) + int(fields[12])      # utime + stime
                try:
                    with open(f"/proc/{name}/io") as f:
                        io = dict(l.split(": ") for l in f.read().splitlines())
                    total += int(io.get("rchar", 0)) + int(io.get("wchar", 0))
                except (OSError, ValueError):
                    pass
            except (OSError, ValueError, IndexError):
                continue
    except OSError:
        return None
    return total


def _limits():
    from . import config
    return config.performance()


def run_cmd(cmd, env=None, job=None, merge=False, on_line=None, cwd=None, stall=None, timeout=None):
    """subprocess.run lookalike. With a job: cancellable. With merge=True stderr is folded into stdout and each line
    is passed to on_line as it arrives. Standard input is always closed (IEDB's bundle otherwise waits for it forever).

    Never hangs: a watchdog kills the program when it shows no CPU, disk or output activity for `stall` seconds
    (Settings > performance.stall_seconds) or runs longer than `timeout`, and raises ToolHung with the reason. Thread
    limits (performance.threads) are put in the environment, and the program runs at a lower priority (performance.nice)."""
    lim = _limits()
    stall = lim["stall_seconds"] if stall is None else stall
    timeout = lim["timeout_seconds"] if timeout is None else timeout
    env = dict(env or os.environ, PYTHONUNBUFFERED="1")     # so progress lines arrive immediately
    threads = THREADS.get() or lim["threads"]
    if threads:
        for v in _THREAD_VARS:
            env[v] = str(threads)
    name = os.path.basename(str(cmd[0]))
    if name in ("bash", "sh") and len(cmd) > 1:
        name = os.path.basename(str(cmd[1]))
    if lim["nice"] and shutil.which("nice") and shutil.which(str(cmd[0])):      # a missing program must still fail as itself
        cmd = ["nice", "-n", str(lim["nice"])] + list(cmd)
    if job:
        job.check()
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT if merge else subprocess.PIPE,
                         stdin=subprocess.DEVNULL, text=True, env=env, bufsize=1, start_new_session=True, cwd=cwd)
    if job:
        job._register(p)
    state = {"why": None, "last_out": time.time()}
    done = threading.Event()

    def watchdog():
        t0 = last_change = time.time()
        sig = _liveness(p.pid)
        tick = max(0.5, min(5.0, (stall or timeout or 20) / 4))
        while not done.wait(tick):
            now = time.time()
            cur = _liveness(p.pid)
            if cur is None or cur != sig or now - state["last_out"] < 5:
                sig, last_change = cur, now
            if timeout and now - t0 > timeout:
                state["why"] = f"ran longer than the {int(timeout)} s limit"
            elif stall and now - last_change > stall and cur is not None:
                state["why"] = f"showed no activity (no CPU, disk or output) for {int(stall)} s"
            if state["why"]:
                for sgn in (signal.SIGTERM, signal.SIGKILL):
                    try:
                        os.killpg(p.pid, sgn)
                    except (ProcessLookupError, PermissionError):
                        break
                    if done.wait(5):
                        return
                return
    wd = threading.Thread(target=watchdog, daemon=True)
    wd.start()
    try:
        if merge:
            out = []
            for line in p.stdout:
                state["last_out"] = time.time()
                out.append(line)
                if on_line:
                    on_line(line.rstrip("\n"))
            p.wait()
            res = subprocess.CompletedProcess(cmd, p.returncode, "".join(out), "")
        else:
            so, se = p.communicate()
            res = subprocess.CompletedProcess(cmd, p.returncode, so, se)
    finally:
        done.set()
        if job:
            job._unregister(p)
    if state["why"]:
        raise ToolHung(f"{name} {state['why']} and was stopped. "
                       "Check the tool in Settings > Test run, or raise performance.stall_seconds if it is only slow.")
    if job:
        job.check()
    return res


class JobManager:
    """One job runs at a time (BigMHC is memory-hungry); the rest wait their turn."""

    def __init__(self, keep=30):
        self._pool = ThreadPoolExecutor(max_workers=1)
        self._jobs, self._keep = {}, keep
        self._lock = threading.Lock()

    def submit(self, fn, kind="pipeline", on_create=None, on_finish=None):
        """on_create(job) runs at once (the run is recorded as queued); on_finish(job) after it ends, whatever the outcome."""
        job = Job(kind)
        if on_create:
            on_create(job)
        with self._lock:
            self._jobs[job.id] = job
            for old in sorted(self._jobs.values(), key=lambda j: j.created)[:-self._keep]:
                if old.state in ("done", "failed", "cancelled"):
                    self._jobs.pop(old.id, None)

        def work():
            if job.cancelled:
                job.state, job.finished, job.started = "cancelled", time.time(), time.time()
                if on_finish:
                    on_finish(job)
                return
            job.state, job.started = "running", time.time()
            try:
                fn(job)
                job.state = "done"
            except Cancelled:
                job.state = "cancelled"
            except Exception as e:                       # surfaced to the page, never lost
                job.state, job.error = "failed", f"{type(e).__name__}: {e}"
                if os.environ.get("ENGINE_DEBUG"):
                    import traceback
                    job.log(traceback.format_exc()[-1500:])
                job.log("ERROR " + job.error)
            finally:
                job.finished = time.time()
                if on_finish:
                    try:
                        on_finish(job)
                    except Exception as e:                   # noqa: BLE001 - recording a run must never hide its result
                        job.log(f"could not record the run: {e}")
        self._pool.submit(work)
        return job

    def get(self, job_id):
        return self._jobs.get(job_id)


MANAGER = JobManager()
