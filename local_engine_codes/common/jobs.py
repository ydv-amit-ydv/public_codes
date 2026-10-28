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
        self._tasks, self._tid = {}, 0          # every program started for this job: what it is working on, since when, how it ended

    # --- progress ---------------------------------------------------------------
    def stage(self, name, label, weight=1.0, total=None, work=None):
        """`work` = how many predictions the stage will make (kept to learn each tool's speed for the next run's progress bar)."""
        with self._lock:
            entry = {"name": name, "label": label, "weight": weight, "done": 0,
                     "total": total, "state": "pending", "msg": "", "work": work, "t0": None, "t1": None}
            for i, s in enumerate(self.stages):
                if s["name"] == name:
                    self.stages[i] = entry                    # registered again later (by the pipeline itself): keep its place
                    return
            self.stages.append(entry)

    def update(self, name, done=None, total=None, msg=None, state=None):
        with self._lock:
            for s in self.stages:
                if s["name"] == name:
                    if done is not None: s["done"] = done
                    if total is not None: s["total"] = total
                    if msg is not None: s["msg"] = msg
                    if state: s["state"] = state
                    elif s["state"] == "pending": s["state"] = "running"
                    if s["state"] == "running" and not s.get("t0"):
                        s["t0"] = time.time()

    def finish_stage(self, name, msg=None):
        with self._lock:
            for s in self.stages:
                if s["name"] == name:
                    s["state"] = "done"
                    s["t1"] = time.time()
                    s["t0"] = s.get("t0") or s["t1"]
                    if s["total"]: s["done"] = s["total"]
                    if msg is not None: s["msg"] = msg

    def log(self, line):
        with self._lock:
            self.logs.append(line)
            del self.logs[:-400]

    def check(self):
        if self.cancelled:
            raise Cancelled()

    # --- tasks: which batch is on which engine ------------------------------------------
    def task_begin(self, model, what, n=None, sample=None, pid=None):
        with self._lock:
            self._tid += 1
            now = time.time()
            self._tasks[self._tid] = {"id": self._tid, "model": model, "what": what, "n": n, "sample": list(sample or [])[:3], "start": now,
                                      "last_active": now, "cpu_s": 0.0, "status": "running", "note": "", "pid": pid, "end": None}
            tid = self._tid
        self.log(f"\u25b6 {model} \u00b7 {what}" + (f" \u2014 {n:,} peptides" if n else "") + (f" (e.g. {', '.join(list(sample)[:3])})" if sample else ""))
        return tid

    def task_touch(self, tid, cpu_s, active):
        with self._lock:
            t = self._tasks.get(tid)
            if t and t["status"] == "running":
                t["cpu_s"] = round(cpu_s, 1)
                if active:
                    t["last_active"] = time.time()

    def task_end(self, tid, status, note=""):
        with self._lock:
            t = self._tasks.get(tid)
            if not t:
                return
            t["status"], t["note"], t["end"] = status, note, time.time()
            dur = t["end"] - t["start"]
        mark = "\u2713" if status == "done" else "\u2717"
        self.log(f"{mark} {t['model']} \u00b7 {t['what']} \u2014 {status} in {dur:.0f} s" + (f": {note}" if note and status != "done" else ""))

    def tasks_snapshot(self):
        if not self._tasks and getattr(self, "_saved_tasks", None):
            return self._saved_tasks
        with self._lock:
            now = time.time()
            rows = [dict(t, elapsed=round((t["end"] or now) - t["start"], 1), quiet=round(now - t["last_active"], 1) if t["status"] == "running" else 0)
                    for t in self._tasks.values()]
        by = {}
        for t in rows:
            m = by.setdefault(t["model"], {"running": 0, "done": 0, "failed": 0, "seconds": 0.0, "slowest": 0.0, "peptides": 0})
            if t["status"] == "running":
                m["running"] += 1
            elif t["status"] == "done":
                m["done"] += 1
                m["seconds"] += t["elapsed"]
                m["slowest"] = max(m["slowest"], t["elapsed"])
                m["peptides"] += t["n"] or 0
            else:
                m["failed"] += 1
        finished = sorted((t for t in rows if t["status"] != "running"), key=lambda t: -(t["end"] or 0))
        return {"running": [t for t in rows if t["status"] == "running"], "recent": finished[:40],
                "problems": [t for t in finished if t["status"] != "done"][:20], "by_model": by}

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
                    "log": list(self.logs[-150:]) if with_log else [],
                    "has_result": self.result is not None, "tasks": self.tasks_snapshot() if (self._tasks or getattr(self, "_saved_tasks", None)) else None,
                    "summary": (self.result or {}).get("summary")}


class ToolHung(RuntimeError):
    """A program stopped making progress (or ran past its time limit) and was killed."""


THREADS = contextvars.ContextVar("engine_threads", default=None)      # per-thread CPU budget set by a parallel stage
_THREAD_VARS = ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                "TF_NUM_INTRAOP_THREADS", "TF_NUM_INTEROP_THREADS", "RAYON_NUM_THREADS")


_TICKS = os.sysconf("SC_CLK_TCK") if hasattr(os, "sysconf") else 100


def _liveness(pgid):
    """(activity number, CPU seconds) of the process group: the number changes whenever it computes, reads or writes (Linux /proc);
    None where unavailable."""
    total, cpu = 0, 0
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
                cpu += int(fields[11]) + int(fields[12])
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
    return total, cpu / _TICKS


def _limits():
    from . import config
    return config.performance()


def run_cmd(cmd, env=None, job=None, merge=False, on_line=None, cwd=None, stall=None, timeout=None, task=None):
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
    tid = job.task_begin(task["model"], task["what"], task.get("n"), task.get("sample"), pid=p.pid) if (job and task and hasattr(job, "task_begin")) else None
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
        if tid is not None:
            tail = ((locals().get("res") and (res.stderr or res.stdout)) or "").strip().splitlines()[-1:] or [""]
            rc = p.returncode
            status = "stalled" if state["why"] and "activity" in state["why"] else "timeout" if state["why"] else "cancelled" if (job and job.cancelled) else "done" if rc == 0 else "failed"
            job.task_end(tid, status, state["why"] or (f"exit {rc}: {tail[0][:160]}" if rc else ""))
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
