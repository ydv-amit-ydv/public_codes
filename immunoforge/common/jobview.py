"""Lets two pipelines share one Job: stage names get a prefix ("c1_netmhcpan") and labels a tag ("MHC I - NetMHCpan").

The pipelines only call stage / update / finish_stage / log / check / process registration and set `.result`; this view
forwards all of those to the real job, so cancelling, progress and the log work for the combined run unchanged.
"""


class JobView:
    def __init__(self, job, prefix, tag, share=1.0):
        self._job, self._p, self._tag, self._share = job, prefix, tag, share
        self.result = None

    @property
    def id(self):
        return f"{self._job.id}_{self._p}"

    @property
    def cancelled(self):
        return self._job.cancelled

    def stage(self, name, label, weight=1.0, total=None, work=None):
        self._job.stage(f"{self._p}_{name}", f"{self._tag}: {label}", weight * self._share, total, work)

    def update(self, name, *a, **k):
        self._job.update(f"{self._p}_{name}", *a, **k)

    def finish_stage(self, name, msg=None):
        self._job.finish_stage(f"{self._p}_{name}", msg)

    def log(self, line):
        self._job.log(f"[{self._tag}] {line}")

    def check(self):
        self._job.check()

    def _register(self, p):
        self._job._register(p)

    def _unregister(self, p):
        self._job._unregister(p)

    def task_begin(self, model, what, n=None, sample=None, pid=None):
        return self._job.task_begin(f"{self._tag} \u00b7 {model}", what, n, sample, pid)

    def task_touch(self, tid, cpu_s, active):
        self._job.task_touch(tid, cpu_s, active)

    def task_end(self, tid, status, note=""):
        self._job.task_end(tid, status, note)
