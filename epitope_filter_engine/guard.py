"""A model that fails must not throw away the models that worked.

try_model() runs one model; if it raises (the tool is broken, hung or returned nothing) the failure is written to the log, the stage is closed with the
reason, a warning is added to the run, and the pipeline goes on without that model. If every selected model fails, the run fails as before.
"""
from common.jobs import Cancelled, ToolHung


def try_model(name, stage, fn, job, warns, failed):
    """-> fn()'s result, or None when it failed (then `name` is appended to `failed`)."""
    try:
        return fn()
    except Cancelled:
        raise
    except (Exception, ToolHung) as e:                       # noqa: BLE001 - reported to the user, never silently dropped
        msg = f"{type(e).__name__}: {e}"
        job.log(f"\u2717 {name} failed and is skipped: {msg[:400]}")
        job.finish_stage(stage, f"FAILED, skipped: {str(e)[:120]}")
        warns.append(f"{name} failed and was left out of this run ({str(e)[:300]}). The other models were not affected.")
        failed.append(name)
        return None
