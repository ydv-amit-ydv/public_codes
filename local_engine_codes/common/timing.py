"""How long each tool takes, learned from the runs on this machine, so the progress bar and time left follow the real work.

A stage's weight in the overall bar is its expected seconds = predictions x seconds per prediction. Defaults are rough; after every run
the measured rate of each tool is stored (average of the last five), so the second run's bar is already close.
"""
import json
import os

from . import config

# seconds per prediction (or per peptide where a tool scores peptides), a fixed start-up cost, on a small machine
DEFAULT = {"netmhcpan": (3.4e-3, 5), "mixmhcpred": (2e-5, 2), "prime": (3e-5, 3), "iedb_imm": (1e-5, 1), "mhcflurry": (1.5e-3, 15),
           "netctlpan": (3e-4, 10), "netchop": (1e-4, 15), "bigmhc_el": (2e-2, 30), "bigmhc_im": (2e-2, 30),
           "netmhciipan": (2e-3, 10), "mixmhc2pred": (2e-5, 2), "cd4episcore": (3e-3, 10), "merge": (2e-6, 5), "regions": (0, 10)}
GATE = 0.03        # share of the pairs BigMHC gets after the funnel (a guess until a run teaches it)


def _path():
    return os.path.join(config.workspace(), "timings.json")


def _load():
    try:
        with open(_path()) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def expected_seconds(model, work):
    rate, fixed = DEFAULT.get(model, (1e-3, 5))
    learned = _load().get(model)
    if learned:
        rate = sum(learned) / len(learned)
    return fixed + rate * (work or 0)


def learn(stages):
    """stages: the finished job's stage dicts (names may carry a class prefix like 'c1_netmhcpan')."""
    data = _load()
    for s in stages:
        if s.get("state") != "done" or not s.get("work") or not s.get("t0") or not s.get("t1"):
            continue
        model = s["name"].split("_", 1)[1] if s["name"][:2] in ("c1", "c2") else s["name"]
        secs = s["t1"] - s["t0"]
        if model not in DEFAULT or secs < 1:
            continue
        fixed = DEFAULT[model][1]
        if (secs - fixed) / s["work"] < 0.05 * DEFAULT[model][0]:
            continue                                     # served from the prediction cache: says nothing about the tool's speed
        data[model] = (data.get(model, []) + [max(0.0, secs - fixed) / s["work"]])[-5:]
    try:
        with open(_path(), "w") as f:
            json.dump(data, f)
    except OSError:
        pass
