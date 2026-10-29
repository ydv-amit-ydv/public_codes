"""Engine settings: one JSON file, re-read on every call so edits made in the GUI
(or by hand) take effect without restarting the server.

Each model runs as a subprocess with its own interpreter / install directory, so
the Flask server never imports torch or any model code.
"""
import copy, json, os, tempfile

from . import discover

ENGINE_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_PATH = os.path.expanduser(
    os.environ.get("ENGINE_CONFIG", os.path.join(ENGINE_ROOT, "engine_config.json")))

# "auto" = look beside this directory (see common/discover.py). Explicit paths win.
DEFAULT = {
    "tmp_dir": "auto",
    "workspace_dir": "auto",
    "performance": {"profile": "auto", "threads": 0, "nice": 5, "stall_seconds": 900, "timeout_minutes": 0, "parallel_models": 2, "cache": 1},
    "population": {"label": "India", "weights_csv": "auto"},
    "models": {
        "netmhcpan": {"enabled": True, "home": "auto", "workers": 0, "alleles_per_call": 3},
        "bigmhc": {"enabled": True, "python": "auto", "home": "auto", "device": "cpu", "jobs": 0},
        "mhcflurry": {"enabled": True, "python": "auto", "command": "auto"},
        "iedb_bundle": {"enabled": True, "home": "auto", "python": "auto"},
        "mixmhc2pred": {"enabled": True, "home": "auto"},
        "immuscope": {"enabled": True, "python": "auto", "home": "auto"},
        "tlimmuno2": {"enabled": True, "python": "auto", "home": "auto"},
        "iedb_bundle_ii": {"enabled": True, "home": "auto", "python": "auto"},
        "mixmhcpred": {"enabled": True, "python": "auto", "home": "auto"},
        "prime": {"enabled": True, "python": "auto", "home": "auto"},
    },
}

_FINDERS = {("netmhcpan", "home"): lambda root: discover.find_netmhcpan(root),
            ("bigmhc", "home"): lambda root: discover.find_bigmhc(root),
            ("bigmhc", "python"): lambda root: discover.find_bigmhc_python(),
            ("mhcflurry", "python"): lambda root: discover.find_mhcflurry_python(),
            ("mixmhc2pred", "home"): lambda root: discover.find_mixmhc2pred(root),
            ("immuscope", "home"): lambda root: discover.find_immuscope(root),
            ("immuscope", "python"): lambda root: discover.find_env_python(discover.find_immuscope(root), "IMMUSCOPE_PYTHON", ("torch", "numpy", "pandas", "h5py", "sklearn", "ruamel.yaml"), "immuscope"),
            ("tlimmuno2", "home"): lambda root: discover.find_tlimmuno2(root),
            ("tlimmuno2", "python"): lambda root: discover.find_env_python(discover.find_tlimmuno2(root), "TLIMMUNO2_PYTHON", ("tensorflow", "numpy", "pandas"), "tlimmuno2"),
            ("iedb_bundle_ii", "home"): lambda root: discover.find_iedb_bundle_ii(root),
            ("iedb_bundle_ii", "python"): lambda root: discover.find_iedb_bundle_ii_python(root),
            ("iedb_bundle", "home"): lambda root: discover.find_iedb_bundle(root),
            ("iedb_bundle", "python"): lambda root: discover.find_iedb_bundle_python(root),
            ("mixmhcpred", "home"): lambda root: discover.find_mixmhcpred(root),
            ("mixmhcpred", "python"): lambda root: discover.find_python(("numpy", "pandas")),
            ("prime", "home"): lambda root: discover.find_prime(root),
            ("prime", "python"): lambda root: discover.find_python(("numpy", "pandas"))}

# key -> accepted type, per known model key
_TYPES = {"enabled": bool, "home": str, "python": str, "device": str, "jobs": int,
          "workers": int, "alleles_per_call": int, "command": str}


class ConfigError(ValueError):
    pass


def _merge(base, over):
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _merge(base[k], v)
        else:
            base[k] = v
    return base


def validate(obj):
    """Return a list of problems (empty if the config is usable)."""
    errs = []
    if not isinstance(obj, dict):
        return ["top level must be a JSON object"]
    for k in ("tmp_dir", "workspace_dir"):
        if k in obj and not isinstance(obj[k], str):
            errs.append(f"{k} must be a string")
    perf = obj.get("performance", {})
    if not isinstance(perf, dict) or any(k != "profile" and (isinstance(v, bool) or not isinstance(v, (int, float))) for k, v in perf.items()) \
            or perf.get("profile", "auto") not in ("auto", "gentle", "balanced", "max"):
        errs.append("performance must be an object of numbers (threads, nice, stall_seconds, timeout_minutes, parallel_models)")
    models = obj.get("models", {})
    if not isinstance(models, dict):
        return errs + ["'models' must be an object"]
    for name, m in models.items():
        if not isinstance(m, dict):
            errs.append(f"models.{name} must be an object")
            continue
        for k, v in m.items():
            want = _TYPES.get(k)
            if want and (not isinstance(v, want) or (want is int and isinstance(v, bool))):
                errs.append(f"models.{name}.{k} must be {want.__name__}")
    return errs


def load():
    cfg = copy.deepcopy(DEFAULT)
    if os.path.exists(CONFIG_PATH):
        try:
            with open(CONFIG_PATH) as f:
                user = json.load(f)
        except json.JSONDecodeError as e:
            raise ConfigError(f"{CONFIG_PATH} is not valid JSON: {e}")
        errs = validate(user)
        if errs:
            raise ConfigError(f"{CONFIG_PATH}: " + "; ".join(errs))
        _merge(cfg, user)
    return cfg


def search_root():
    return discover.models_dir(ENGINE_ROOT)


def model(name):
    """Settings for one model with '~' expanded and 'auto'/empty paths resolved by
    discovery. `auto` lists the keys that were discovered. Unknown models give
    {'enabled': False}."""
    m = dict(load()["models"].get(name, {"enabled": False}))
    auto = []
    for k in ("home", "python", "command"):
        v = m.get(k)
        finder = _FINDERS.get((name, k))
        if finder and (not v or v == "auto"):
            m[k] = finder(search_root()) or ""
            auto.append(k)
        elif v:
            m[k] = os.path.expanduser(v)
    m["auto"] = auto
    return m


def tmp_dir():
    v = load()["tmp_dir"]
    d = os.path.join(tempfile.gettempdir(), "epitope_engine") if v in ("", "auto") else os.path.expanduser(v)
    os.makedirs(d, exist_ok=True)
    return d


PROFILES = {"gentle": 0.25, "balanced": 0.5, "max": 1.0}
_ACTIVE = {"profile": None}


def set_profile(name):
    """The profile of the job that is running now (one job runs at a time); None = the one in Settings."""
    _ACTIVE["profile"] = name if name in PROFILES else None


def machine():
    """CPU threads this process may use and the memory that is free (cgroup limits respected)."""
    try:
        cpus = len(os.sched_getaffinity(0))
    except AttributeError:
        cpus = os.cpu_count() or 1
    mem = 0.0
    try:
        with open("/proc/meminfo") as f:
            info = {l.split(":")[0]: int(l.split()[1]) for l in f if ":" in l}
        mem = info.get("MemAvailable", info.get("MemFree", 0)) / 1048576
        for path in ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
            try:
                v = open(path).read().strip()
                if v.isdigit() and int(v) < 1 << 50:
                    mem = min(mem, int(v) / 2 ** 30)
            except OSError:
                pass
    except (OSError, ValueError):
        pass
    return cpus, round(mem, 1)


def performance():
    """Resource budget applied to every program the engine starts, from the profile (gentle 1/4 of the cores, balanced 1/2,
    max all but one on a big machine) and the limits in Settings:
      workers          single-threaded programs run side by side (NetMHCpan, PRIME, ...)
      threads          threads each multi-threaded program may use (BigMHC, MHCflurry); Settings.threads > 0 overrides
      shards           BigMHC / MHCflurry processes run side by side on slices of the work (memory permitting)
      parallel_models  models of one funnel stage that run at the same time
      nice             operating-system priority; stall_seconds / timeout_seconds: watchdog limits (0 = off)."""
    perf = {**DEFAULT["performance"], **(load().get("performance") or {})}
    cpus, mem = machine()
    profile = _ACTIVE["profile"] or perf.get("profile", "auto")
    if profile not in PROFILES:                     # "auto": a big, idle machine runs flat out; a small one leaves room for the user
        profile = "max" if (cpus >= 16 and mem >= 16) else "balanced"
    cores = max(1, cpus - 1 if (profile == "max" and cpus > 8) else round(cpus * PROFILES[profile]))
    mem_cap = max(1, int(mem // 3)) if mem else cores
    pm = max(1, min(8, int(perf["parallel_models"])))
    if profile == "gentle":
        pm = 1
    elif profile == "max":
        pm = max(pm, min(4, cores // 6))
    shards = {"gentle": 1, "balanced": max(1, min(4, cores // 8)), "max": max(1, min(8, cores // 4))}[profile]
    return {"profile": profile, "cpus": cpus, "mem_gb": mem, "cores": cores, "workers": max(1, min(cores, mem_cap * 4)),
            "threads": max(0, int(perf["threads"])) or cores, "shards": max(1, min(shards, mem_cap)),
            "nice": {"gentle": 10, "max": 0}.get(profile, min(19, max(0, int(perf["nice"])))),
            "stall_seconds": max(0, int(perf["stall_seconds"])), "timeout_seconds": max(0, int(perf["timeout_minutes"])) * 60,
            "parallel_models": pm, "cache": 1 if perf.get("cache", 1) else 0}


def workspace():
    """Where every run is kept (inputs, settings, tables, reports): ENGINE_WORKSPACE, else workspace_dir, else
    <models folder>/engine_workspace, so it survives replacing the code folder."""
    v = os.environ.get("ENGINE_WORKSPACE") or load().get("workspace_dir", "auto")
    d = os.path.join(search_root(), "engine_workspace") if v in ("", "auto") else os.path.expanduser(v)
    os.makedirs(os.path.join(d, "runs"), exist_ok=True)
    return d


def run_dir(run_id):
    """Folder of one run; `run_id` may be a job id or a view id like '<job>_c1'."""
    d = os.path.join(workspace(), "runs", str(run_id).split("_")[0])
    os.makedirs(d, exist_ok=True)
    return d


def effective():
    """What the server will actually use right now, for display in Settings."""
    out = {"search_root": search_root(), "tmp_dir": tmp_dir(), "workspace": workspace(), "models": {}}
    for name in load()["models"]:
        m = model(name)
        out["models"][name] = {k: (v if v != "" else None) for k, v in m.items()}
    return out


def raw_text():
    """What the GUI editor shows: the user's file if present, else the defaults."""
    if os.path.exists(CONFIG_PATH):
        with open(CONFIG_PATH) as f:
            return f.read()
    return json.dumps(DEFAULT, indent=2) + "\n"


def save_text(text):
    """Validate then atomically write. Raises ConfigError and leaves the file alone on a bad edit."""
    try:
        obj = json.loads(text)
    except json.JSONDecodeError as e:
        raise ConfigError(f"invalid JSON: {e}")
    errs = validate(obj)
    if errs:
        raise ConfigError("; ".join(errs))
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(CONFIG_PATH), suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        f.write(json.dumps(obj, indent=2) + "\n")
    os.replace(tmp, CONFIG_PATH)
    return obj
