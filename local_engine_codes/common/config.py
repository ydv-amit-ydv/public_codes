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
    "performance": {"threads": 0, "nice": 5, "stall_seconds": 900, "timeout_minutes": 0, "parallel_models": 2},
    "population": {"label": "India", "weights_csv": "auto"},
    "models": {
        "netmhcpan": {"enabled": True, "home": "auto", "workers": 0, "alleles_per_call": 3},
        "bigmhc": {"enabled": True, "python": "auto", "home": "auto", "device": "cpu", "jobs": 0},
        "mhcflurry": {"enabled": True, "python": "auto", "command": "auto"},
        "iedb_bundle": {"enabled": True, "home": "auto", "python": "auto"},
        "mixmhc2pred": {"enabled": True, "home": "auto"},
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
    if not isinstance(perf, dict) or any(isinstance(v, bool) or not isinstance(v, (int, float)) for v in perf.values()):
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


def performance():
    """Resource limits applied to every program the engine starts. threads 0 = let each tool use all cores; nice 0-19 lowers
    the priority so the machine stays responsive; stall_seconds 0 = never kill a quiet program; timeout_minutes 0 = no limit;
    parallel_models = how many models of one funnel stage run at the same time (they share the thread budget)."""
    perf = {**DEFAULT["performance"], **(load().get("performance") or {})}
    cpus = os.cpu_count() or 1
    return {"threads": max(0, int(perf["threads"])), "nice": min(19, max(0, int(perf["nice"]))), "cpus": cpus,
            "stall_seconds": max(0, int(perf["stall_seconds"])), "timeout_seconds": max(0, int(perf["timeout_minutes"])) * 60,
            "parallel_models": max(1, min(8, int(perf["parallel_models"])))}


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
