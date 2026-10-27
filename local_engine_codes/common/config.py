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
    "models": {
        "netmhcpan": {"enabled": True, "home": "auto", "workers": 0, "alleles_per_call": 3},
        "bigmhc": {"enabled": True, "python": "auto", "home": "auto", "device": "cpu", "jobs": 0},
        "mhcflurry": {"enabled": True, "python": "auto", "command": "auto"},
        "iedb_bundle": {"enabled": True, "home": "auto", "python": "auto"},
        "mixmhcpred": {"enabled": True, "python": "auto", "home": "auto"},
        "prime": {"enabled": True, "python": "auto", "home": "auto"},
    },
}

_FINDERS = {("netmhcpan", "home"): lambda root: discover.find_netmhcpan(root),
            ("bigmhc", "home"): lambda root: discover.find_bigmhc(root),
            ("bigmhc", "python"): lambda root: discover.find_bigmhc_python(),
            ("mhcflurry", "python"): lambda root: discover.find_mhcflurry_python(),
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
    if "tmp_dir" in obj and not isinstance(obj["tmp_dir"], str):
        errs.append("tmp_dir must be a string")
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


def effective():
    """What the server will actually use right now, for display in Settings."""
    out = {"search_root": search_root(), "tmp_dir": tmp_dir(), "models": {}}
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
