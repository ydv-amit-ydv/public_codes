"""Find the models without configuration.

Layout the engine expects (names are searched, not fixed):

    ~/local_engines/
        runner_codes/   <- this directory (anything; discovery starts at its parent)
        netmhc/         <- somewhere below it: a `netMHCpan` script beside `Linux_x86_64/`
        bigmhc/         <- somewhere below it: `src/predict.py` beside `models/`

Set ENGINE_MODELS_DIR to search elsewhere. Explicit values in engine_config.json
always win over discovery.
"""
import glob, hashlib, os, re, subprocess, sys

_SKIP = {".git", "data", "models", "nb", "test", "node_modules", "__pycache__", ".venv",
         "lib", "code", "proteome", "pwm", "temp", "input", "output",
         "method", "test_data", "docker", "deploy", "examples", "PWMdef", "nnetsDef", "rpep"}      # inside the tools, never a tool root
_MAX_DEPTH = 6


def models_dir(engine_root):
    return os.path.expanduser(os.environ.get("ENGINE_MODELS_DIR") or os.path.dirname(engine_root))


def _walk(root):
    root = os.path.abspath(root)
    base = root.count(os.sep)
    for dp, dns, fns in os.walk(root):
        if dp.count(os.sep) - base >= _MAX_DEPTH:
            dns[:] = []
        dns[:] = sorted(d for d in dns if d not in _SKIP)
        yield dp, dns, fns


def find_netmhcpan(root):
    """Directory holding the netMHCpan script and its Linux_x86_64 folder (shallowest wins)."""
    hits = [dp for dp, dns, fns in _walk(root)
            if "netMHCpan" in fns and os.path.isdir(os.path.join(dp, "Linux_x86_64"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_bigmhc(root):
    """BigMHC clone: the directory that has src/predict.py and models/."""
    hits = [dp for dp, dns, fns in _walk(root)
            if os.path.isfile(os.path.join(dp, "src", "predict.py"))
            and os.path.isdir(os.path.join(dp, "models"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_mixmhcpred(root):
    """MixMHCpred 3.x clone: a `MixMHCpred` script beside lib/ and code/main.py."""
    hits = [dp for dp, dns, fns in _walk(root)
            if "MixMHCpred" in fns and os.path.isfile(os.path.join(dp, "code", "main.py"))
            and os.path.isdir(os.path.join(dp, "lib"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_prime(root):
    """PRIME 2.x clone: a `PRIME` script beside lib/PRIME.cc."""
    hits = [dp for dp, dns, fns in _walk(root)
            if "PRIME" in fns and os.path.isfile(os.path.join(dp, "lib", "PRIME.cc"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_mixmhc2pred(root):
    """MixMHC2pred clone/release: the Linux binary MixMHC2pred_unix beside nnetsDef/."""
    hits = [dp for dp, dns, fns in _walk(root) if "MixMHC2pred_unix" in fns and os.path.isdir(os.path.join(dp, "nnetsDef"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_immuscope(root):
    """ImmuScope clone with its weights and pseudosequence table: ImmuScope/models/ImmuScope.py, weights/IM/, data/raw/pseudosequence.2023.dat."""
    hits = [dp for dp, dns, fns in _walk_all(root)
            if os.path.isfile(os.path.join(dp, "ImmuScope", "models", "ImmuScope.py"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_tlimmuno2(root):
    """TLimmuno2 clone: Python/TLimmuno2.py beside Python/model/TLimmuno2."""
    hits = [dp for dp, dns, fns in _walk_all(root)
            if os.path.isfile(os.path.join(dp, "Python", "TLimmuno2.py")) and os.path.isdir(os.path.join(dp, "Python", "model", "TLimmuno2"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def _walk_all(root):
    """Like _walk but also enters folders named data/models (they are part of these two repositories); virtualenvs are not entered."""
    root = os.path.abspath(root)
    base = root.count(os.sep)
    skip = _SKIP - {"data", "models"} | {"site-packages", "envs", "env"}
    for dp, dns, fns in os.walk(root):
        if dp.count(os.sep) - base >= 4:
            dns[:] = []
        dns[:] = sorted(d for d in dns if d not in skip and not d.startswith("."))
        yield dp, dns, fns


def find_env_python(home, forced_var, mods, hint):
    """Interpreter for a tool that has its own environment: <VAR> wins, then a venv inside/beside the clone or a conda env named
    like the tool, then the server's own interpreter, that can import `mods`."""
    forced = os.environ.get(forced_var)
    if forced:
        return os.path.expanduser(forced)
    cands = []
    if home:
        for base in (home, os.path.dirname(home)):
            for sub in (".venv", "venv", "env"):
                cands.append(os.path.join(base, sub, "bin", "python"))
    envs = _conda_pythons()
    envs.sort(key=lambda p: (hint not in p.lower(), p))
    for py in cands + envs + [sys.executable]:
        if os.path.isfile(py) and _can_import(py, mods):
            return py
    return None


def find_iedb_bundle_ii(root):
    """IEDB next-generation T Cell Class II standalone: src/tcell_mhcii.py beside pyproject.toml."""
    hits = [dp for dp, dns, fns in _walk(root)
            if "pyproject.toml" in fns and os.path.isfile(os.path.join(dp, "src", "tcell_mhcii.py"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_iedb_bundle_ii_python(root):
    home = find_iedb_bundle_ii(root)
    py = os.path.join(home, ".venv", "bin", "python") if home else None
    return py if py and os.path.isfile(py) else None


def find_iedb_bundle(root):
    """IEDB next-generation T Cell Class I standalone: a folder with src/tcell_mhci.py and pyproject.toml."""
    hits = [dp for dp, dns, fns in _walk(root)
            if "pyproject.toml" in fns and os.path.isfile(os.path.join(dp, "src", "tcell_mhci.py"))]
    return min(hits, key=lambda p: p.count(os.sep)) if hits else None


def find_iedb_bundle_python(root):
    """The bundle's own uv environment (./.venv created by `uv sync`)."""
    home = find_iedb_bundle(root)
    py = os.path.join(home, ".venv", "bin", "python") if home else None
    return py if py and os.path.isfile(py) else None


_py_cache = {}


def _can_import(py, mods):
    key = (py, tuple(mods))
    if key not in _py_cache:
        try:
            r = subprocess.run([py, "-c", "import " + ", ".join(mods)], capture_output=True, timeout=90)
            _py_cache[key] = r.returncode == 0
        except (OSError, subprocess.TimeoutExpired):
            _py_cache[key] = False
    return _py_cache[key]


def _conda_pythons():
    bases = [os.environ.get("CONDA_PREFIX_1"), os.environ.get("CONDA_PREFIX")]
    exe = os.environ.get("CONDA_EXE")
    if exe:
        bases.append(os.path.dirname(os.path.dirname(exe)))
    bases += [os.path.expanduser(p) for p in
              ("~/miniconda3", "~/anaconda3", "~/miniforge3", "~/mambaforge", "~/micromamba", "~/.conda", "/opt/conda")]
    envs = []
    for b in filter(None, bases):
        b = b.split(os.sep + "envs" + os.sep)[0]
        envs += glob.glob(os.path.join(b, "envs", "*", "bin", "python"))
    return list(dict.fromkeys(envs))


def find_python(mods=("numpy", "pandas")):
    """An interpreter that can import `mods`: the server's own first, then conda envs."""
    for py in [sys.executable] + _conda_pythons():
        if os.path.isfile(py) and _can_import(py, mods):
            return py
    return None


def find_mhcflurry_python():
    """An interpreter that can import mhcflurry. MHCFLURRY_PYTHON wins; an env named like
    'mhcflurry' is tried first."""
    forced = os.environ.get("MHCFLURRY_PYTHON")
    if forced:
        return os.path.expanduser(forced)
    envs = _conda_pythons()
    envs.sort(key=lambda p: ("mhcflurry" not in p.lower(), p))
    for py in [sys.executable] + envs:
        if os.path.isfile(py) and _can_import(py, ("mhcflurry",)):
            return py
    return None


def find_bigmhc_python():
    """An interpreter that can import torch. BIGMHC_PYTHON wins, then conda envs
    (an env named like 'bigmhc' first), then the interpreter running the server."""
    forced = os.environ.get("BIGMHC_PYTHON")
    if forced:
        return os.path.expanduser(forced)
    envs = _conda_pythons()
    envs.sort(key=lambda p: ("bigmhc" not in p.lower(), p))
    for py in envs + [sys.executable]:
        if os.path.isfile(py) and _can_import(py, ("torch", "numpy", "pandas", "psutil")):
            return py
    return None


_NMHOME = re.compile(r"^(setenv\s+NMHOME\s+)(\S+)", re.M)


def runnable_netmhcpan(home, workdir):
    """Path of a netMHCpan script whose NMHOME points at `home`.

    DTU's script ships with a placeholder home that has to be edited. Instead of
    modifying your copy, run a patched copy from workdir when (and only when) needed.
    """
    script = os.path.join(home, "netMHCpan")
    try:
        text = open(script).read()
    except OSError:
        return script
    m = _NMHOME.search(text)
    if not m or m.group(2).rstrip("/") == home.rstrip("/"):
        return script
    patched = _NMHOME.sub(lambda mm: mm.group(1) + home, text, count=1)
    tag = hashlib.md5(home.encode()).hexdigest()[:8]
    out = os.path.join(workdir, f"netMHCpan.{tag}")
    if not os.path.exists(out) or open(out).read() != patched:
        with open(out, "w") as f:
            f.write(patched)
    os.chmod(out, 0o755)
    return out
