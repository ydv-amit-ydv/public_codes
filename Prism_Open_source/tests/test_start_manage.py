"""start.py server management: export / import / reset ×3 / selfpurge (run as a real command)."""
import json
import os
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from helpers import H, env, make_apps, name_ct, register

HERE = Path(__file__).resolve().parents[1]


def sh(args, data, cwd=HERE, answer=""):
    return subprocess.run([sys.executable, str(cwd / "start.py"), *args], capture_output=True, text=True, input=answer,
                          env={**os.environ, "OPENPRISM_DATA": str(data)}, cwd=str(cwd))


def seed(data: Path):
    app, _, store = make_apps(data)
    u = register(app, "alice")
    register(app, "bob")
    assert u.post("/api/projects", json={"name": name_ct(), "data": env()}, headers=H).status_code == 200


def test_export_reset_import_cycle(tmp_path):
    data = tmp_path / "data"
    seed(data)
    out = tmp_path / "backup.zip"
    r = sh(["--export", "--output", str(out)], data)
    assert r.returncode == 0 and "2 account(s) and 1 project(s)" in r.stdout
    with zipfile.ZipFile(out) as z:
        assert {"manifest.json", "openprism.db"} <= set(z.namelist())
        assert json.loads(z.read("manifest.json"))["users"] == 2
    r = sh(["--reset"], data)                                   # one flag is not enough
    assert r.returncode != 0 and (data / "openprism.db").exists()
    r = sh(["--reset", "--reset"], data)
    assert r.returncode != 0 and (data / "openprism.db").exists()
    r = sh(["--reset", "--reset", "--reset"], data)
    assert r.returncode == 0 and not (data / "openprism.db").exists()
    r = sh(["--import", "--input", str(out)], data, answer="no\n")      # confirmation required
    assert r.returncode != 0 and "Cancelled" in r.stderr and not (data / "openprism.db").exists()
    r = sh(["--import", "--input", str(out)], data, answer="IMPORT\n")
    assert r.returncode == 0 and "Restored 2 account(s) and 1 project(s)" in r.stdout
    r = sh(["--import", str(out), "--yes"], data)                 # path right after --import; current data saved first
    assert r.returncode == 0 and "REPLACES everything" in r.stdout and list((data / "backups").glob("before-import-*.db"))


def test_import_refuses_incompatible_exports(tmp_path):
    import sqlite3
    data = tmp_path / "data"
    seed(data)
    out = tmp_path / "backup.zip"
    assert sh(["--export", "--output", str(out)], data).returncode == 0
    newer = tmp_path / "newer.zip"                                 # made by a future version with a newer schema
    with zipfile.ZipFile(out) as z:
        z.extractall(tmp_path / "x")
    con = sqlite3.connect(tmp_path / "x" / "openprism.db")
    con.execute("UPDATE settings SET value = '999' WHERE key = 'schema'"); con.commit(); con.close()
    with zipfile.ZipFile(newer, "w") as z:
        z.write(tmp_path / "x" / "openprism.db", "openprism.db")
        z.write(tmp_path / "x" / "manifest.json", "manifest.json")
    r = sh(["--import", str(newer), "--yes"], data)
    assert r.returncode != 0 and "newer OpenPrism" in r.stderr
    other = tmp_path / "other.zip"                                 # unknown export format
    with zipfile.ZipFile(other, "w") as z:
        z.write(tmp_path / "x" / "openprism.db", "openprism.db")
        z.writestr("manifest.json", json.dumps({"format": "something-else-9"}))
    r = sh(["--import", str(other), "--yes"], data)
    assert r.returncode != 0 and "incompatible" in r.stderr
    assert sh(["--export", "--output", str(tmp_path / "still.zip")], data).returncode == 0   # data untouched


def test_refusals(tmp_path):
    data = tmp_path / "data"
    seed(data)
    (tmp_path / "junk.txt").write_text("hello")
    assert sh(["--import", "--input", str(tmp_path / "junk.txt")], data).returncode != 0
    assert sh(["--export"], data).returncode != 0                 # no output path
    assert "replaced" in sh(["--update"], data).stderr
    assert sh(["--export", "--output", str(tmp_path / "a.zip"), "--port", "9"], data).returncode != 0
    (data / "server.pid").write_text(json.dumps({"pid": os.getpid(), "port": 1}))   # a live process "uses" the data
    r = sh(["--reset", "--reset", "--reset"], data)
    assert r.returncode != 0 and "running" in r.stderr and (data / "openprism.db").exists()
    (data / "server.pid").write_text(json.dumps({"pid": 999999, "port": 1}))        # stale marker is ignored
    assert sh(["--export", "--output", str(tmp_path / "b.zip")], data).returncode == 0


def test_selfpurge(tmp_path):
    app_dir = tmp_path / "Prism_Open_source"
    app_dir.mkdir()
    shutil.copy2(HERE / "start.py", app_dir / "start.py")
    (app_dir / "openprism").mkdir()
    (app_dir / "openprism" / "__init__.py").write_text('__version__ = "test"\n')
    data = app_dir / "projects"
    seed(data)
    assert sh(["--selfpurge"], data, cwd=app_dir).returncode != 0                       # export is required
    assert sh(["--selfpurge", "--export", str(app_dir / "x.zip")], data, cwd=app_dir).returncode != 0  # not inside
    out = tmp_path / "final.zip"
    r = sh(["--selfpurge", "--export", str(out)], data, cwd=app_dir)
    assert r.returncode == 0 and not app_dir.exists() and zipfile.is_zipfile(out)
    # a data folder outside the installation (OPENPRISM_DATA / --data) is removed as well
    app2 = tmp_path / "Prism2"
    app2.mkdir()
    shutil.copy2(HERE / "start.py", app2 / "start.py")
    (app2 / "openprism").mkdir()
    (app2 / "openprism" / "__init__.py").write_text('__version__ = "test"\n')
    ext = tmp_path / "srv-data"
    seed(ext)
    r = sh(["--selfpurge", "--export", str(tmp_path / "final2.zip")], ext, cwd=app2)
    assert r.returncode == 0 and not app2.exists() and not ext.exists() and zipfile.is_zipfile(tmp_path / "final2.zip")
    assert "--import --input" in r.stdout
