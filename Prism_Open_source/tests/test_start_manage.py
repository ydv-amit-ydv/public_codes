"""start.py server commands: --reset ×3, --reset ×3 --loadolddata, --selfpurge ×3 --output (run as real commands)."""
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import zipfile
from pathlib import Path

from helpers import H, env, make_apps, name_ct, register

HERE = Path(__file__).resolve().parents[1]
R3 = ["--reset", "--reset", "--reset"]
P3 = ["--selfpurge", "--selfpurge", "--selfpurge"]


def sh(args, data, cwd=HERE):
    return subprocess.run([sys.executable, str(cwd / "start.py"), *args], capture_output=True, text=True,
                          env={**os.environ, "OPENPRISM_DATA": str(data), "OPENPRISM_NO_START": "1"}, cwd=str(cwd))


def seed(data: Path, users=("alice", "bob")):
    app, _, _ = make_apps(data)
    first = None
    for u in users:
        c = register(app, u)
        first = first or c
    assert first.post("/api/projects", json={"name": name_ct(), "data": env()}, headers=H).status_code == 200


def accounts(data: Path):
    con = sqlite3.connect(data / "openprism.db")
    try:
        return sorted(r[0] for r in con.execute("SELECT username FROM users"))
    finally:
        con.close()


def fake_install(root: Path) -> Path:
    """A minimal copy of the installation (start.py + version + schema) to purge."""
    app = root / "Prism_Open_source"
    (app / "openprism").mkdir(parents=True)
    shutil.copy2(HERE / "start.py", app / "start.py")
    shutil.copy2(HERE / "openprism" / "__init__.py", app / "openprism" / "__init__.py")
    shutil.copy2(HERE / "openprism" / "store.py", app / "openprism" / "store.py")
    return app


def test_reset_needs_three_flags_and_deletes_everything(tmp_path):
    data = tmp_path / "data"
    seed(data)
    for n in (1, 2):
        r = sh(["--reset"] * n, data)
        assert r.returncode != 0 and "three times" in r.stderr and accounts(data) == ["alice", "bob"]
    r = sh(R3, data)
    assert r.returncode == 0 and "deleted" in r.stdout and not (data / "openprism.db").exists()


def test_selfpurge_then_loadolddata_in_a_new_copy(tmp_path):
    old = fake_install(tmp_path / "old")
    data = old / "projects"
    seed(data)
    backups = tmp_path / "backups"
    assert "three times" in sh(["--selfpurge"], data, cwd=old).stderr
    assert sh(P3, data, cwd=old).returncode != 0                                    # --output is required
    assert sh(P3 + ["--output", str(old / "inside")], data, cwd=old).returncode != 0  # not inside the installation
    r = sh(P3 + ["--output", str(backups)], data, cwd=old)
    assert r.returncode == 0 and not old.exists(), r.stderr
    zips = list(backups.glob("openprism-backup-*.zip"))
    assert len(zips) == 1 and "--loadolddata" in r.stdout
    with zipfile.ZipFile(zips[0]) as z:
        assert json.loads(z.read("manifest.json"))["users"] == 2
    # a new copy that already has some other data: --loadolddata replaces it with the backup
    new = fake_install(tmp_path / "new")
    ndata = new / "projects"
    seed(ndata, users=("carol",))
    assert sh(["--loadolddata", str(backups)], ndata, cwd=new).returncode != 0      # needs --reset ×3
    r = sh(R3 + ["--loadolddata", str(backups)], ndata, cwd=new)                     # the folder → newest backup inside
    assert r.returncode == 0 and "OK" in r.stdout, r.stderr
    assert accounts(ndata) == ["alice", "bob"]


def test_bad_or_newer_backup_deletes_nothing(tmp_path):
    data = tmp_path / "data"
    seed(data)
    (tmp_path / "junk.zip").write_text("not a zip")
    r = sh(R3 + ["--loadolddata", str(tmp_path / "junk.zip")], data)
    assert r.returncode != 0 and accounts(data) == ["alice", "bob"]
    # a backup made by a future version (newer data format) is refused before anything is deleted
    other = tmp_path / "other"
    seed(other, users=("zoe",))
    con = sqlite3.connect(other / "openprism.db")
    con.execute("UPDATE settings SET value = '999' WHERE key = 'schema'")
    con.commit()
    con.close()
    newer = tmp_path / "newer.zip"
    with zipfile.ZipFile(newer, "w") as z:
        z.write(other / "openprism.db", "openprism.db")
    r = sh(R3 + ["--loadolddata", str(newer)], data)
    assert r.returncode != 0 and "newer OpenPrism" in r.stderr and accounts(data) == ["alice", "bob"]


def test_old_commands_explain_the_new_ones(tmp_path):
    data = tmp_path / "data"
    seed(data)
    for old in (["--update"], ["--export", "--output", "x.zip"], ["--import", "x.zip"], ["--data", "x"]):
        r = sh(old, data)
        assert r.returncode != 0 and "--loadolddata" in r.stderr
    assert sh(["--output", "x"], data).returncode != 0
    (data / "server.pid").write_text(json.dumps({"pid": os.getpid(), "port": 1}))     # a running server uses the data
    r = sh(R3, data)
    assert r.returncode != 0 and "running" in r.stderr and accounts(data) == ["alice", "bob"]
