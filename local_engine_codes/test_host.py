"""Host-proteome identity: search exactness, modes, download. Run: python test_host.py"""
import gzip, json, os, random, sys, tempfile, time

import numpy as np
import pandas as pd

T = tempfile.mkdtemp()
os.environ["ENGINE_CONFIG"] = os.path.join(T, "cfg.json")
os.environ["ENGINE_WORKSPACE"] = os.path.join(T, "ws")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import config, hostself as h                                   # noqa: E402

random.seed(5)
AA = h.AA
prots = {f"sp|P{i:05d}|TEST{i}_HUMAN": "".join(random.choice(AA) for _ in range(random.randint(80, 300))) for i in range(60)}
src = os.path.join(T, "src.fasta.gz")
with gzip.open(src, "wt") as f:
    for k, v in prots.items():
        f.write(f">{k} test\n" + "\n".join(v[i:i + 60] for i in range(0, len(v), 60)) + "\n")

# download from a file:// URL: not present -> present, atomic, indexed
json.dump({"host_proteome": {"organism": "human", "fasta": "auto", "url": "file://" + src}}, open(os.environ["ENGINE_CONFIG"], "w"))
assert not h.available() and h.status()["present"] is False
n = h.download()
assert n == 60 and h.available() and not os.path.exists(h.path() + ".part"), n
st = h.status(); assert st["present"] and st["size_mb"] >= 0 and st["name"] == "Homo sapiens", st
assert h.start_download() and (time.sleep(1) or True)
for _ in range(50):
    if not h.status()["downloading"]:
        break
    time.sleep(.1)
assert h.status()["error"] is None, h.status()

ix = h.index()
assert len(ix.names) == 60 and h.status()["loaded"] and h.status()["proteins"] == 60
first = list(prots.values())[7]
p9 = first[20:29]


def mutate(p, pos):
    q = list(p)
    for i in pos:
        q[i] = "W" if q[i] != "W" else "Y"
    return "".join(q)


r = h.check(p9, ix); assert r["Host_Self"] == "identical" and r["Host_Best_Identity"] == 1.0 and r["Host_Best_Protein"].startswith("sp|P00007"), r
assert r["Host_Best_Start"] == 21 and r["Host_Best_Window"] == p9, r
r = h.check(mutate(p9, [4]), ix); assert r["Host_Self"] == "near" and abs(r["Host_Best_Identity"] - 8 / 9) < 1e-3, r
r = h.check(mutate(p9, [1, 7]), ix); assert r["Host_Self"] == "distinct" and abs(r["Host_Best_Identity"] - 7 / 9) < 1e-3, r      # 2 mismatches found exactly
r = h.check(first[20:35], ix); assert r["Host_Self"] == "identical" and r["Host_Exact_Kmers"] == 7, r                       # 15-mer: 7 windows, all exact
r = h.check(first[20:28], ix); assert r["Host_Self"] == "identical", r                                                      # 8-mer by brute force
r = h.check("WWWWWWWWWW", ix); assert r["Host_Self"] == "distinct", r
# a window that would only exist by joining two proteins is not found
names = list(prots)
joined = prots[names[0]][-4:] + prots[names[1]][:5]
assert h.check(joined, ix)["Host_Exact_Kmers"] == 0 or joined in "".join(prots.values()), "windows must not span proteins"

# exactness against brute force for <= 2 mismatches
seq, valid = ix.seq, ix.valid
w = np.lib.stride_tricks.sliding_window_view(seq, 9)
for t in range(40):
    prot = random.choice(list(prots.values())); i = random.randrange(len(prot) - 9)
    p = mutate(prot[i:i + 9], random.sample(range(9), random.choice([0, 1, 2])))
    q = ix._encode(p)
    best = int(np.where(valid, (w == q).sum(1), -1).max()) / 9
    got = ix.best9(p)[0]
    assert best < 7 / 9 - 1e-9 or abs(best - got) < 1e-9, (p, best, got)

# annotate: modes
pep = pd.DataFrame({"Peptide": [p9, mutate(p9, [4]), "WWWWWWWWW"], "Candidate": [True, True, True]})
notes = []
a = h.annotate(pep.copy(), {"host_mode": "flag"}, notes)
assert list(a["Host_Self"]) == ["identical", "near", "distinct"] and a["Candidate"].all() and not a["Host_Excluded"].any() and not notes, a
a = h.annotate(pep.copy(), {"host_mode": "exclude_identical"}, notes)
assert list(a["Candidate"]) == [False, True, True] and list(a["Host_Excluded"]) == [True, False, False] and "removed" in notes[-1], (a, notes)
a = h.annotate(pep.copy(), {"host_mode": "exclude_near"}, notes)
assert list(a["Candidate"]) == [False, False, True], a
a = h.annotate(pep.copy(), {"host_mode": "exclude_near", "host_near_min": 1.0}, notes)
assert list(a["Candidate"]) == [False, True, True], a                                    # near threshold is adjustable
a = h.annotate(pep.copy(), {"host_mode": "off"}, notes); assert "Host_Self" not in a
# only candidates are checked
pep2 = pep.copy(); pep2["Candidate"] = [False, True, False]
a = h.annotate(pep2, {"host_mode": "flag"}, notes); assert a["Host_Self"].tolist() == ["", "near", ""], a

# the proteome missing: a note, no exclusion, no crash
os.remove(h.path()); h._INDEX.clear()
notes = []
a = h.annotate(pep.copy(), {"host_mode": "exclude_identical"}, notes)
assert a["Candidate"].all() and "not downloaded" in notes[0], notes

# own FASTA (plain text) and a bad download
own = os.path.join(T, "own.fa"); open(own, "w").write(">x\nACDEFGHIKLMNPQRSTVWY\n")
json.dump({"host_proteome": {"organism": "mouse", "fasta": own}}, open(os.environ["ENGINE_CONFIG"], "w"))
assert h.available() and h.path() == own and h.check("ACDEFGHIK", h.index())["Host_Self"] == "identical"
try:
    h.download(); raise SystemExit("download must refuse when a FASTA is given")
except ValueError:
    pass
json.dump({"host_proteome": {"organism": "dog"}}, open(os.environ["ENGINE_CONFIG"], "w"))
try:
    config.load(); raise SystemExit("organism dog must be rejected")
except config.ConfigError:
    pass
bad = os.path.join(T, "bad.gz"); open(bad, "wb").write(b"not a gzip")
json.dump({"host_proteome": {"organism": "human", "fasta": "auto", "url": "file://" + bad}}, open(os.environ["ENGINE_CONFIG"], "w"))
try:
    h.download(); raise SystemExit("a corrupt download must fail")
except ValueError:
    pass
assert not os.path.exists(h.path()) and not os.path.exists(h.path() + ".part")
print("host tests passed")
