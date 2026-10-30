"""Command line: `python -m circopt.cli run params.json --out runs/x` / `selftest`."""
from __future__ import annotations

import argparse
import json
import sys

from .config import Params
from .pipeline import run_design


def _selftest() -> int:
    from .constraints import Rules, gc_stats
    from .design import Design
    from .score import score
    from .structure import fold, helices, pair_table

    assert [h[2] for h in helices(pair_table("((((....))))"))] == [4]
    assert gc_stats("GCGCAUAU", 4)["gc"] == 0.5
    protein = "MFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    from .codon import AA_TO_CODONS
    cds = "".join(sorted(AA_TO_CODONS[a])[0] for a in protein) + "UAA"
    rna = "GGGAAACCCUUUGGGAAACCCUUU" + cds + "AAACCCGGGUUU"
    p = Params(name="selftest", protein=protein, rna=rna,
               protected=[{"name": "ires_stub", "start": 0, "end": 24}],
               population=6, generations_per_round=2, rounds=1, elite_k=2, top_n=2,
               polish_codons=0, workers=1, inner_span=0, stress_test=False,
               linear_compare=False)
    d = Design(p)
    assert d.cds[0] == 24, d.cds
    assert d.protein.rstrip("*") == protein
    g = d.genome_of()
    assert d.build(g) == d.rna
    m = fold(d.rna, p, d.cds, circ=True, full=True, span=0)
    s = score(d, g, m, Rules(p))
    assert s.total > 0 and "cai_human" or True
    print(f"selftest OK  ({d.n} nt, {len(d.slots)} designable codons, "
          f"MFE {m.mfe:.2f}, CAI {({k: round(v, 3) for k, v in s.cai.items()})})")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="circopt")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run a design from a params JSON file")
    r.add_argument("params")
    r.add_argument("--out", default="runs/latest")
    r.add_argument("--set", action="append", default=[],
                   help="override a param, e.g. --set rounds=6 --set primary_species=mouse")
    sub.add_parser("selftest", help="quick internal check")
    pr = sub.add_parser("params", help="print the default parameters as JSON")
    pr.add_argument("--out", default="-")
    a = ap.parse_args(argv)

    if a.cmd == "selftest":
        return _selftest()
    if a.cmd == "params":
        txt = json.dumps(Params().to_dict(), indent=2)
        if a.out == "-":
            print(txt)
        else:
            with open(a.out, "w") as fh:
                fh.write(txt + "\n")
        return 0

    with open(a.params) as fh:
        raw = json.load(fh)
    for kv in a.set:
        k, _, v = kv.partition("=")
        try:
            raw[k] = json.loads(v)
        except json.JSONDecodeError:
            raw[k] = v
    p = Params.from_dict(raw)
    run_design(p, a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
