"""Helpers shared by the MixMHCpred 3.0 and PRIME 2.1 wrappers (Gfeller lab tools).

Both are shell scripts that call Python/Perl/C++ underneath, take simple allele names
(A0201) and write a tab-separated file with '#' comment lines. Neither fails loudly on
a bad allele (MixMHCpred prints a message, writes nothing and exits 0), so alleles are
checked against each tool's own list first and the output file is checked afterwards.
"""
import os

import pandas as pd

from . import alleles as allele_lib, envclean


class ToolError(RuntimeError):
    pass


def tool_env(python):
    """Clean environment whose first `python3` is the chosen interpreter (the scripts call python3)."""
    return envclean.clean(os.path.dirname(python))


def read_allele_list(path, skip_header=False):
    with open(path) as f:
        lines = [l.split()[0] for l in f if l.strip()]
    return set(lines[1:] if skip_header else lines)


def split_alleles(alleles, supported, tool):
    """-> ({simple_name: canonical}, [notes]). Skips alleles the tool does not know."""
    ok, bad = {}, []
    for a in alleles:
        g = allele_lib.gfeller_name(a)
        if g and g in supported:
            ok[g] = allele_lib.canonical(a)
        else:
            bad.append(a)
    notes = [f"{tool} does not support: " + ", ".join(bad)] if bad else []
    return ok, notes


def chunks(items, n):
    for i in range(0, len(items), n):
        yield items[i:i + n]


def read_table(path):
    return pd.read_csv(path, sep="\t", comment="#")


def require_output(path, result, tool):
    if not os.path.exists(path):
        raise ToolError(f"{tool} wrote no output:\n{(result.stdout or '')[-1500:]}\n{(result.stderr or '')[-800:]}")
