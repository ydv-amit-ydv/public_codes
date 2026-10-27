"""Allele presets and spelling conversion shared by every model.

Internally alleles are canonical ("HLA-A*02:01"). NetMHCpan wants the star
removed ("HLA-A02:01"); BigMHC and IEDB accept the canonical form.
"""
import re

_HLA = re.compile(r"^(?:HLA-)?([ABCEFG])\*?(\d{2,3}):(\d{2,3})", re.I)

_IEDB27 = [
    "A*01:01", "A*02:01", "A*02:03", "A*02:06", "A*03:01", "A*11:01", "A*23:01", "A*24:02",
    "A*26:01", "A*30:01", "A*30:02", "A*31:01", "A*32:01", "A*33:01", "A*68:01", "A*68:02",
    "B*07:02", "B*08:01", "B*15:01", "B*35:01", "B*40:01", "B*44:02", "B*44:03", "B*51:01",
    "B*53:01", "B*57:01", "B*58:01"]
_SUPERTYPES = ["A*01:01", "A*02:01", "A*03:01", "A*24:02", "B*07:02", "B*27:05",
               "B*44:03", "B*58:01", "B*15:01"]

PRESETS = {
    "IEDB human reference 27": {"human": True, "alleles": ["HLA-" + a for a in _IEDB27]},
    "HLA supertype representatives": {"human": True, "alleles": ["HLA-" + a for a in _SUPERTYPES]},
    "Mouse H-2 (all)": {"human": False, "alleles": [
        "H-2-Kb", "H-2-Db", "H-2-Kd", "H-2-Dd", "H-2-Ld", "H-2-Kk", "H-2-Dk"]},
    "Mouse C57BL/6 (Kb + Db)": {"human": False, "alleles": ["H-2-Kb", "H-2-Db"]},
    "Mouse BALB/c (Kd + Dd + Ld)": {"human": False, "alleles": ["H-2-Kd", "H-2-Dd", "H-2-Ld"]},
}


def canonical(a):
    a = a.strip()
    m = _HLA.match(a)
    return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}" if m else a


def netmhcpan_name(a):
    c = canonical(a)
    return c.replace("*", "") if c.startswith("HLA-") else c


def gfeller_name(a):
    """HLA-A*02:01 -> A0201 (MixMHCpred / PRIME naming). None for non-HLA alleles."""
    m = _HLA.match(canonical(a))
    return f"{m.group(1).upper()}{m.group(2)}{m.group(3)}" if m else None
