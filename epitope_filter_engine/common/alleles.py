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


# ---------------------------------------------------------------- MHC class II ----------------
# Canonical class II names: HLA-DRB1*01:01 (DR, beta chain only), HLA-DQA1*05:01/DQB1*02:01 and
# HLA-DPA1*01:03/DPB1*04:01 (alpha/beta pairs), H2-IAb (mouse). This is the spelling IEDB's tools and
# MixMHC2pred's allele tables ("OfficialName") use.
_CHAIN = r"(DRB[1345]|DQA1|DQB1|DPA1|DPB1)[\*_](\d{2,3})[:_](\d{2,3})"
_II_SINGLE = re.compile(rf"^(?:HLA-)?{_CHAIN}$", re.I)
_II_PAIR = re.compile(rf"^(?:HLA-)?{_CHAIN}(?:/|__)(DQB1|DPB1)[\*_](\d{{2,3}})[:_](\d{{2,3}})$", re.I)
_MOUSE_II = re.compile(r"^H-?2-?(I[AE][a-z0-9]*)$", re.I)


def canonical_ii(a):
    a = a.strip()
    m = _II_PAIR.match(a)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}/{m.group(4).upper()}*{m.group(5)}:{m.group(6)}"
    m = _II_SINGLE.match(a)
    if m:
        return f"HLA-{m.group(1).upper()}*{m.group(2)}:{m.group(3)}"
    m = _MOUSE_II.match(a)
    if m:
        g = m.group(1)
        return "H2-I" + g[1].upper() + g[2:]
    return a


def is_class_ii(a):
    return bool(_II_SINGLE.match(a) or _II_PAIR.match(a) or _MOUSE_II.match(a)) or a.upper().startswith(("HLA-DR", "HLA-DQ", "HLA-DP", "H2-I"))


def mix2_name(a):
    """HLA-DQA1*05:01/DQB1*02:01 -> DQA1_05_01__DQB1_02_01 (MixMHC2pred's AlleleName); None if not HLA."""
    c = canonical_ii(a)
    if not (_II_SINGLE.match(c) or _II_PAIR.match(c)):
        return None
    return c[4:].replace("*", "_").replace(":", "_").replace("/", "__")


_REF27 = ["HLA-DRB1*01:01", "HLA-DRB1*03:01", "HLA-DRB1*04:01", "HLA-DRB1*04:05", "HLA-DRB1*07:01", "HLA-DRB1*08:02",
          "HLA-DRB1*09:01", "HLA-DRB1*11:01", "HLA-DRB1*12:01", "HLA-DRB1*13:02", "HLA-DRB1*15:01",
          "HLA-DRB3*01:01", "HLA-DRB3*02:02", "HLA-DRB4*01:01", "HLA-DRB5*01:01",
          "HLA-DQA1*01:01/DQB1*05:01", "HLA-DQA1*01:02/DQB1*06:02", "HLA-DQA1*03:01/DQB1*03:02",
          "HLA-DQA1*04:01/DQB1*04:02", "HLA-DQA1*05:01/DQB1*02:01", "HLA-DQA1*05:01/DQB1*03:01",
          "HLA-DPA1*01:03/DPB1*02:01", "HLA-DPA1*01:03/DPB1*04:01", "HLA-DPA1*02:01/DPB1*01:01",
          "HLA-DPA1*02:01/DPB1*05:01", "HLA-DPA1*02:01/DPB1*14:01", "HLA-DPA1*03:01/DPB1*04:02"]
_DR7 = ["HLA-DRB1*03:01", "HLA-DRB1*07:01", "HLA-DRB1*15:01", "HLA-DRB3*01:01", "HLA-DRB3*02:02",
        "HLA-DRB4*01:01", "HLA-DRB5*01:01"]

def pseudo_name(a):
    """Name in the NetMHCIIpan pseudosequence tables used by TLimmuno2 and ImmuScope: DRB1_0101, HLA-DQA10501-DQB10201,
    HLA-DPA10103-DPB10401, H-2-IAb. None when the allele is not a DR chain or a DQ/DP pair (the tables have no single DQ/DP chain)."""
    c = canonical_ii(a)
    m = re.match(r"^HLA-(DRB[1345])\*(\d{2,3}):(\d{2,3})$", c)
    if m:
        return f"{m.group(1)}_{m.group(2)}{m.group(3)}"
    m = _II_PAIR.match(c)
    if m:
        return f"HLA-{m.group(1).upper()}{m.group(2)}{m.group(3)}-{m.group(4).upper()}{m.group(5)}{m.group(6)}"
    if c.startswith("H2-I"):
        return "H-2-" + c[3:]
    return None


PRESETS_II = {
    "IEDB class II reference panel 27": {"human": True, "alleles": _REF27},
    "DR only (15 of the 27)": {"human": True, "alleles": [a for a in _REF27 if "DRB" in a]},
    "DR 7-allele set (CD4episcore)": {"human": True, "alleles": _DR7},
    "Mouse I-A / I-E": {"human": False, "alleles": ["H2-IAb", "H2-IAd", "H2-IAg7", "H2-IEd", "H2-IEk"]},
}
