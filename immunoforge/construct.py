"""Multi-epitope vaccine construct: assembly, properties, junction screening, VaxiJen.

Pure functions over strings and tables, like comparative.py: no GUI, no jobs, so the same
code backs the web tab, the headless runner and the tests.

The part that is not in the usual assembler is **junction screening**. Joining epitopes with
linkers creates peptides that exist in neither the epitopes nor the linkers: every window
that straddles a seam is new sequence, invented by the assembly. Those junctional peptides
can be immunogenic (competing with the epitopes you chose), allergenic, or toxic, and they
are routinely left unchecked. `junctions()` enumerates exactly the windows that cross a seam
so the engine's own predictors can score them, which is the whole reason this lives inside
the engine rather than in a separate tool.

Physicochemistry follows the published definitions, with the references in PROPERTY_NOTES.
"""
import re

AA = set("ACDEFGHIKLMNPQRSTVWY")

# ---------------------------------------------------------------- linkers
#
# Sequence, what it joins, why, and the paper it comes from. The notes are shown in the page,
# so a construct can be justified without leaving the tool.
LINKERS = {
    "adjuvant_to_epitope": ("EAAAK", "adjuvant -> first epitope",
                            "Rigid alpha-helix; keeps the adjuvant domain from interfering with the epitope region.",
                            "Arai R et al. (2001) Protein Eng 14(8):529-532"),
    "CTL_intra":           ("AAY", "between CTL epitopes",
                            "Proteasomal cleavage signal; Y is favoured at P1 of the chymotryptic-like site, so the "
                            "upstream epitope is liberated cleanly.",
                            "Tenzer S et al. (2005) Nat Immunol 6(9):928-937"),
    "HTL_intra":           ("GPGPG", "between HTL epitopes",
                            "Flexible, breaks secondary structure and keeps the MHC-II groove accessible.",
                            "Ojha R, Prajapati VK (2023) Front Immunol 14:1101679"),
    "Bcell_intra":         ("KK", "between B-cell epitopes",
                            "Charged dibasic spacer; favours an extended, surface-accessible conformation. Also a "
                            "cathepsin B cleavage site.",
                            "Nezafat N et al. (2016) Chem Biol Drug Des 87(5):678-686"),
    "block_to_block":      ("GGGS", "between epitope blocks",
                            "Short flexible glycine-serine junction between structurally different blocks.",
                            "Chen X, Zaro JL, Shen WC (2013) Adv Drug Deliv Rev 65(10):1357-1369"),
    "epitope_to_Ctag":     ("GGGS", "last epitope -> C-terminal tag",
                            "Flexible junction before the purification tag.",
                            "Chen X, Zaro JL, Shen WC (2013) Adv Drug Deliv Rev 65(10):1357-1369"),
}

ADJUVANTS = {
    "RS09":        ("APPHALS", "TLR4 agonist peptide", "Shanmugam A et al. (2012) PLoS One 7(2):e30839"),
    "PADRE":       ("AKFVAAWTLKAAA", "universal pan-DR helper epitope", "Alexander J et al. (1994) Immunity 1(9):751-761"),
    "beta-defensin-3": ("GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK",
                        "human beta-defensin 3; TLR1/2 agonist, links innate and adaptive",
                        "Ferris LK et al. (2013) J Invest Dermatol 133(2):460-468"),
    "flagellin-D1": ("MAQVINTNSLSLLTQNNLNKSQSALGTAIERLSSGLRINSAKDDAAGQAIANRFTANIKGLTQASRNANDGISIAQTTEGALNEINNNLQRVRELAVQSAN",
                     "Salmonella flagellin D1 domain; TLR5 agonist", "Mizel SB, Bates JT (2010) J Immunol 185(10):5677-5682"),
}

CTAGS = {
    "His6":    ("HHHHHH", "Ni-NTA affinity purification"),
    "His8":    ("HHHHHHHH", "Ni-NTA, tighter binding"),
    "Strep-II": ("WSHPQFEK", "gentle elution; good for labile proteins"),
    "FLAG":    ("DYKDDDDK", "high-affinity antibody available"),
}

# ---------------------------------------------------------------- physicochemistry

# free amino acid molecular weights (Da); a peptide loses one water per bond
AA_MW = {'A': 89.09, 'R': 174.20, 'N': 132.12, 'D': 133.10, 'C': 121.16, 'E': 147.13,
         'Q': 146.15, 'G': 75.07, 'H': 155.16, 'I': 131.17, 'L': 131.17, 'K': 146.19,
         'M': 149.21, 'F': 165.19, 'P': 115.13, 'S': 105.09, 'T': 119.12, 'W': 204.23,
         'Y': 181.19, 'V': 117.15}
WATER = 18.02

AA_PKA = {'N_term': 8.00, 'C_term': 3.10, 'D': 3.65, 'E': 4.25, 'H': 6.00,
          'C': 8.18, 'Y': 10.07, 'K': 10.53, 'R': 12.48}

KD = {'I': 4.5, 'V': 4.2, 'L': 3.8, 'F': 2.8, 'C': 2.5, 'M': 1.9, 'A': 1.8, 'G': -0.4,
      'T': -0.7, 'S': -0.8, 'W': -0.9, 'Y': -1.3, 'P': -1.6, 'H': -3.2, 'D': -3.5,
      'E': -3.5, 'N': -3.5, 'Q': -3.5, 'K': -3.9, 'R': -4.5}

# Dipeptide Instability Weight Values, all 400 pairs as published.
# Guruprasad K, Reddy BVB, Pandit MW (1990) Protein Eng 4(2):155-161.
DIWV = {
    "AA": 1.0, "AC": 44.94, "AD": -7.49, "AE": 1.0, "AF": 1.0, "AG": 1.0, "AH": 1.0, "AI": 1.0, "AK": 1.0, "AL":
    1.0, "AM": 1.0, "AN": 1.0, "AP": 20.26, "AQ": 1.0, "AR": 1.0, "AS": 1.0, "AT": 1.0, "AV": 1.0, "AW": 1.0,
    "AY": 1.0, "CA": 1.0, "CC": 1.0, "CD": 1.0, "CE": 1.0, "CF": 1.0, "CG": 1.0, "CH": 1.0, "CI": 1.0, "CK":
    65.17, "CL": 1.0, "CM": 33.6, "CN": 1.0, "CP": 1.0, "CQ": 1.0, "CR": 1.0, "CS": 1.0, "CT": 33.6, "CV": 1.0,
    "CW": 24.68, "CY": 1.0, "DA": 1.0, "DC": 1.0, "DD": 1.0, "DE": 1.0, "DF": 1.0, "DG": 1.0, "DH": 1.0, "DI":
    1.0, "DK": 1.0, "DL": 1.0, "DM": 1.0, "DN": 1.0, "DP": 1.0, "DQ": 1.0, "DR": 1.0, "DS": 20.26, "DT": 1.0,
    "DV": 1.0, "DW": 1.0, "DY": 1.0, "EA": 1.0, "EC": 44.94, "ED": 1.0, "EE": 1.0, "EF": 1.0, "EG": 1.0, "EH":
    1.0, "EI": 20.26, "EK": 1.0, "EL": 1.0, "EM": 1.0, "EN": 1.0, "EP": 20.26, "EQ": 20.26, "ER": 1.0, "ES":
    1.0, "ET": 1.0, "EV": 1.0, "EW": 1.0, "EY": 1.0, "FA": 1.0, "FC": 1.0, "FD": 13.34, "FE": 1.0, "FF": 1.0,
    "FG": 1.0, "FH": 1.0, "FI": 1.0, "FK": 1.0, "FL": 1.0, "FM": 1.0, "FN": 1.0, "FP": 20.26, "FQ": 1.0, "FR":
    1.0, "FS": 1.0, "FT": 1.0, "FV": 1.0, "FW": 1.0, "FY": 33.6, "GA": 1.0, "GC": 1.0, "GD": 1.0, "GE": 1.0,
    "GF": 1.0, "GG": 13.34, "GH": 1.0, "GI": 1.0, "GK": 1.0, "GL": 1.0, "GM": 1.0, "GN": 1.0, "GP": 1.0, "GQ":
    1.0, "GR": 1.0, "GS": 1.0, "GT": -7.49, "GV": 1.0, "GW": 13.34, "GY": 1.0, "HA": 1.0, "HC": 1.0, "HD": 1.0,
    "HE": 1.0, "HF": 1.0, "HG": 1.0, "HH": 1.0, "HI": 44.94, "HK": 1.0, "HL": 1.0, "HM": 1.0, "HN": 1.0, "HP":
    1.0, "HQ": 1.0, "HR": 1.0, "HS": 1.0, "HT": 1.0, "HV": 1.0, "HW": 1.0, "HY": 44.94, "IA": 1.0, "IC": 1.0,
    "ID": 1.0, "IE": 44.94, "IF": 1.0, "IG": 1.0, "IH": 44.94, "II": 1.0, "IK": 1.0, "IL": 1.0, "IM": 1.0, "IN":
    1.0, "IP": 1.0, "IQ": 1.0, "IR": 1.0, "IS": 1.0, "IT": 1.0, "IV": 1.0, "IW": 1.0, "IY": 1.0, "KA": 1.0,
    "KC": 1.0, "KD": 1.0, "KE": 1.0, "KF": 1.0, "KG": -7.49, "KH": 1.0, "KI": 1.0, "KK": 1.0, "KL": 1.0, "KM":
    33.6, "KN": 1.0, "KP": 1.0, "KQ": 24.68, "KR": 33.6, "KS": 1.0, "KT": 1.0, "KV": 1.0, "KW": 1.0, "KY": 1.0,
    "LA": 1.0, "LC": 1.0, "LD": 1.0, "LE": 1.0, "LF": 1.0, "LG": 1.0, "LH": 1.0, "LI": 1.0, "LK": 1.0, "LL":
    1.0, "LM": 1.0, "LN": 1.0, "LP": 20.26, "LQ": 33.6, "LR": 20.26, "LS": 1.0, "LT": 1.0, "LV": 1.0, "LW":
    24.68, "LY": 1.0, "MA": 13.34, "MC": 1.0, "MD": 1.0, "ME": 1.0, "MF": 1.0, "MG": 1.0, "MH": 58.28, "MI":
    1.0, "MK": 1.0, "ML": 1.0, "MM": 1.0, "MN": 1.0, "MP": 44.94, "MQ": 1.0, "MR": 1.0, "MS": 1.0, "MT": 1.0,
    "MV": 1.0, "MW": 1.0, "MY": 24.68, "NA": 1.0, "NC": 1.0, "ND": 1.0, "NE": 1.0, "NF": 1.0, "NG": -14.03,
    "NH": 1.0, "NI": 44.94, "NK": 1.0, "NL": 1.0, "NM": 1.0, "NN": 1.0, "NP": 1.0, "NQ": 1.0, "NR": 1.0, "NS":
    1.0, "NT": -14.03, "NV": 1.0, "NW": -9.37, "NY": 1.0, "PA": 20.26, "PC": 1.0, "PD": -6.54, "PE": 18.38,
    "PF": 1.0, "PG": 1.0, "PH": 1.0, "PI": 1.0, "PK": 1.0, "PL": 1.0, "PM": 1.0, "PN": 1.0, "PP": 1.0, "PQ":
    20.26, "PR": 1.0, "PS": 20.26, "PT": 1.0, "PV": 1.0, "PW": 1.0, "PY": 1.0, "QA": 1.0, "QC": 1.0, "QD": 1.0,
    "QE": 20.26, "QF": 1.0, "QG": 1.0, "QH": 1.0, "QI": 1.0, "QK": 1.0, "QL": 1.0, "QM": 1.0, "QN": 1.0, "QP":
    20.26, "QQ": 20.26, "QR": 1.0, "QS": 44.94, "QT": 1.0, "QV": 1.0, "QW": 1.0, "QY": 33.6, "RA": 1.0, "RC":
    1.0, "RD": 1.0, "RE": 1.0, "RF": 1.0, "RG": 1.0, "RH": 1.0, "RI": 1.0, "RK": 1.0, "RL": 1.0, "RM": 1.0,
    "RN": 1.0, "RP": 20.26, "RQ": 20.26, "RR": 1.0, "RS": 44.94, "RT": 1.0, "RV": 1.0, "RW": 1.0, "RY": 1.0,
    "SA": 1.0, "SC": 33.6, "SD": 1.0, "SE": 20.26, "SF": 1.0, "SG": 1.0, "SH": 1.0, "SI": 1.0, "SK": 1.0, "SL":
    1.0, "SM": 1.0, "SN": 1.0, "SP": 44.94, "SQ": 20.26, "SR": 1.0, "SS": 20.26, "ST": 1.0, "SV": 1.0, "SW":
    1.0, "SY": 1.0, "TA": 1.0, "TC": 1.0, "TD": 1.0, "TE": 20.26, "TF": 1.0, "TG": -7.49, "TH": 1.0, "TI": 1.0,
    "TK": 1.0, "TL": 1.0, "TM": 1.0, "TN": 1.0, "TP": 1.0, "TQ": 1.0, "TR": 1.0, "TS": 1.0, "TT": 1.0, "TV":
    1.0, "TW": 1.0, "TY": 1.0, "VA": 1.0, "VC": 1.0, "VD": 1.0, "VE": 1.0, "VF": 1.0, "VG": -7.49, "VH": 1.0,
    "VI": 1.0, "VK": 1.0, "VL": 1.0, "VM": 1.0, "VN": 1.0, "VP": 20.26, "VQ": 1.0, "VR": 1.0, "VS": 1.0, "VT":
    1.0, "VV": 1.0, "VW": 1.0, "VY": 1.0, "WA": 1.0, "WC": 1.0, "WD": 1.0, "WE": 1.0, "WF": 1.0, "WG": 13.34,
    "WH": 1.0, "WI": 1.0, "WK": 1.0, "WL": 13.34, "WM": 24.68, "WN": 13.34, "WP": 1.0, "WQ": 1.0, "WR": 1.0,
    "WS": 1.0, "WT": 1.0, "WV": 1.0, "WW": 1.0, "WY": 1.0, "YA": 1.0, "YC": 1.0, "YD": 1.0, "YE": 1.0, "YF":
    1.0, "YG": -6.54, "YH": 13.34, "YI": 1.0, "YK": 1.0, "YL": 1.0, "YM": 44.94, "YN": 1.0, "YP": 13.34, "YQ":
    1.0, "YR": 1.0, "YS": 1.0, "YT": 1.0, "YV": 1.0, "YW": 1.0, "YY": 1.0
}

PROPERTY_NOTES = {
    "mw": "Sum of residue weights minus one water per peptide bond. Fasman GD (1989) CRC Handbook.",
    "pi": "pH where the net charge is zero, by bisection on the Henderson-Hasselbalch sum. "
          "Free-amino-acid pKa values, so it is accurate to about 0.1-0.2 pH units; a folded "
          "protein shifts individual pKa values. Bjellqvist B et al. (1993) Electrophoresis 14:1023.",
    "gravy": "Mean Kyte-Doolittle hydropathy. Below 0 is hydrophilic and usually what a soluble "
             "construct wants. Kyte J, Doolittle RF (1982) J Mol Biol 157:105.",
    "instability": "(10/N) x sum of dipeptide instability weights. Under 40 predicts stability in "
                   "vitro. Trained on 44 proteins, so it is a heuristic: disordered and short "
                   "sequences often exceed 40 while being perfectly usable. "
                   "Guruprasad K et al. (1990) Protein Eng 4:155.",
    "aliphatic": "Relative volume of aliphatic side chains (A, V, I, L). Higher means more "
                 "thermostable; above about 80 is comfortable. Ikai A (1980) J Biochem 88:1895.",
}


def mw(seq):
    if not seq:
        return 0.0
    return round(sum(AA_MW.get(a, 110.0) for a in seq) - WATER * (len(seq) - 1), 2)


def pi(seq):
    """Isoelectric point by bisection on the net-charge curve."""
    if not seq:
        return 0.0
    n = {a: seq.count(a) for a in "DEHCYKR"}

    def charge(ph):
        q = 1.0 / (1.0 + 10 ** (ph - AA_PKA["N_term"]))
        q -= 1.0 / (1.0 + 10 ** (AA_PKA["C_term"] - ph))
        for a in ("K", "R", "H"):
            q += n.get(a, 0) / (1.0 + 10 ** (ph - AA_PKA[a]))
        for a in ("D", "E", "C", "Y"):
            q -= n.get(a, 0) / (1.0 + 10 ** (AA_PKA[a] - ph))
        return q

    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        if charge(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2, 2)


def gravy(seq):
    return round(sum(KD.get(a, 0.0) for a in seq) / len(seq), 3) if seq else 0.0


def aliphatic_index(seq):
    """Ikai's aliphatic index: 100 x (A% + 2.9 V% + 3.9 (I% + L%))."""
    if not seq:
        return 0.0
    n = len(seq)
    a, v, i, l = (seq.count(x) / n for x in "AVIL")
    return round(100 * (a + 2.9 * v + 3.9 * (i + l)), 2)


def instability(seq, diwv):
    """Guruprasad instability index. `diwv` is the 400-entry dipeptide weight table."""
    if len(seq) < 2:
        return 0.0
    return round((10.0 / len(seq)) * sum(diwv.get(seq[i:i + 2], 1.0) for i in range(len(seq) - 1)), 2)


def properties(seq, diwv=None):
    out = {"length": len(seq), "mw": mw(seq), "pi": pi(seq), "gravy": gravy(seq),
           "aliphatic_index": aliphatic_index(seq)}
    if diwv:
        out["instability"] = instability(seq, diwv)
    return out


def verdicts(p):
    """The design checklist, as (ok, message) pairs. Thresholds from the usual design papers."""
    out = []
    if "instability" in p:
        out.append((p["instability"] < 40,
                    f"Instability index {p['instability']}: "
                    + ("under 40, predicted stable in vitro" if p["instability"] < 40
                       else "40 or more, may degrade; re-ordering epitopes or changing linkers often fixes it")))
    out.append((p["gravy"] < 0, f"GRAVY {p['gravy']}: "
                + ("negative, good aqueous solubility expected" if p["gravy"] < 0
                   else "positive, may be hard to express solubly")))
    out.append((4.5 <= p["pi"] <= 9.5, f"pI {p['pi']}: "
                + ("in the workable 4.5-9.5 range" if 4.5 <= p["pi"] <= 9.5
                   else "extreme; solubility is worst at the pI, so keep buffers away from it")))
    out.append((p["aliphatic_index"] >= 80, f"Aliphatic index {p['aliphatic_index']}: "
                + ("80 or above, thermostable" if p["aliphatic_index"] >= 80 else "below 80, less thermostable")))
    out.append((20 <= p["mw"] / 1000 <= 150, f"Molecular weight {p['mw'] / 1000:.1f} kDa: "
                + ("in the usual 20-150 kDa range" if 20 <= p["mw"] / 1000 <= 150
                   else "outside the usual 20-150 kDa range; check the construct is complete")))
    return out


# ---------------------------------------------------------------- assembly

PREFERRED_C_TERM = set("LVIFMA")      # favoured at P1 by the proteasome's chymotryptic-like site


def order_ctl(epitopes):
    """CTL epitopes ending in a favoured proteasomal P1 residue first, order otherwise kept."""
    return sorted(epitopes, key=lambda e: (not e or e[-1] not in PREFERRED_C_TERM))


def assemble(ctl=(), htl=(), bcell=(), adjuvant="", ctag="", optimise_ctl=True, linkers=None):
    """Build the construct. Returns (sequence, parts) where parts are the labelled pieces in order.

    Each part is {kind, label, seq, start, end} with 1-based inclusive coordinates, so the page
    can colour the construct and every downstream table can say which piece a position is in.
    """
    L = {k: (linkers or {}).get(k, v[0]) for k, v in LINKERS.items()}
    ctl = order_ctl([e for e in ctl if e]) if optimise_ctl else [e for e in ctl if e]
    htl = [e for e in htl if e]
    bcell = [e for e in bcell if e]

    parts, blocks = [], []
    for kind, eps, intra in (("CTL", ctl, "CTL_intra"), ("HTL", htl, "HTL_intra"), ("Bcell", bcell, "Bcell_intra")):
        if not eps:
            continue
        block = []
        for i, e in enumerate(eps, 1):
            if i > 1:
                block.append(("linker", intra, L[intra]))
            block.append(("epitope", f"{kind} {i}", e))
        blocks.append((kind, block))

    seq_parts = []
    if adjuvant:
        seq_parts.append(("adjuvant", "adjuvant", adjuvant))
        if blocks or ctag:
            seq_parts.append(("linker", "adjuvant_to_epitope", L["adjuvant_to_epitope"]))
    for i, (_kind, block) in enumerate(blocks):
        if i:
            seq_parts.append(("linker", "block_to_block", L["block_to_block"]))
        seq_parts.extend(block)
    if ctag:
        if blocks or adjuvant:
            seq_parts.append(("linker", "epitope_to_Ctag", L["epitope_to_Ctag"]))
        seq_parts.append(("tag", "C-terminal tag", ctag))

    seq, pos = [], 1
    for kind, label, s in seq_parts:
        parts.append({"kind": kind, "label": label, "seq": s, "start": pos, "end": pos + len(s) - 1})
        seq.append(s)
        pos += len(s)
    return "".join(seq), parts


# ---------------------------------------------------------------- junction screening

def junctions(parts, lengths=(9, 15)):
    """Peptides created by the joining itself: every window that crosses a seam.

    A window lying entirely inside one epitope is something you already chose and already
    screened. A window spanning a seam is new sequence that the assembly invented, and it is
    what has to be checked for unintended immunogenicity, allergenicity and toxicity.

    Returns [{peptide, length, start, end, spans, seam}] de-duplicated by (peptide, start).
    """
    seq = "".join(p["seq"] for p in parts)
    if not seq or not parts:
        return []
    # which part each 0-based position belongs to
    owner = []
    for i, p in enumerate(parts):
        owner.extend([i] * len(p["seq"]))

    # Sequence that is already in the construct as a part is not new, however the window
    # falls: an epitope flanked by a linker that repeats its residues (AAAAAAAAA + AAY) can
    # be re-spelled by a window that crosses the seam. It was chosen and screened already.
    known = {p["seq"] for p in parts if p["kind"] != "linker"}

    out, seen = [], set()
    for k in sorted({int(x) for x in lengths if int(x) > 1}):
        for s in range(0, len(seq) - k + 1):
            span = owner[s:s + k]
            if span[0] == span[-1]:
                continue                                    # inside a single part
            idx = sorted(set(span))
            if all(parts[i]["kind"] == "linker" for i in idx):
                continue                                    # pure linker, nothing new to present
            if seq[s:s + k] in known:
                continue                                    # not new sequence
            key = (seq[s:s + k], s)
            if key in seen:
                continue
            seen.add(key)
            names = [parts[i]["label"] for i in idx]
            out.append({"peptide": seq[s:s + k], "length": k, "start": s + 1, "end": s + k,
                        "spans": names, "seam": " | ".join(names)})
    return out


def junction_summary(parts, lengths=(9, 15)):
    j = junctions(parts, lengths)
    return {"n": len(j), "unique_peptides": len({x["peptide"] for x in j}),
            "by_length": {k: sum(1 for x in j if x["length"] == k) for k in sorted({x["length"] for x in j})}}


# ---------------------------------------------------------------- VaxiJen (online only)

def to_fasta(seqs, prefix="seq"):
    """Peptides (or (name, seq) pairs) -> FASTA text, for pasting into a web service."""
    lines = []
    for i, item in enumerate(seqs, 1):
        name, s = item if isinstance(item, (tuple, list)) else (f"{prefix}_{i}", item)
        s = str(s).strip().upper()
        if not s:
            continue
        lines.append(f">{name}")
        lines.extend(s[j:j + 60] for j in range(0, len(s), 60))
    return "\n".join(lines) + ("\n" if lines else "")


# VaxiJen's wording has changed between versions and between the protein and tumour models,
# so several shapes are accepted rather than one. Score may be a bare number or a percentage.
_VJ_PATTERNS = (
    # "Overall Prediction for the Protective Antigen = 0.4521 ( Probable NON-ANTIGEN )"
    re.compile(r"(?:=|:)\s*([-+]?\d*\.?\d+)\s*\(\s*Probable\s+(ANTIGEN|NON-?ANTIGEN|IMMUNOGEN|NON-?IMMUNOGEN)",
               re.I),
    # "Probable IMMUNOGEN with a probability of 87%"
    re.compile(r"Probable\s+(ANTIGEN|NON-?ANTIGEN|IMMUNOGEN|NON-?IMMUNOGEN)\s*,?\s*with a probability of\s*"
               r"([-+]?\d*\.?\d+)\s*%?", re.I),
    # bare "... = 0.4521" with the verdict elsewhere on the line
    re.compile(r"Probable\s+(ANTIGEN|NON-?ANTIGEN|IMMUNOGEN|NON-?IMMUNOGEN)", re.I),
)
_VJ_NAME = re.compile(r">\s*(\S+)|\bseq[_\s]?(\d+)\b", re.I)


def parse_vaxijen(text, threshold=0.4):
    """VaxiJen output (pasted or saved) -> [{name, antigen, score, verdict}] plus warnings.

    Handles the formats VaxiJen has used: "= 0.4521 ( Probable ANTIGEN )", "Probable IMMUNOGEN
    with a probability of 87%", and a verdict with no number. A percentage is divided by 100 so
    `score` is always on VaxiJen's 0-1 scale, comparable with the threshold (0.4 for the default
    protein model; 0.5 for tumour and some virus models, so it is a parameter).
    """
    warns, rows = [], []
    blocks = re.split(r"(?=^\s*>)", text or "", flags=re.M)
    if len(blocks) <= 1:                                     # no FASTA headers: split per line
        blocks = [b for b in (text or "").splitlines() if b.strip()]

    auto = 0
    for block in blocks:
        if not block.strip():
            continue
        score, label = None, None
        for pat in _VJ_PATTERNS:
            m = pat.search(block)
            if not m:
                continue
            groups = [g for g in m.groups() if g is not None]
            for g in groups:
                if re.fullmatch(r"[-+]?\d*\.?\d+", g):
                    score = float(g)
                    if "%" in block[m.start():m.end() + 1] or score > 1.0:
                        score /= 100.0
                else:
                    label = g.upper().replace("NONANTIGEN", "NON-ANTIGEN").replace("NONIMMUNOGEN", "NON-IMMUNOGEN")
            break
        if label is None and score is None:
            continue
        nm = _VJ_NAME.search(block)
        if nm and nm.group(1):
            name = nm.group(1)
        elif nm and nm.group(2):
            name = f"seq_{nm.group(2)}"
        else:
            auto += 1
            name = f"seq_{auto}"
        antigen = None
        if label:
            antigen = not label.startswith("NON")
        elif score is not None:
            antigen = score >= threshold
        rows.append({"name": name, "antigen": bool(antigen),
                     "score": None if score is None else round(score, 4),
                     "verdict": label or ("ANTIGEN" if antigen else "NON-ANTIGEN")})
    if not rows:
        warns.append("Nothing recognisable in that text. Paste VaxiJen's result page, including the "
                     "'Probable ANTIGEN / NON-ANTIGEN' lines.")
    missing = [r["name"] for r in rows if r["score"] is None]
    if missing:
        warns.append(f"{len(missing)} entr(y/ies) had a verdict but no score; the threshold could not be "
                     f"re-applied to them ({', '.join(missing[:4])}).")
    return rows, warns


# ---------------------------------------------------------------- population coverage

def parse_epitope_alleles(text):
    """Lines of "PEPTIDE  allele, allele" -> ({peptide: [alleles]}, [peptides in order], warnings).

    The alleles may be left off, which is normal when the epitopes came from the literature
    rather than from a run: the peptide still goes into the construct, it just cannot
    contribute to coverage, and the page says how many are in that position.
    """
    mapping, order, warns = {}, [], []
    for n, line in enumerate((text or "").splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[\s,;]+", line, maxsplit=1)
        pep = parts[0].strip().upper()
        if not pep:
            continue
        bad = set(pep) - AA
        if bad:
            warns.append(f"line {n}: '{pep[:18]}' has non-standard residues {''.join(sorted(bad))}; skipped")
            continue
        alleles = []
        if len(parts) > 1:
            alleles = [a.strip() for a in re.split(r"[,;|]| and ", parts[1]) if a.strip()]
        if pep in mapping:
            mapping[pep] = sorted(set(mapping[pep]) | set(alleles))
            continue
        mapping[pep] = alleles
        order.append(pep)
    return mapping, order, warns


def coverage_report(mapping, pop):
    """Population coverage of a construct, per region, the way IEDB's tool reports it.

    A region's coverage is the share of people carrying at least one allele that presents at
    least one epitope in the construct, so the union of all the epitopes' alleles is what
    counts. Per-epitope coverage is reported beside it to show which ones carry the construct
    and which are only there for one region.
    """
    if pop is None:
        return {"error": "No population file is loaded, so coverage cannot be computed."}
    all_alleles = sorted({a for alleles in mapping.values() for a in alleles})
    typed = [a for a in all_alleles if pop.resolve(a).source != "none"]
    untyped = [a for a in all_alleles if a not in typed]
    out = {"regions": [], "national": 0.0, "alleles": all_alleles, "alleles_without_data": untyped,
           "n_epitopes": len(mapping),
           "n_epitopes_without_alleles": sum(1 for v in mapping.values() if not v)}
    if not typed:
        out["error"] = ("None of the alleles given have frequency data, so coverage is 0. "
                        "HLA-DP and mouse alleles are not typed in this file.")
        return out
    per = pop.per_region(typed)
    out["national"] = round(float(per.get("national", 0.0)), 4)
    # float(): these come back as numpy scalars, which json cannot serialise
    out["regions"] = [{"region": r, "coverage": round(float(v), 4)}
                      for r, v in per.items() if r != "national"]
    out["per_epitope"] = sorted(
        ({"peptide": p, "alleles": a,
          "national": round(float(pop.per_region([x for x in a if pop.resolve(x).source != "none"]).get("national", 0.0)), 4)
          if any(pop.resolve(x).source != "none" for x in a) else 0.0}
         for p, a in mapping.items()),
        key=lambda d: -d["national"])
    # what one more epitope could add: the biggest single gain from an allele not yet covered
    out["weakest_region"] = min(out["regions"], key=lambda d: d["coverage"]) if out["regions"] else None
    return out
