#!/usr/bin/env python3
"""
Harvest real immunological data from the IEDB Query API (PostgREST).

Two datasets:
  1. BINDING   - peptide-MHC affinity, quantitative IC50/KD in nM  (mhc_search)
  2. IMMUNOGEN - T cell response outcome, Positive/Negative        (tcell_search)

Output: tidy CSV in ./data/
"""
import csv
import html
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://query-api.iedb.org"
OUT = Path(__file__).parent / "data"
OUT.mkdir(exist_ok=True)

PAGE = 5000          # rows per request
MAX_PAGES = 200      # safety stop
PAUSE = 0.34         # be polite to a public academic API

# HLA alleles to harvest. The first three match the C-ImmSim web form defaults.
ALLELES = [
    "HLA-A*02:01", "HLA-B*07:02", "HLA-DRB1*01:01",
    "HLA-A*01:01", "HLA-A*03:01", "HLA-A*11:01", "HLA-A*24:02",
    "HLA-B*08:01", "HLA-B*35:01", "HLA-B*44:02",
    "HLA-DRB1*04:01", "HLA-DRB1*07:01", "HLA-DRB1*15:01",
]

AA = set("ACDEFGHIKLMNPQRSTVWY")


def fetch(endpoint, params):
    """One GET against the query API. Returns parsed JSON list."""
    qs = urllib.parse.urlencode(params, quote_via=urllib.parse.quote)
    url = f"{API}/{endpoint}?{qs}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if attempt == 3:
                print(f"    FAILED {endpoint} {params.get('mhc_allele_name')}: {e}",
                      file=sys.stderr)
                return []
            time.sleep(2 * (attempt + 1))
    return []


# The API requires a stable sort key whenever offset paging is used.
ORDER_KEY = {"mhc_search": "elution_id.asc", "tcell_search": "tcell_id.asc",
             "bcell_search": "bcell_id.asc"}


def paginate(endpoint, params):
    """Yield all rows, following limit/offset until a short page comes back."""
    offset = 0
    for _ in range(MAX_PAGES):
        p = dict(params, limit=PAGE, offset=offset,
                 order=ORDER_KEY.get(endpoint, "structure_id.asc"))
        rows = fetch(endpoint, p)
        if not rows:
            return
        yield from rows
        if len(rows) < PAGE:
            return
        offset += PAGE
        time.sleep(PAUSE)


QUANT = re.compile(r"^\s*([<>=]*)\s*([0-9.eE+-]+)\s*(\S+)\s*$")


def parse_quant(raw):
    """'=&nbsp;123.4&nbsp;nM' -> ('=', 123.4, 'nM'). Returns None if unparseable."""
    if not raw:
        return None
    s = html.unescape(raw).replace("\xa0", " ").strip()
    m = QUANT.match(s)
    if not m:
        return None
    op, val, unit = m.groups()
    try:
        return (op or "="), float(val), unit
    except ValueError:
        return None


def clean_peptide(seq):
    """Keep only plain unmodified amino acid strings."""
    if not seq:
        return None
    s = seq.strip().upper()
    if not s or set(s) - AA:
        return None
    return s


def harvest_binding():
    """Quantitative peptide-MHC affinity, nM."""
    cols = ["linear_sequence", "mhc_allele_name", "mhc_class", "quantitative_measure",
            "qualitative_measure", "assay_names", "parent_source_antigen_name",
            "source_organism_name", "host_organism_name", "reference_id", "pubmed_id"]
    seen, out = set(), []

    for allele in ALLELES:
        n0 = len(out)
        params = {
            "select": ",".join(cols),
            "assay_names": "like.*nM*",
            "mhc_allele_name": f"eq.{allele}",
        }
        for r in paginate("mhc_search", params):
            pep = clean_peptide(r.get("linear_sequence"))
            if not pep or not (6 <= len(pep) <= 30):
                continue
            q = parse_quant(r.get("quantitative_measure"))
            if not q or q[2] != "nM":
                continue
            op, nm, _ = q
            if not (0 < nm < 1e7):
                continue
            key = (pep, allele, round(nm, 3), r.get("reference_id"))
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "peptide": pep,
                "length": len(pep),
                "allele": allele,
                "mhc_class": r.get("mhc_class"),
                "operator": op,
                "affinity_nM": nm,
                "qualitative": r.get("qualitative_measure"),
                "assay": (r.get("assay_names") or "").split("|")[0],
                "antigen": r.get("parent_source_antigen_name"),
                "source_organism": r.get("source_organism_name"),
                "host_organism": r.get("host_organism_name"),
                "pubmed_id": r.get("pubmed_id"),
            })
        print(f"  {allele:16s} {len(out) - n0:7d} binding records")

    return out


def harvest_immunogenicity():
    """T cell assay outcomes - the response actually observed in a host."""
    cols = ["linear_sequence", "mhc_allele_name", "mhc_class", "qualitative_measure",
            "assay_names", "parent_source_antigen_name", "source_organism_name",
            "host_organism_name", "disease_names", "reference_id", "pubmed_id"]
    seen, out = set(), []

    for allele in ALLELES:
        n0 = len(out)
        params = {
            "select": ",".join(cols),
            "mhc_allele_name": f"eq.{allele}",
            "qualitative_measure": "not.is.null",
        }
        for r in paginate("tcell_search", params):
            pep = clean_peptide(r.get("linear_sequence"))
            if not pep or not (6 <= len(pep) <= 30):
                continue
            qual = (r.get("qualitative_measure") or "").strip()
            if not qual:
                continue
            key = (pep, allele, qual, r.get("reference_id"))
            if key in seen:
                continue
            seen.add(key)
            out.append({
                "peptide": pep,
                "length": len(pep),
                "allele": allele,
                "mhc_class": r.get("mhc_class"),
                "qualitative": qual,
                "positive": int(qual.lower().startswith("positive")),
                "assay": (r.get("assay_names") or "").split("|")[0],
                "antigen": r.get("parent_source_antigen_name"),
                "source_organism": r.get("source_organism_name"),
                "host_organism": r.get("host_organism_name"),
                "disease": r.get("disease_names"),
                "pubmed_id": r.get("pubmed_id"),
            })
        print(f"  {allele:16s} {len(out) - n0:7d} T cell records")

    return out


def write_csv(rows, path):
    if not rows:
        print(f"  nothing to write for {path}", file=sys.stderr)
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {len(rows):,} rows -> {path}")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"

    if which in ("all", "binding"):
        print("BINDING (peptide-MHC affinity, nM)")
        write_csv(harvest_binding(), OUT / "binding_affinity.csv")

    if which in ("all", "tcell"):
        print("IMMUNOGENICITY (T cell assays)")
        write_csv(harvest_immunogenicity(), OUT / "tcell_immunogenicity.csv")
