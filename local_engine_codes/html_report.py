"""One self-contained HTML page for a finished run: regions with labels, a map of the protein, the best candidates,
the funnel, the minimum panel and exactly which rules and tools produced it. No scripts, no external files, so it can
be mailed or archived as it is.
"""
import html
import time

import pandas as pd

E = html.escape
TYPE_CLASS = {"CTL+HTL": "dual", "CTL only": "ctl", "HTL only": "htl"}
CSS = """
:root{--bg:#fff;--fg:#1c2330;--muted:#5b6678;--card:#f5f7fa;--line:#d5dbe5;--ctl:#2f6fdb;--htl:#d9822b;--dual:#1e8a4c}
@media (prefers-color-scheme:dark){:root{--bg:#141a26;--fg:#e6eaf2;--muted:#9aa6bb;--card:#1c2433;--line:#32405a;--ctl:#6b95ff;--htl:#f0a35a;--dual:#4ccf86}}
body{font:14px/1.5 system-ui,sans-serif;margin:0;background:var(--bg);color:var(--fg)}
main{max-width:1180px;margin:0 auto;padding:24px 16px 60px}
h1{font-size:22px;margin:0 0 2px}h2{font-size:16px;margin:28px 0 8px;border-bottom:1px solid var(--line);padding-bottom:4px}
.muted{color:var(--muted)}.kpis{display:flex;gap:10px;flex-wrap:wrap;margin:14px 0}
.kpi{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:8px 14px;min-width:130px}.kpi b{display:block;font-size:22px}.kpi span{color:var(--muted);font-size:12px}
table{border-collapse:collapse;font-size:13px;margin:6px 0}th,td{border:1px solid var(--line);padding:3px 8px;text-align:left;vertical-align:top}th{background:var(--card)}
.tag{display:inline-block;border-radius:10px;padding:0 8px;font-size:12px;color:#fff}.tag.dual{background:var(--dual)}.tag.ctl{background:var(--ctl)}.tag.htl{background:var(--htl)}
.seq{font-family:ui-monospace,monospace;word-break:break-all}.bar{display:inline-block;width:80px;height:8px;background:var(--line);border-radius:4px;vertical-align:middle;margin-right:5px;overflow:hidden}.bar i{display:block;height:100%;background:var(--ctl)}
.warn{color:#a86a00}.box{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:10px 14px}svg{max-width:100%;height:auto}
@media print{body{background:#fff;color:#000}}
"""


def _pct(x):
    return "" if x is None or pd.isna(x) else f"{100 * float(x):.1f}%"


def _num(x, d=3):
    return "" if x is None or (isinstance(x, float) and pd.isna(x)) else (f"{x:.{d}g}" if isinstance(x, float) else E(str(x)))


def _cell(v, f):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    return f(v) if f else _num(v)


def _table(df, cols, fmt=None):
    fmt = fmt or {}
    head = "".join(f"<th>{E(c)}</th>" for c in cols)
    body = "".join("<tr>" + "".join(f"<td>{_cell(r[c], fmt.get(c))}</td>" for c in cols) + "</tr>" for _, r in df.iterrows())
    return f"<table><tr>{head}</tr>{body}</table>"


def protein_map(regions, seqs):
    """Inline SVG: one lane per protein, regions drawn where they sit, coloured by type."""
    if regions.empty or not seqs:
        return ""
    W, pad, row = 1000, 90, 34
    out = [f'<svg viewBox="0 0 {W + pad + 10} {row * len(seqs) + 24}" xmlns="http://www.w3.org/2000/svg" font-size="11" fill="currentColor">']
    for k, (name, seq) in enumerate(seqs.items()):
        y = 12 + k * row
        sc = W / max(len(seq), 1)
        out.append(f'<text x="0" y="{y + 14}">{E(name[:14])}</text><rect x="{pad}" y="{y + 6}" width="{W}" height="12" fill="var(--line)" rx="3"/>')
        for _, r in regions[regions["Protein"] == name].iterrows():
            col = "var(--" + TYPE_CLASS[r["Type"]] + ")"
            x = pad + (r["Start"] - 1) * sc
            w = max((r["End"] - r["Start"] + 1) * sc, 2)
            out.append(f'<rect x="{x:.1f}" y="{y}" width="{w:.1f}" height="24" fill="{col}" opacity=".85"><title>{E(r["Region_ID"])} {E(r["Type"])} {r["Start"]}-{r["End"]}</title></rect>')
        out.append(f'<text x="{pad}" y="{y + 34}" fill="var(--muted)">1</text><text x="{pad + W}" y="{y + 34}" text-anchor="end" fill="var(--muted)">{len(seq)}</text>')
    out.append("</svg>")
    return "".join(out)


def build(run):
    """run: dict with title, info (inputs), regions (DataFrame), classes {1|2: {peptides DataFrame, consensus dict, details dict}},
    seqs {protein: seq}, warnings [str], population info."""
    reg, seqs, now = run["regions"], run.get("seqs", {}), time.strftime("%Y-%m-%d %H:%M")
    h = [f"<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>"
         f"<title>{E(run['title'])}</title><style>{CSS}</style></head><body><main>",
         f"<h1>{E(run['title'])}</h1><div class='muted'>{now} &middot; job {E(run['job'])} &middot; strategy <b>{E(run['strategy'])}</b></div>"]

    n_dual = int((reg["Type"] == "CTL+HTL").sum()) if len(reg) else 0
    k = [(len(reg), "epitope regions"), (n_dual, "CTL+HTL regions"), (int((reg["Type"] == "CTL only").sum()) if len(reg) else 0, "CTL only"),
         (int((reg["Type"] == "HTL only").sum()) if len(reg) else 0, "HTL only")]
    for c, tag in ((1, "CTL"), (2, "HTL")):
        res = run["classes"].get(c)
        if res and res["consensus"]["panel"]:
            p = res["consensus"]["panel"]
            k.append((_pct(p[-1]["cumulative_coverage"]), f"{tag} panel coverage ({len(p)} peptides)"))
    h.append("<div class='kpis'>" + "".join(f"<div class='kpi'><b>{E(str(v))}</b><span>{E(l)}</span></div>" for v, l in k) + "</div>")
    if run["warnings"]:
        h.append("<div class='box warn'>" + "<br>".join(E(w) for w in run["warnings"]) + "</div>")

    h.append("<h2>Epitope regions</h2>")
    h.append("<p class='muted'>A region is a stretch of the protein where predicted epitopes overlap. <b>CTL+HTL</b>: both a class I (CD8) and a class II (CD4) candidate sit there. "
             "Ranked CTL+HTL first, then by region score (mean priority score of the best peptides, +10 for CTL+HTL). Coverage is the share of the population presenting the best peptide.</p>")
    h.append(protein_map(reg, seqs))
    if len(reg):
        show = reg.head(60).copy()
        tagf = lambda v: f"<span class='tag {TYPE_CLASS.get(v, '')}'>{E(v)}</span>"
        seqf = lambda v: f"<span class='seq'>{E(v)}</span>"
        cov = lambda v: "" if pd.isna(v) else f"<span class='bar'><i style='width:{100 * float(v):.0f}%'></i></span>{_pct(v)}"
        h.append(_table(show, ["Region_ID", "Type", "Protein", "Start", "End", "Sequence", "Region_Score", "CTL_Best_Peptide", "CTL_Coverage",
                               "HTL_Best_Peptide", "HTL_Coverage", "In_Panel", "Synthesis_Notes"],
                        {"Type": tagf, "Sequence": seqf, "CTL_Best_Peptide": seqf, "HTL_Best_Peptide": seqf, "CTL_Coverage": cov, "HTL_Coverage": cov,
                         "In_Panel": lambda v: "yes" if v is True or v == "True" else ""}))
        if len(reg) > 60:
            h.append(f"<p class='muted'>Top 60 of {len(reg)} regions shown; all are in regions.csv.</p>")
    else:
        h.append("<p>No region passed the rules. See the suggestions below.</p>")

    for c, tag, title in ((1, "CTL", "MHC class I (CTL)"), (2, "HTL", "MHC class II (HTL)")):
        res = run["classes"].get(c)
        if not res:
            continue
        cons, pep = res["consensus"], res["peptides"]
        h.append(f"<h2>{title}</h2>")
        h.append("<div class='box'>" + E(", ".join(t["label"] for t in cons["tools"])) + f"<br><span class='muted'>{len(res['details']['input']['alleles'])} alleles, "
                 f"{res['details']['input']['peptides']:,} unique peptides, coverage by {E(cons['coverage_metric'])}</span></div>")
        h.append("<h3>Funnel</h3><table>" + "".join(f"<tr><td>{E(f['label'])}</td><td style='text-align:right'>{f['n']:,}</td></tr>" for f in cons["funnel"]) + "</table>")
        top = pep[pep["Is_Representative"].astype(bool)].head(25)
        if len(top):
            cols = [c2 for c2 in ["Priority_Rank", "Peptide", "Region_ID", "Region", "Best_Allele", "N_Alleles_Candidate", "Coverage_Candidate",
                                  "Agree_Fraction_PopW", "Host_Self", "Bind_Consensus_Best", "Imm_TopPct_Best", "Cluster_Size", "Cluster_Members", "In_Panel"] if c2 in top.columns]
            h.append("<h3>Best candidates</h3>" + _table(top, cols, {"Peptide": lambda v: f"<span class='seq'>{E(v)}</span>", "Coverage_Candidate": _pct,
                                                                   "Agree_Fraction_PopW": _pct, "Cluster_Members": lambda v: f"<span class='seq'>{E(v if len(v) < 120 else v[:117] + '...')}</span>", "In_Panel": lambda v: "yes" if v is True or v == "True" else ""}))
        if cons["panel"]:
            h.append("<h3>Minimum panel</h3>" + _table(pd.DataFrame(cons["panel"]).assign(new_alleles=lambda d: d["new_alleles"].map(", ".join)),
                                                      ["order", "peptide", "region", "gain", "cumulative_coverage", "new_alleles"],
                                                      {"peptide": lambda v: f"<span class='seq'>{E(v)}</span>", "gain": _pct, "cumulative_coverage": _pct}))
        sg = cons["suggestions"]
        if sg.get("advice"):
            h.append("<h3>Suggestions</h3>" + "".join(f"<div>&bull; {E(a)}</div>" for a in sg["advice"]))
        if sg.get("alleles", {}).get("add"):
            h.append("<div class='muted'>Common alleles not selected: " + E(", ".join(f"{a['allele']} ({_pct(a['frequency'])})" for a in sg["alleles"]["add"])) + "</div>")

    h.append("<h2>How this was produced</h2><div class='box'><b>Input</b><br>" + E(run["info_text"]) + "<br><br><b>Rules</b><br>"
             + "<br>".join(f"{E(k2)} = {E(str(v))}" for k2, v in run["rules"].items() if k2 != "rank_weights")
             + f"<br>priority weights: {E(str(run['rules'].get('rank_weights')))}<br><br><b>Population</b><br>{E(run['population_text'])}</div>")
    h.append("<p class='muted'>Votes: percentile tools vote at their cut-off; BigMHC-EL at its score cut-off; immunogenicity tools without an absolute scale vote for the top fraction of plausible binders. "
             "Agreement = votes / tools that scored the pair. Coverage uses Hardy-Weinberg per region with loci treated as independent. The numbers are predictions, not measurements.</p>")
    h.append("</main></body></html>")
    return "".join(h)
