"""BigMHC (EL + IM) Flask page. Standalone: python3 bigmhc_app.py; also mounted by engine_server.py."""
import os, sys

from flask import Flask, jsonify, request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import bigmhc_runner as bm                                  # noqa: E402
from common import alleles as allele_lib, config            # noqa: E402

app = Flask(__name__)


@app.get("/api/presets")
def api_presets():
    # BigMHC's allele panel is human HLA-I; mouse presets are not offered.
    return jsonify({k: v["alleles"] for k, v in allele_lib.PRESETS.items() if v["human"]})


@app.post("/api/run")
def api_run():
    d = request.get_json(force=True)
    alleles = [a for a in d.get("alleles", []) if a.strip()]
    models = [m for m in d.get("models", ["el", "im"]) if m in ("el", "im")]
    if not alleles or not models:
        return jsonify(error="Select at least one allele and one model."), 400
    lengths = sorted({int(x) for x in d.get("lengths", [9]) if 8 <= int(x) <= 14})
    text = d.get("text", "")
    if d.get("mode", "scan") == "peptides":
        peps = bm.given_peptides([t for t in text.replace(">", " ").split()
                                  if 8 <= len(t) <= 14 and not set(t.upper()) - bm.AA])
    else:
        peps = bm.windows(bm.parse_fasta_text(text), lengths)
    if peps.empty:
        return jsonify(error="No valid peptides."), 400
    try:
        grid = bm.predict(peps, alleles, models)
    except (bm.BigMHCError, config.ConfigError) as e:
        return jsonify(error=str(e)), 500
    grid = grid.astype(object).where(grid.notna(), None)
    return jsonify(rows=grid.to_dict("records"), columns=list(grid.columns),
                   n_peptides=int(peps["peptide"].nunique()))


PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>BigMHC</title>
<style>
body{font:14px system-ui,sans-serif;margin:0;display:flex;height:100vh}
#side{width:400px;padding:14px;overflow:auto;background:#f4f6f8;border-right:1px solid #ccd}
#main{flex:1;padding:14px;overflow:auto}
textarea{width:100%;font-family:monospace} h3{margin:14px 0 6px}
label.i{display:block;margin:3px 0} .presets button{margin:2px}
table{border-collapse:collapse;font-size:13px} th,td{border:1px solid #ccd;padding:2px 6px}
th{background:#eef;cursor:pointer;position:sticky;top:0}
</style></head><body>
<div id="side">
<h2>BigMHC - MHC-I EL / IM</h2>
<h3>1. Sequence / peptides</h3>
<label class="i"><input type="radio" name="mode" value="scan" checked> Scan FASTA/protein into k-mers</label>
<label class="i"><input type="radio" name="mode" value="peptides"> Use peptides as given</label>
<textarea id="text" rows="8">&gt;demo
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF</textarea>
<div>Lengths: <span id="lens"></span></div>
<h3>2. Alleles (human HLA-I)</h3>
<div class="presets" id="presets"></div>
<textarea id="alleles" rows="5" placeholder="one per line, e.g. HLA-A*02:01"></textarea>
<h3>3. Models</h3>
<label class="i"><input type="checkbox" id="m_el" checked> EL (presentation)</label>
<label class="i"><input type="checkbox" id="m_im" checked> IM (immunogenicity)</label>
<p><button id="go" onclick="run()" style="font-size:16px;padding:6px 18px">Run</button> <span id="status"></span></p>
<small>CPU runs can take minutes for large scans.</small>
</div>
<div id="main"><div id="summary"></div><button onclick="csv()">Download CSV</button><div id="tbl"></div></div>
<script>
let DATA=null,sortKey=null,sortDir=1,PRE={};const $=id=>document.getElementById(id);
$("lens").innerHTML=[8,9,10,11,12,13,14].map(n=>`<label><input type="checkbox" class="len" value="${n}" ${n==9?"checked":""}>${n}</label> `).join("");
fetch("api/presets").then(r=>r.json()).then(d=>{PRE=d;$("presets").innerHTML=Object.keys(d).map(k=>`<button onclick='$("alleles").value=PRE[${JSON.stringify(k)}].join("\\n")'>${k} (${d[k].length})</button>`).join("");});
async function run(){
 $("go").disabled=true;$("status").textContent="running...";
 const body={text:$("text").value,mode:document.querySelector("input[name=mode]:checked").value,
  lengths:[...document.querySelectorAll(".len:checked")].map(x=>+x.value),
  alleles:$("alleles").value.split(/[\s,]+/).filter(Boolean),
  models:["el","im"].filter(m=>$("m_"+m).checked)};
 const r=await fetch("api/run",{method:"POST",headers:{"Content-Type":"application/json"},body:JSON.stringify(body)});
 const d=await r.json();$("go").disabled=false;$("status").textContent="";
 if(d.error){$("summary").innerHTML="<pre style='color:#b00;white-space:pre-wrap'>"+d.error.replace(/</g,"&lt;")+"</pre>";return;}
 DATA=d;render();}
function render(){let rows=DATA.rows.slice();
 if(sortKey)rows.sort((a,b)=>(a[sortKey]>b[sortKey]?1:a[sortKey]<b[sortKey]?-1:0)*sortDir);
 $("summary").innerHTML=`<b>${DATA.n_peptides}</b> unique peptides, ${DATA.rows.length} peptide-allele rows`;
 $("tbl").innerHTML="<table><tr>"+DATA.columns.map(c=>`<th onclick="srt('${c}')">${c}</th>`).join("")+"</tr>"+
  rows.map(r=>"<tr>"+DATA.columns.map(c=>`<td>${r[c]??""}</td>`).join("")+"</tr>").join("")+"</table>";}
function srt(k){sortDir=(sortKey==k)?-sortDir:1;sortKey=k;render();}
function csv(){if(!DATA)return;const l=[DATA.columns.join(",")].concat(DATA.rows.map(r=>DATA.columns.map(c=>r[c]??"").join(",")));
 const a=document.createElement("a");a.href=URL.createObjectURL(new Blob([l.join("\n")]));a.download="bigmhc_results.csv";a.click();}
</script></body></html>"""


@app.get("/")
def index():
    return PAGE


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(os.environ.get("PORT", 5001)), threaded=True)
