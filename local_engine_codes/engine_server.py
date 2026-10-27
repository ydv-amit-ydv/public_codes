"""Epitope Filter Engine server: one page with a tab per model plus a Settings tab
whose JSON (engine_config.json) can be edited live.

    python3 engine_server.py          # http://127.0.0.1:5000

Each model keeps its own interpreter/install (see common/config.py); this process
only needs flask and pandas. Settings are re-read on every request, so saving in
the GUI takes effect immediately.
"""
import os, subprocess, sys

from flask import Flask, jsonify, request
from werkzeug.middleware.dispatcher import DispatcherMiddleware
from werkzeug.serving import run_simple

HERE = os.path.dirname(os.path.abspath(__file__))
for sub in ("MHC1/BA_EL_scripts", "MHC1/Immunogenicity_scripts"):
    sys.path.insert(0, os.path.join(HERE, sub))
sys.path.insert(0, HERE)

from common import config                                   # noqa: E402
import netmhcpan_app, bigmhc_app, bigmhc_runner as bm       # noqa: E402

main = Flask(__name__)


def _check(name, ok, detail=""):
    return {"name": name, "ok": bool(ok), "detail": str(detail)[:600]}


def _auto(m, key):
    return "  (auto-detected)" if key in m.get("auto", []) else ""


def health_netmhcpan(deep):
    m = config.model("netmhcpan")
    script = netmhcpan_app.netmhc_path()
    checks = [_check("install dir", os.path.isdir(m.get("home") or ""),
                     (m.get("home") or f"not found under {config.search_root()}") + _auto(m, "home"))]
    checks.append(_check("netMHCpan script executable", os.access(script, os.X_OK), script))
    if not checks[-1]["ok"]:
        return checks
    netmhcpan_app._list_alleles.cache_clear()
    n = len(netmhcpan_app.list_alleles())
    checks.append(_check("-listMHC", n > 0, f"{n} alleles" if n else "no output: is tcsh installed and NMHOME set in the script?"))
    if deep and n:
        res, _ = netmhcpan_app.run_netmhcpan(["GILGFVFTL"], ["HLA-A02:01"],
                                             dict(rankS=0.5, rankW=2.0, ba=True, t=None, raw=""))
        hit = res[res["Peptide"] == "GILGFVFTL"] if not res.empty else res
        got = (not hit.empty) and "SB" in str(hit["NMP_BindLevel"].iloc[0])
        checks.append(_check("GILGFVFTL / A02:01 -> strong binder", got,
                             "ok" if got else "unexpected or empty result"))
    return checks


def health_bigmhc(deep):
    m = config.model("bigmhc")
    py, home = m.get("python", ""), m.get("home", "")
    checks = [_check("python (has torch)" if "python" in m["auto"] else "python", os.path.isfile(py or ""),
                     (py or "no env with torch found; create the bigmhc env or set models.bigmhc.python") + _auto(m, "python")),
              _check("predict.py", os.path.isfile(os.path.join(home, "src", "predict.py")),
                     (home or f"BigMHC clone not found under {config.search_root()}") + _auto(m, "home"))]
    if not all(c["ok"] for c in checks):
        return checks
    r = subprocess.run([py, "-c", "import torch, numpy, pandas, psutil; print('torch', torch.__version__)"],
                       capture_output=True, text=True)
    checks.append(_check("imports (torch, numpy, pandas, psutil)", r.returncode == 0,
                         r.stdout.strip() or r.stderr.strip().splitlines()[-1:]))
    if deep and checks[-1]["ok"]:
        try:
            g = bm.predict(bm.given_peptides(["GILGFVFTL"]), ["HLA-A*02:01"], ["el", "im"])
            checks.append(_check("GILGFVFTL / A*02:01 prediction", True,
                                 f"EL {g['BigMHC_EL'].iloc[0]:.3f}  IM {g['BigMHC_IM'].iloc[0]:.3f}"))
        except Exception as e:
            checks.append(_check("GILGFVFTL / A*02:01 prediction", False, e))
    return checks


HEALTH = {"netmhcpan": health_netmhcpan, "bigmhc": health_bigmhc}


@main.get("/api/config")
def get_config():
    try:
        config.load()
        err = None
    except config.ConfigError as e:
        err = str(e)
    eff = None
    if not err:
        eff = config.effective()
    return jsonify(path=config.CONFIG_PATH, text=config.raw_text(), error=err,
                   models=list(HEALTH), effective=eff)


@main.put("/api/config")
def put_config():
    try:
        config.save_text(request.get_json(force=True).get("text", ""))
    except config.ConfigError as e:
        return jsonify(error=str(e)), 400
    return jsonify(ok=True, path=config.CONFIG_PATH)


@main.post("/api/health/<name>")
def health(name):
    if name not in HEALTH:
        return jsonify(error="unknown model"), 404
    try:
        cfg = config.model(name)
        if not cfg.get("enabled", True):
            return jsonify(ok=False, checks=[_check("enabled", False, "disabled in settings")])
        checks = HEALTH[name](request.args.get("deep") == "1")
    except config.ConfigError as e:
        checks = [_check("settings", False, e)]
    except Exception as e:                                  # health must never 500 the page
        checks = [_check("unexpected error", False, e)]
    return jsonify(ok=all(c["ok"] for c in checks), checks=checks)


PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>Epitope Filter Engine</title>
<style>
body{font:14px system-ui,sans-serif;margin:0;height:100vh;display:flex;flex-direction:column}
nav{background:#223;color:#fff;padding:6px 10px;display:flex;gap:6px;align-items:center}
nav b{margin-right:12px} nav button{background:#445;color:#fff;border:0;padding:6px 14px;border-radius:4px;cursor:pointer}
nav button.on{background:#69f} .dot{display:inline-block;width:9px;height:9px;border-radius:50%;background:#888;margin-right:5px}
.dot.ok{background:#3c3}.dot.bad{background:#e44}
iframe{flex:1;border:0;width:100%} .pane{flex:1;display:none;min-height:0} .pane.on{display:flex}
#settings{padding:14px;flex-direction:column;gap:8px;overflow:auto}
#cfg{width:100%;height:300px;font-family:monospace;font-size:13px}
.card{border:1px solid #ccd;border-radius:6px;padding:8px;margin:6px 0;background:#f8f9fb}
.card .c{margin-left:14px} .ok{color:#287} .bad{color:#b00}
</style></head><body>
<nav><b>Epitope Filter Engine</b>
<button data-t="netmhcpan"><span class="dot" id="d_netmhcpan"></span>NetMHCpan BA/EL</button>
<button data-t="bigmhc"><span class="dot" id="d_bigmhc"></span>BigMHC EL/IM</button>
<button data-t="settings">Settings</button></nav>
<div class="pane" id="p_netmhcpan"><iframe src="netmhcpan/"></iframe></div>
<div class="pane" id="p_bigmhc"><iframe src="bigmhc/"></iframe></div>
<div class="pane" id="settings">
 <div>Config file: <code id="path"></code> (optional; <code>"auto"</code> = find it beside this folder). Saved changes apply immediately, no restart.</div>
 <textarea id="cfg" spellcheck="false"></textarea>
 <div><button onclick="save()">Save</button> <button onclick="loadCfg()">Revert</button> <span id="msg"></span></div>
 <h3>In use right now</h3>
 <pre id="eff" style="background:#f4f6f8;padding:8px;margin:0"></pre>
 <h3>Model checks</h3>
 <div id="cards"></div>
</div>
<script>
const $=id=>document.getElementById(id), MODELS=["netmhcpan","bigmhc"];
function show(t){document.querySelectorAll("nav button").forEach(b=>b.classList.toggle("on",b.dataset.t==t));
 document.querySelectorAll(".pane").forEach(p=>p.classList.toggle("on",p.id=="p_"+t||p.id==t));}
document.querySelectorAll("nav button").forEach(b=>b.onclick=()=>show(b.dataset.t));
async function loadCfg(){const d=await (await fetch("api/config")).json();
 $("path").textContent=d.path;$("cfg").value=d.text;
 $("eff").textContent=d.effective?JSON.stringify(d.effective,null,2):"(fix the config to see this)";$("msg").innerHTML=d.error?`<span class="bad">${d.error}</span>`:"";}
async function save(){const r=await fetch("api/config",{method:"PUT",headers:{"Content-Type":"application/json"},body:JSON.stringify({text:$("cfg").value})});
 const d=await r.json();$("msg").innerHTML=r.ok?'<span class="ok">saved</span>':`<span class="bad">${d.error}</span>`;
 if(r.ok){loadCfg();checkAll(false);}}
async function check(m,deep){const d=await (await fetch(`api/health/${m}`+(deep?"?deep=1":""),{method:"POST"})).json();
 $("d_"+m).className="dot "+(d.ok?"ok":"bad");return d;}
async function checkAll(deep){$("cards").innerHTML="";
 for(const m of MODELS){const d=await check(m,deep);
  $("cards").innerHTML+=`<div class="card"><b>${m}</b> <span class="${d.ok?"ok":"bad"}">${d.ok?"ready":"not ready"}</span>
   <button onclick="deepTest('${m}')">Test run</button>`+d.checks.map(c=>`<div class="c ${c.ok?"ok":"bad"}">${c.ok?"&#10003;":"&#10007;"} ${c.name}${c.detail?" - "+c.detail.replace(/</g,"&lt;"):""}</div>`).join("")+"</div>";}}
async function deepTest(m){$("msg").textContent="testing "+m+"...";const d=await check(m,true);
 $("msg").textContent="";await checkAll(false);alert(m+(d.ok?" test run OK":" test run FAILED")+"\n"+d.checks.map(c=>(c.ok?"ok  ":"FAIL ")+c.name+" "+c.detail).join("\n"));}
loadCfg();checkAll(false);show("settings");
</script></body></html>"""


@main.get("/")
def index():
    return PAGE


application = DispatcherMiddleware(main, {"/netmhcpan": netmhcpan_app.app, "/bigmhc": bigmhc_app.app})

if __name__ == "__main__":
    run_simple("127.0.0.1", int(os.environ.get("PORT", 5000)), application, threaded=True)
