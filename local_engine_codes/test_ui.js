// Browser smoke test: node test_ui.js http://127.0.0.1:5000/   (needs playwright; test_engine.py runs it when node and playwright are present)
const { chromium } = require('playwright');
(async () => {
  const url = process.argv[2], errs = [], fail = m => { console.log('UI TEST FAILED: ' + m); process.exitCode = 1; };
  const b = await chromium.launch(); const p = await b.newPage({ viewport: { width: 1400, height: 900 } });
  await p.addInitScript(() => { window.__rej = []; window.addEventListener('unhandledrejection', e => window.__rej.push(String(e.reason && e.reason.stack || e.reason))); });
  p.on('pageerror', e => errs.push(String(e))); p.on('console', m => { if (m.type() == 'error') errs.push(m.text()); });
  await p.goto(url); await p.evaluate(() => localStorage.clear()); await p.reload(); await p.waitForTimeout(2500);
  const n = id => p.evaluate(i => document.getElementById(i) ? document.getElementById(i).options.length : -1, id);
  if ((await n('profile')) < 4) fail('the Performance dropdown is empty');
  if ((await n('pf_profile')) < 4) fail('the Settings profile dropdown is empty');
  if ((await n('strat')) < 3) fail('the strategy dropdown is empty');
  if ((await n('s_model')) < 10) fail('the single-model dropdown is empty');
  for (const t of ['single', 'projects', 'settings', 'pipeline']) { await p.click(`nav button[data-t=${t}]`); await p.waitForTimeout(300); }
  await p.selectOption('#fmode1', 'manual'); await p.selectOption('#ftpl1', 'recommended'); await p.waitForTimeout(500);
  if ((await p.locator('#fbuild1 .arm[data-i]').count()) < 3) fail('the manual funnel builder shows no stages');
  const prm = await p.evaluate(() => params().c1.funnel.stages.length); if (prm < 3) fail('params() lost the funnel');
  await p.click('nav button[data-t=projects]'); await p.waitForTimeout(600);
  if (!(await p.evaluate(() => document.getElementById('pj_tbl').innerText.length))) fail('the Projects table did not draw');
  const rej = await p.evaluate(() => window.__rej); if (rej.length) fail('unhandled errors: ' + rej.join(' | ').slice(0, 400));
  if (errs.length) fail('page errors: ' + errs.join(' | ').slice(0, 400));
  await b.close(); if (!process.exitCode) console.log('ui ok');
})();
