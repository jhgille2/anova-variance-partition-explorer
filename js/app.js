/* app.js — UI wiring and canvas rendering for the ANOVA variance-partition
 * explorer. Depends on the globals from js/anova-core.js. */
"use strict";

const $ = (id) => document.getElementById(id);
const state = { p: null, anova: null, mc: null, curves: null };
let runTimer = null, mcTimer = null;

const C_TRUE = "#0072b2", C_EST = "#e69f00", C_GRID = "#e5e7eb",
      C_INK = "#374151", C_MUT = "#9ca3af";

/* ------------------------------------------------------------------ */
/* Small utilities                                                     */
/* ------------------------------------------------------------------ */

function fmt(x, d) {
  d = (d === undefined) ? 2 : d;
  if (x === null || x === undefined || !Number.isFinite(x)) return "—";
  return x.toFixed(d);
}

function fmtP(p) {
  if (p === null || !Number.isFinite(p)) return "—";
  if (p < 0.0005) return "<0.001";
  return p.toFixed(3);
}

function fitCanvas(cv, hCss) {
  const w = Math.max(80, cv.clientWidth || 600);
  const dpr = window.devicePixelRatio || 1;
  cv.width = Math.round(w * dpr);
  cv.height = Math.round(hCss * dpr);
  cv.style.height = hCss + "px";
  const ctx = cv.getContext("2d");
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, w, hCss);
  return { ctx, w, h: hCss };
}

function axes(ctx, x0, y0, x1, y1) {
  ctx.strokeStyle = C_INK; ctx.lineWidth = 1;
  ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y0); ctx.lineTo(x1, y1); ctx.stroke();
}

/* ------------------------------------------------------------------ */
/* Parameters                                                          */
/* ------------------------------------------------------------------ */

function readParams() {
  const v = (id) => parseFloat($(id).value);
  return {
    g: Math.max(2, Math.round(v("ctl-g"))),
    r: Math.max(1, Math.round(v("ctl-r"))),
    e: Math.max(2, Math.round(v("ctl-e"))),
    mu: v("ctl-mu"), s2g: v("ctl-s2g"), s2e: v("ctl-s2e"),
    s2ge: v("ctl-s2ge"), s2rep: v("ctl-s2rep"), s2err: v("ctl-s2err"),
    sims: Math.max(10, Math.round(v("ctl-sims"))),
    seed: Math.max(1, Math.round(v("ctl-seed")) || 1),
  };
}

function updateLabels() {
  const p = readParams();
  $("val-g").textContent = p.g;
  $("val-r").textContent = p.r;
  $("val-e").textContent = p.e;
  $("val-mu").textContent = fmt(p.mu, 0);
  $("val-s2g").textContent = fmt(p.s2g, 1);
  $("val-s2e").textContent = fmt(p.s2e, 1);
  $("val-s2ge").textContent = fmt(p.s2ge, 1);
  $("val-s2rep").textContent = fmt(p.s2rep, 1);
  $("val-s2err").textContent = fmt(p.s2err, 1);
  $("val-sims").textContent = p.sims;
  const h2 = trueHeritability(p);
  $("true-h2").textContent = Number.isFinite(h2) ? h2.toFixed(3) : "—";
  $("design-note").textContent =
    "Design: " + p.g + " lines × " + p.e + " environments × " + p.r +
    " replications = " + (p.g * p.e * p.r).toLocaleString() + " plots." +
    (p.r === 1 ? " With r = 1, genotype × environment and residual variance are confounded." : "");
}

/* ------------------------------------------------------------------ */
/* Rendering                                                           */
/* ------------------------------------------------------------------ */

function renderTable() {
  const a = state.anova;
  let html = "";
  for (const row of a.rows) {
    const vcHtml = (row.vc === null || Number.isNaN(row.vc))
      ? "—"
      : (row.vc < 0
          ? '<span class="neg">' + fmt(row.vc) + "</span>"
          : fmt(row.vc));
    html += "<tr><td>" + row.source + "</td><td>" + row.df + "</td>" +
            "<td>" + fmt(row.ss, 1) + "</td><td>" + fmt(row.ms, 2) + "</td>" +
            "<td>" + (row.F === null ? "—" : fmt(row.F, 2)) + "</td>" +
            "<td>" + fmtP(row.p) + "</td><td>" + vcHtml + "</td></tr>";
  }
  $("anova-body").innerHTML = html;
}

function componentList() {
  const p = state.p, vc = state.anova.vc;
  const L = [
    { label: "σ²_G", t: p.s2g, e: vc.s2g },
    { label: "σ²_E", t: p.s2e, e: vc.s2e },
  ];
  if (state.anova.r >= 2) {
    L.push({ label: "σ²_GE", t: p.s2ge, e: vc.s2ge });
    L.push({ label: "σ²_R", t: p.s2rep, e: vc.s2rep });
    L.push({ label: "σ²_ε", t: p.s2err, e: vc.s2err });
  } else {
    L.push({ label: "σ²_GE+σ²_ε", t: p.s2ge + p.s2err, e: state.anova.rows[2].ms });
  }
  return L;
}

function renderBars() {
  const { ctx, w, h } = fitCanvas($("cv-bars"), 300);
  const comps = componentList();
  const padL = 46, padR = 12, padT = 26, padB = 34;
  const iw = w - padL - padR, ih = h - padT - padB;
  let ymax = 0;
  for (const c of comps) ymax = Math.max(ymax, c.t, Math.max(0, c.e));
  if (!(ymax > 0)) ymax = 1;
  ymax *= 1.15;
  const y = (v) => padT + ih - (Math.max(0, v) / ymax) * ih;

  // gridlines + y ticks
  ctx.font = "11px sans-serif"; ctx.fillStyle = C_MUT; ctx.textAlign = "right";
  const step = niceStep(ymax / 4);
  for (let v = 0; v <= ymax * 1.001; v += step) {
    ctx.strokeStyle = C_GRID; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(padL, y(v)); ctx.lineTo(w - padR, y(v)); ctx.stroke();
    ctx.fillText(trimNum(v), padL - 6, y(v) + 4);
  }
  axes(ctx, padL, padT + ih, w - padR, padT);

  const gw = iw / comps.length, bw = Math.min(34, gw * 0.28);
  ctx.textAlign = "center";
  comps.forEach((c, i) => {
    const cx = padL + gw * (i + 0.5);
    // true bar
    ctx.fillStyle = C_TRUE;
    ctx.fillRect(cx - bw - 2, y(c.t), bw, padT + ih - y(c.t));
    // estimated bar (truncated at 0)
    ctx.fillStyle = C_EST;
    const ev = Math.max(0, c.e);
    ctx.fillRect(cx + 2, y(ev), bw, padT + ih - y(ev));
    if (c.e < 0) { // mark truncated negatives
      ctx.fillStyle = "#b91c1c";
      ctx.fillText("−0", cx + 2 + bw / 2, padT + ih - 6);
    }
    ctx.fillStyle = C_INK;
    ctx.fillText(c.label, cx, h - 12);
  });
  // legend
  ctx.textAlign = "left";
  ctx.fillStyle = C_TRUE; ctx.fillRect(padL, 6, 12, 12);
  ctx.fillStyle = C_INK; ctx.fillText("True (hidden)", padL + 17, 16);
  ctx.fillStyle = C_EST; ctx.fillRect(padL + 110, 6, 12, 12);
  ctx.fillStyle = C_INK; ctx.fillText("ANOVA estimate", padL + 127, 16);

  const nTrunc = comps.filter((c) => c.e < 0).length;
  $("bars-cap").textContent =
    "Blue: hidden true components used to simulate the trait. Orange: method-of-moments " +
    "estimates from this trial's ANOVA." +
    (nTrunc > 0
      ? " " + nTrunc + " estimate(s) were negative and are shown truncated at zero " +
        "(untruncated values appear in the ANOVA table)."
      : " All estimates are non-negative in this trial.");
}

function drawH2Curve(cv, pts, curX, xLabel) {
  const { ctx, w, h } = fitCanvas(cv, 240);
  const padL = 40, padR = 12, padT = 12, padB = 32;
  const iw = w - padL - padR, ih = h - padT - padB;
  const xs = pts.map((p) => p.x);
  const x0 = Math.min.apply(null, xs), x1 = Math.max.apply(null, xs);
  const X = (x) => padL + ((x - x0) / Math.max(1e-9, x1 - x0)) * iw;
  const Y = (v) => padT + ih - Math.min(1, Math.max(0, v)) * ih;

  ctx.font = "11px sans-serif"; ctx.fillStyle = C_MUT;
  ctx.textAlign = "right";
  [0, 0.25, 0.5, 0.75, 1].forEach((v) => {
    ctx.strokeStyle = C_GRID;
    ctx.beginPath(); ctx.moveTo(padL, Y(v)); ctx.lineTo(w - padR, Y(v)); ctx.stroke();
    ctx.fillText(v.toFixed(2), padL - 6, Y(v) + 4);
  });
  axes(ctx, padL, padT + ih, w - padR, padT);
  ctx.textAlign = "center"; ctx.fillStyle = C_MUT;
  xs.forEach((x) => ctx.fillText(String(x), X(x), h - 14));
  ctx.fillText(xLabel, padL + iw / 2, h - 1);

  // curve
  ctx.strokeStyle = C_TRUE; ctx.lineWidth = 2;
  ctx.beginPath();
  pts.forEach((p, i) => {
    const px = X(p.x), py = Y(Number.isFinite(p.h2) ? p.h2 : 0);
    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  });
  ctx.stroke();
  // current design marker
  ctx.strokeStyle = C_EST; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
  ctx.beginPath(); ctx.moveTo(X(curX), padT); ctx.lineTo(X(curX), padT + ih); ctx.stroke();
  ctx.setLineDash([]);
  const cur = pts.find((p) => p.x === curX);
  if (cur && Number.isFinite(cur.h2)) {
    ctx.fillStyle = C_EST;
    ctx.beginPath(); ctx.arc(X(curX), Y(cur.h2), 4.5, 0, 2 * Math.PI); ctx.fill();
    ctx.fillStyle = C_INK; ctx.textAlign = "left";
    ctx.fillText("h²=" + cur.h2.toFixed(3), Math.min(X(curX) + 8, w - 70), Y(cur.h2) - 8);
  }
}

function renderCurves() {
  drawH2Curve($("cv-h2r"), state.curves.vsR, state.p.r, "replications (r)");
  drawH2Curve($("cv-h2e"), state.curves.vsE, state.p.e, "environments (e)");
}

function renderH2Stats() {
  const hTrue = trueHeritability(state.p);
  $("stat-true-h2").textContent = Number.isFinite(hTrue) ? hTrue.toFixed(3) : "—";
  $("stat-est-h2").textContent = Number.isFinite(state.anova.h2est)
    ? state.anova.h2est.toFixed(3) : "—";
}

function renderMC() {
  const mc = state.mc;
  if (!mc) return;
  $("stat-mc").textContent = Number.isFinite(mc.mean)
    ? mc.mean.toFixed(3) + " ± " + mc.sd.toFixed(3) : "—";
  $("stat-trunc").textContent = (100 * mc.fracTruncated).toFixed(1) + "%";

  const { ctx, w, h } = fitCanvas($("cv-hist"), 260);
  const padL = 40, padR = 12, padT = 14, padB = 32;
  const iw = w - padL - padR, ih = h - padT - padB;
  const NB = 24, bins = new Array(NB).fill(0);
  for (const v of mc.h2) {
    if (!Number.isFinite(v)) continue;
    bins[Math.min(NB - 1, Math.max(0, Math.floor(v * NB)))]++;
  }
  const bmax = Math.max.apply(null, bins.concat([1]));
  const X = (v) => padL + v * iw, Y = (c) => padT + ih - (c / bmax) * ih * 0.92;

  ctx.font = "11px sans-serif"; ctx.fillStyle = C_MUT; ctx.textAlign = "center";
  [0, 0.25, 0.5, 0.75, 1].forEach((v) => {
    ctx.fillText(v.toFixed(2), X(v), h - 14);
  });
  axes(ctx, padL, padT + ih, w - padR, padT);
  ctx.fillText("estimated entry-mean h²", padL + iw / 2, h - 1);

  const bw = iw / NB;
  ctx.fillStyle = "#93b8d4";
  for (let i = 0; i < NB; i++) {
    if (bins[i] === 0) continue;
    ctx.fillRect(padL + i * bw + 1, Y(bins[i]), bw - 2, padT + ih - Y(bins[i]));
  }
  const hTrue = trueHeritability(state.p);
  if (Number.isFinite(hTrue)) {
    ctx.strokeStyle = C_TRUE; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(X(hTrue), padT); ctx.lineTo(X(hTrue), padT + ih); ctx.stroke();
    ctx.fillStyle = C_TRUE; ctx.textAlign = "left";
    ctx.fillText("true h²", Math.min(X(hTrue) + 6, w - 60), padT + 12);
  }
  if (Number.isFinite(mc.mean)) {
    ctx.strokeStyle = C_EST; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(X(mc.mean), padT); ctx.lineTo(X(mc.mean), padT + ih); ctx.stroke();
    ctx.setLineDash([]);
  }
  $("hist-cap").textContent =
    "Distribution of the entry-mean ĥ² estimate across " + mc.h2.length +
    " Monte Carlo repetitions of the whole trial. Blue line: true h². Orange dashed line: " +
    "MC mean. " + (100 * mc.fracTruncated).toFixed(1) + "% of trials produced a negative " +
    "σ̂²_G (truncated at zero for ĥ²).";
}

/* Nice axis step ~ target. */
function niceStep(target) {
  if (!(target > 0)) return 1;
  const mag = Math.pow(10, Math.floor(Math.log10(target)));
  const n = target / mag;
  return (n >= 5 ? 5 : n >= 2 ? 2 : 1) * mag;
}
function trimNum(v) {
  if (v === 0) return "0";
  const s = v.toPrecision(3);
  return parseFloat(s).toString();
}

/* ------------------------------------------------------------------ */
/* Main flow                                                           */
/* ------------------------------------------------------------------ */

function renderAll() {
  renderTable();
  renderBars();
  renderH2Stats();
  renderCurves();
}

function run() {
  state.p = readParams();
  const t = simulateTrial(state.p);
  state.anova = anovaTable(t.y, t.g, t.r, t.e);
  state.curves = heritabilityCurves(state.p);
  renderAll();
  $("mc-status").textContent = "Running Monte Carlo…";
  clearTimeout(mcTimer);
  mcTimer = setTimeout(() => {
    try {
      state.mc = monteCarlo(state.p, state.p.sims);
    } catch (err) {
      state.mc = null;
    }
    renderMC();
    $("mc-status").textContent = "";
  }, 30);
}

function scheduleRun() {
  clearTimeout(runTimer);
  runTimer = setTimeout(run, 350);
}

function init() {
  updateLabels();
  const ids = ["ctl-g", "ctl-r", "ctl-e", "ctl-mu", "ctl-s2g", "ctl-s2e",
               "ctl-s2ge", "ctl-s2rep", "ctl-s2err", "ctl-sims"];
  for (const id of ids) {
    $(id).addEventListener("input", () => { updateLabels(); scheduleRun(); });
  }
  $("ctl-seed").addEventListener("change", () => { updateLabels(); run(); });
  $("btn-resim").addEventListener("click", run);
  $("btn-solve").addEventListener("click", () => {
    const p = readParams();
    const target = parseFloat($("target-h2").value);
    if (!Number.isFinite(target) || target <= 0 || target >= 1) {
      $("mc-status").textContent = "Target h² must be between 0 and 1.";
      return;
    }
    const s = solveS2gForHeritability(p, target);
    $("ctl-s2g").value = Math.min(100, Math.max(0, s));
    $("mc-status").textContent =
      s > 100 ? "σ²_G = " + s.toFixed(1) + " needed; capped at slider maximum 100."
              : "";
    updateLabels();
    run();
  });
  let rzTimer = null;
  window.addEventListener("resize", () => {
    clearTimeout(rzTimer);
    rzTimer = setTimeout(() => { renderAll(); renderMC(); }, 200);
  });
  run();
}

init();
