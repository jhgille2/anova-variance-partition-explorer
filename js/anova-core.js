/* anova-core.js
 *
 * Statistical core for the ANOVA variance-partition explorer.
 *
 * Balanced multi-environment agricultural trial:
 *     g genotypes (lines) x e environments x r replicates,
 * laid out as a randomized complete block design within each environment
 * (replicates nested in environments).
 *
 * Random-effects model (all effects independent):
 *     y_ijk = mu + G_i + E_j + GE_ij + R(E)_jk + eps_ijk
 *     G_i    ~ N(0, s2g)     genotype
 *     E_j    ~ N(0, s2e)     environment
 *     GE_ij  ~ N(0, s2ge)    genotype x environment interaction
 *     R(E)_jk~ N(0, s2rep)   replicate nested in environment
 *     eps_ijk~ N(0, s2err)   residual
 *
 * ANOVA is computed in closed form (no matrix algebra needed for the
 * balanced case). Variance components are estimated by the ANOVA
 * method-of-moments from the expected mean squares. Heritability is
 * reported on an entry-mean (line-mean) basis:
 *     h^2 = s2g / ( s2g + s2ge/e + s2err/(r*e) )
 * the standard multi-environment form (Holland, Nyquist & Cervantes-Martinez,
 * 2003, Plant Breeding Reviews 22:9-112).
 *
 * Special case r = 1 (unreplicated within environments): GE and residual are
 * confounded, so MS_GE serves as the residual. Genotype is tested against
 * MS_GE and the pooled (GE + residual) variance MS_GE/e enters heritability.
 *
 * Pure functions, no DOM. Also loadable in Node for testing.
 */
"use strict";

/* ------------------------------------------------------------------ */
/* Random number generation                                            */
/* ------------------------------------------------------------------ */

function mulberry32(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function makeRandn(rng) {
  let spare = null;
  return function () {
    if (spare !== null) { const v = spare; spare = null; return v; }
    let u = 0, v = 0;
    while (u === 0) u = rng();
    v = rng();
    const m = Math.sqrt(-2.0 * Math.log(u));
    spare = m * Math.sin(2.0 * Math.PI * v);
    return m * Math.cos(2.0 * Math.PI * v);
  };
}

/* ------------------------------------------------------------------ */
/* F-distribution CDF via the regularized incomplete beta function     */
/* ------------------------------------------------------------------ */

function lgamma(x) {
  // Lanczos approximation (g = 7, 9 coefficients).
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - lgamma(1 - x);
  x -= 1;
  let a = c[0];
  for (let i = 1; i < 9; i++) a += c[i] / (x + i);
  const t = x + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

function betacf(a, b, x) {
  // Continued fraction for the incomplete beta function (Numerical Recipes).
  const MAXIT = 200, EPS = 3.0e-14, FPMIN = 1.0e-300;
  let qab = a + b, qap = a + 1.0, qam = a - 1.0;
  let c = 1.0, d = 1.0 - qab * x / qap;
  if (Math.abs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  let h = d;
  for (let m = 1; m <= MAXIT; m++) {
    const m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (Math.abs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (Math.abs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    const del = d * c;
    h *= del;
    if (Math.abs(del - 1.0) < EPS) break;
  }
  return h;
}

function betai(a, b, x) {
  // Regularized incomplete beta function I_x(a, b).
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) +
                      a * Math.log(x) + b * Math.log(1 - x));
  if (x < (a + 1) / (a + b + 2))
    return bt * betacf(a, b, x) / a;
  return 1 - bt * betacf(b, a, 1 - x) / b;
}

/* P(F(d1, d2) <= x). */
function fCdf(x, d1, d2) {
  if (!(x > 0)) return 0;
  if (!Number.isFinite(x)) return 1;
  const z = (d1 * x) / (d1 * x + d2);
  return betai(d1 / 2, d2 / 2, z);
}

/* ------------------------------------------------------------------ */
/* Trial simulation                                                    */
/* ------------------------------------------------------------------ */

/* p = { g, r, e, mu, s2g, s2e, s2ge, s2rep, s2err, seed }.
 * Returns { y: Float64Array, g, r, e } with index (i, j, k) -> (i*e + j)*r + k. */
function simulateTrial(p) {
  const g = p.g | 0, r = p.r | 0, e = p.e | 0;
  const rng = mulberry32((p.seed >>> 0) || 1);
  const randn = makeRandn(rng);
  const sd = (v) => Math.sqrt(Math.max(0, v));

  const G = new Float64Array(g), E = new Float64Array(e),
        GE = new Float64Array(g * e), R = new Float64Array(e * r);
  for (let i = 0; i < g; i++) G[i] = sd(p.s2g) * randn();
  for (let j = 0; j < e; j++) E[j] = sd(p.s2e) * randn();
  for (let m = 0; m < g * e; m++) GE[m] = sd(p.s2ge) * randn();
  for (let m = 0; m < e * r; m++) R[m] = sd(p.s2rep) * randn();

  const y = new Float64Array(g * e * r);
  const sErr = sd(p.s2err);
  for (let i = 0; i < g; i++)
    for (let j = 0; j < e; j++)
      for (let k = 0; k < r; k++)
        y[(i * e + j) * r + k] =
          p.mu + G[i] + E[j] + GE[i * e + j] + R[j * r + k] + sErr * randn();
  return { y, g, r, e };
}

/* ------------------------------------------------------------------ */
/* ANOVA for the balanced g x e x r design                             */
/* ------------------------------------------------------------------ */

/* Returns { rows, vc, h2est, r }.
 * rows: [{ source, df, ss, ms, F, p, den }] with F/p null when undefined.
 * vc:   raw (untruncated) method-of-moments variance component estimates
 *       { s2g, s2e, s2ge, s2rep, s2err } (NaN where not estimable, r = 1).
 * h2est: entry-mean heritability estimate (components truncated at 0). */
function anovaTable(y, g, r, e) {
  const n = g * e * r;
  let grand = 0;
  for (let t = 0; t < n; t++) grand += y[t];
  grand /= n;

  const gMean = new Float64Array(g), eMean = new Float64Array(e),
        geMean = new Float64Array(g * e), reMean = new Float64Array(e * r);
  for (let i = 0; i < g; i++)
    for (let j = 0; j < e; j++)
      for (let k = 0; k < r; k++) {
        const v = y[(i * e + j) * r + k];
        gMean[i] += v; eMean[j] += v;
        geMean[i * e + j] += v; reMean[j * r + k] += v;
      }
  for (let i = 0; i < g; i++) gMean[i] /= e * r;
  for (let j = 0; j < e; j++) eMean[j] /= g * r;
  for (let m = 0; m < g * e; m++) geMean[m] /= r;
  for (let m = 0; m < e * r; m++) reMean[m] /= g;

  let ssG = 0, ssE = 0, ssGE = 0, ssRE = 0, ssErr = 0;
  for (let i = 0; i < g; i++) { const d = gMean[i] - grand; ssG += d * d; }
  ssG *= r * e;
  for (let j = 0; j < e; j++) { const d = eMean[j] - grand; ssE += d * d; }
  ssE *= r * g;
  for (let i = 0; i < g; i++)
    for (let j = 0; j < e; j++) {
      const d = geMean[i * e + j] - gMean[i] - eMean[j] + grand;
      ssGE += d * d;
    }
  ssGE *= r;
  for (let j = 0; j < e; j++)
    for (let k = 0; k < r; k++) {
      const d = reMean[j * r + k] - eMean[j];
      ssRE += d * d;
    }
  ssRE *= g;
  for (let i = 0; i < g; i++)
    for (let j = 0; j < e; j++)
      for (let k = 0; k < r; k++) {
        const d = y[(i * e + j) * r + k] - geMean[i * e + j] - reMean[j * r + k] + eMean[j];
        ssErr += d * d;
      }

  const dfG = g - 1, dfE = e - 1, dfGE = (g - 1) * (e - 1),
        dfRE = e * (r - 1), dfErr = e * (r - 1) * (g - 1);
  const msG = ssG / dfG, msE = ssE / dfE, msGE = ssGE / dfGE;
  const msRE = dfRE > 0 ? ssRE / dfRE : NaN;
  const msErr = dfErr > 0 ? ssErr / dfErr : NaN;

  // Expected mean squares (random-effects model):
  //   E[MS_G]    = s2err + r*s2ge        + r*e*s2g
  //   E[MS_E]    = s2err + r*s2ge        + r*g*s2e
  //   E[MS_GE]   = s2err + r*s2ge
  //   E[MS_R(E)] = s2err + g*s2rep
  //   E[MS_err]  = s2err
  const vc = {
    s2g:   (msG - msGE) / (r * e),
    s2e:   (msE - msGE) / (r * g),
    s2ge:  dfErr > 0 ? (msGE - msErr) / r : NaN,
    s2rep: dfErr > 0 ? (msRE - msErr) / g : NaN,
    s2err: dfErr > 0 ? msErr : NaN,
  };

  const fRow = (F, df1, df2) =>
    (F === null || !Number.isFinite(F) || df2 <= 0)
      ? { F: null, p: null }
      : { F, p: 1 - fCdf(F, df1, df2) };

  const rows = [];
  const rG = fRow(msG / msGE, dfG, dfGE);
  rows.push({ source: "Genotype", df: dfG, ss: ssG, ms: msG,
              F: rG.F, p: rG.p, den: "MS_GE", vc: vc.s2g, vcLabel: "σ²_G" });
  const rE = fRow(msE / msGE, dfE, dfGE);
  rows.push({ source: "Environment", df: dfE, ss: ssE, ms: msE,
              F: rE.F, p: rE.p, den: "MS_GE", vc: vc.s2e, vcLabel: "σ²_E" });
  if (dfErr > 0) {
    const rGE = fRow(msGE / msErr, dfGE, dfErr);
    rows.push({ source: "Genotype × Environment", df: dfGE, ss: ssGE, ms: msGE,
                F: rGE.F, p: rGE.p, den: "MS_error",
                vc: vc.s2ge, vcLabel: "σ²_GE" });
    const rRE = fRow(msRE / msErr, dfRE, dfErr);
    rows.push({ source: "Rep(Environment)", df: dfRE, ss: ssRE, ms: msRE,
                F: rRE.F, p: rRE.p, den: "MS_error", vc: vc.s2rep, vcLabel: "σ²_R" });
    rows.push({ source: "Residual", df: dfErr, ss: ssErr, ms: msErr,
                F: null, p: null, den: "—", vc: vc.s2err, vcLabel: "σ²_ε" });
  } else {
    // r = 1: GE and residual are confounded; MS_GE serves as the residual.
    rows.push({ source: "Genotype × Environment (= residual, r = 1)",
                df: dfGE, ss: ssGE, ms: msGE, F: null, p: null, den: "—",
                vc: NaN, vcLabel: "σ²_GE + σ²_ε" });
  }

  // Entry-mean heritability estimate (negative components truncated at 0,
  // the standard practice for method-of-moments estimators).
  const t = (v) => Math.max(0, v);
  let h2est;
  if (dfErr > 0) {
    const vg = t(vc.s2g), vge = t(vc.s2ge), ve = t(vc.s2err);
    h2est = vg / (vg + vge / e + ve / (r * e));
  } else {
    const vg = t(vc.s2g);
    h2est = vg / (vg + msGE / e); // pooled GE + residual variance
  }
  if (!Number.isFinite(h2est)) h2est = NaN;

  return { rows, vc, h2est, r };
}

/* ------------------------------------------------------------------ */
/* Heritability from true (hidden) parameters                          */
/* ------------------------------------------------------------------ */

function trueVarianceComponents(p) {
  return { s2g: p.s2g, s2e: p.s2e, s2ge: p.s2ge, s2rep: p.s2rep, s2err: p.s2err };
}

/* Entry-mean heritability implied by the true variance components. */
function trueHeritability(p) {
  const den = p.s2g + p.s2ge / p.e + p.s2err / (p.r * p.e);
  return den > 0 ? p.s2g / den : NaN;
}

/* Theoretical h^2 as a function of replications (r = 1..8) at the
 * current number of environments, and as a function of environments
 * (e = 1..8) at the current number of replications. */
function heritabilityCurves(p) {
  const vsR = [], vsE = [];
  for (let r = 1; r <= 8; r++)
    vsR.push({ x: r, h2: p.s2g / (p.s2g + p.s2ge / p.e + p.s2err / (r * p.e)) });
  for (let e = 1; e <= 8; e++)
    vsE.push({ x: e, h2: p.s2g / (p.s2g + p.s2ge / e + p.s2err / (p.r * e)) });
  return { vsR, vsE };
}

/* Solve for the genetic variance that yields a target entry-mean h^2
 * under the current design and non-genetic variance components. */
function solveS2gForHeritability(p, h2target) {
  const h = Math.min(0.99, Math.max(0.01, h2target));
  const rest = p.s2ge / p.e + p.s2err / (p.r * p.e);
  return (h / (1 - h)) * rest;
}

/* ------------------------------------------------------------------ */
/* Monte Carlo: repeat the whole trial to show estimator behavior      */
/* ------------------------------------------------------------------ */

/* Returns { h2: Float64Array, s2g: Float64Array, mean, sd, fracTruncated }. */
function monteCarlo(p, nSim) {
  const h2 = new Float64Array(nSim), s2g = new Float64Array(nSim);
  let nTrunc = 0;
  for (let s = 0; s < nSim; s++) {
    const q = Object.assign({}, p, { seed: ((p.seed * 2654435761 + s * 40503) >>> 0) || 1 });
    const t = simulateTrial(q);
    const a = anovaTable(t.y, t.g, t.r, t.e);
    h2[s] = a.h2est;
    s2g[s] = a.vc.s2g;
    if (a.vc.s2g < 0) nTrunc++;
  }
  let mean = 0;
  for (let s = 0; s < nSim; s++) mean += h2[s];
  mean /= nSim;
  let sd = 0;
  for (let s = 0; s < nSim; s++) { const d = h2[s] - mean; sd += d * d; }
  sd = Math.sqrt(sd / (nSim - 1));
  return { h2, s2g, mean, sd, fracTruncated: nTrunc / nSim };
}

/* Node export shim (browsers use the globals directly). */
if (typeof module !== "undefined" && module.exports) {
  module.exports = {
    mulberry32, makeRandn, lgamma, betai, fCdf,
    simulateTrial, anovaTable, trueVarianceComponents,
    trueHeritability, heritabilityCurves, solveS2gForHeritability, monteCarlo,
  };
}
