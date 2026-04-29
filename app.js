(function () {
  "use strict";

  // AWG -> conductor diameter (inches), bare copper
  const AWG = {
    10: 0.1019,
    12: 0.0808,
    14: 0.0641,
    16: 0.0508,
    18: 0.0403,
    20: 0.0320,
  };

  const FT_PER_M = 3.28084;
  const QUARTER_FT = (f) => 234 / f;
  const HF_F_MIN = 1.8;
  const HF_F_MAX = 54;

  // Hall (1974) loading-coil formula, adapted for a quarter-wave vertical
  // over an assumed perfect ground plane (image theory: A_dipole/2 = H_vertical).
  //
  //   L_uH = 1e6 / (68*pi^2*f^2) * {
  //              [ln(24(234/f - B)/D) - 1] * [(1 - fB/234)^2 - 1] / (234/f - B)
  //            - [ln(24(H - B)/D) - 1]    * [(f(H-B)/234)^2 - 1] / (H - B)
  //          }
  //
  // H, B in feet; f in MHz; D in inches; L in microhenries.
  function forwardL(H, B, f, D) {
    if (!(H > 0 && B > 0 && f > 0 && D > 0)) {
      throw new RangeError("non-physical");
    }
    if (B >= H) throw new RangeError("non-physical");
    const qw = QUARTER_FT(f);
    if (H >= qw) throw new RangeError("non-physical");

    const top1 = qw - B;
    const top2 = H - B;
    const ln1 = Math.log((24 * top1) / D) - 1;
    const ln2 = Math.log((24 * top2) / D) - 1;
    const sq1 = Math.pow(1 - (f * B) / 234, 2) - 1;
    const sq2 = Math.pow((f * top2) / 234, 2) - 1;
    const coef = 1e6 / (68 * Math.PI * Math.PI * f * f);
    return coef * ((ln1 * sq1) / top1 - (ln2 * sq2) / top2);
  }

  // Generic bisection. fn must change sign across [lo, hi]. Returns NaN on failure.
  function bisect(fn, lo, hi, tol, maxIter) {
    tol = tol || 1e-6;
    maxIter = maxIter || 80;
    let flo, fhi;
    try {
      flo = fn(lo);
      fhi = fn(hi);
    } catch (e) {
      return NaN;
    }
    if (!isFinite(flo) || !isFinite(fhi)) return NaN;
    if (flo === 0) return lo;
    if (fhi === 0) return hi;
    if (flo * fhi > 0) return NaN;

    for (let i = 0; i < maxIter; i++) {
      const mid = 0.5 * (lo + hi);
      let fmid;
      try {
        fmid = fn(mid);
      } catch (e) {
        return NaN;
      }
      if (!isFinite(fmid)) return NaN;
      if (Math.abs(fmid) < tol || (hi - lo) / 2 < tol) return mid;
      if (flo * fmid < 0) {
        hi = mid;
        fhi = fmid;
      } else {
        lo = mid;
        flo = fmid;
      }
    }
    return 0.5 * (lo + hi);
  }

  // Solve for `target` ('H'|'B'|'L'|'f') given the other three values + wire
  // diameter. When solving for f, an optional band restricts the search to that
  // band's edges so the solution stays where the user expects.
  function solve(target, vals, D, band) {
    if (target === "L") {
      return forwardL(vals.H, vals.B, vals.f, D);
    }
    if (target === "B") {
      const g = (B) => forwardL(vals.H, B, vals.f, D) - vals.L;
      return bisect(g, 0.01, vals.H - 0.01);
    }
    if (target === "H") {
      const lo = Math.max(vals.B + 0.01, 1);
      const hi = QUARTER_FT(vals.f) - 0.01;
      if (hi <= lo) return NaN;
      const g = (H) => forwardL(H, vals.B, vals.f, D) - vals.L;
      return bisect(g, lo, hi);
    }
    if (target === "f") {
      const fHiByGeom = 234 / vals.H - 1e-3; // need H < 234/f
      const lo = band ? band.lo : HF_F_MIN;
      const hi = Math.min(band ? band.hi : HF_F_MAX, fHiByGeom);
      if (hi <= lo) return NaN;
      const g = (f) => forwardL(vals.H, vals.B, f, D) - vals.L;
      return bisect(g, lo, hi);
    }
    return NaN;
  }

  // ---------- SWR & bandwidth ----------
  // USA general-license HF band edges (MHz). Calculator auto-detects which band
  // the design frequency falls in; if outside any band, falls back to ±5%.
  const HAM_BANDS = [
    { name: "160m", lo: 1.8, hi: 2.0 },
    { name: "80m", lo: 3.5, hi: 4.0 },
    { name: "60m", lo: 5.33, hi: 5.405 },
    { name: "40m", lo: 7.0, hi: 7.3 },
    { name: "30m", lo: 10.1, hi: 10.15 },
    { name: "20m", lo: 14.0, hi: 14.35 },
    { name: "17m", lo: 18.068, hi: 18.168 },
    { name: "15m", lo: 21.0, hi: 21.45 },
    { name: "12m", lo: 24.89, hi: 24.99 },
    { name: "10m", lo: 28.0, hi: 29.7 },
    { name: "6m", lo: 50.0, hi: 54.0 },
  ];

  // Per-band starting geometry. Roughly ~60% of free-space quarter-wave at the
  // band midpoint with the coil at H/2. User overrides are remembered per-band
  // so switching bands and back recovers what they had.
  const BAND_DEFAULTS = {
    "160m": { H: 30, B: 15 },
    "80m": { H: 25, B: 12.5 },
    "60m": { H: 20, B: 10 },
    "40m": { H: 15, B: 7.5 },
    "30m": { H: 10, B: 5 },
    "20m": { H: 8, B: 4 },
    "17m": { H: 6, B: 3 },
    "15m": { H: 5, B: 2.5 },
    "12m": { H: 4, B: 2 },
    "10m": { H: 3.5, B: 1.75 },
    "6m": { H: 2.5, B: 1.25 },
  };
  const SWR_Z0 = 50;
  const SWR_THRESHOLD = 2.0;

  function findBand(f) {
    for (const b of HAM_BANDS) if (f >= b.lo && f <= b.hi) return b;
    return { name: "± 5%", lo: f * 0.95, hi: f * 1.05, custom: true };
  }

  // Antenna feedpoint reactance using the lumped Hall model. By construction
  // X(f_design) = 0 (the chosen L equals L_required at the design frequency).
  // At neighboring frequencies, the chosen L differs from what would be needed
  // for resonance, and the residual reactance is 2π·f·(L − L_req(f)).
  function reactance(H, B, L, f, D) {
    try {
      const Lreq = forwardL(H, B, f, D);
      return 2 * Math.PI * f * (L - Lreq);
    } catch (e) {
      return NaN;
    }
  }

  function swrFromZ(R, X, Z0) {
    Z0 = Z0 || SWR_Z0;
    if (!(R > 0)) return Infinity;
    const num = (R - Z0) * (R - Z0) + X * X;
    const den = (R + Z0) * (R + Z0) + X * X;
    if (den <= 0) return Infinity;
    const g2 = num / den;
    if (g2 >= 1) return Infinity;
    const g = Math.sqrt(Math.max(0, g2));
    return (1 + g) / (1 - g);
  }

  function sweep(H, B, L, D, R, band, N) {
    N = N || 240;
    const out = [];
    for (let i = 0; i <= N; i++) {
      const f = band.lo + ((band.hi - band.lo) * i) / N;
      const X = reactance(H, B, L, f, D);
      const s = isFinite(X) ? swrFromZ(R, X, SWR_Z0) : Infinity;
      out.push({ f: f, X: X, swr: s });
    }
    return out;
  }

  // Locate the SWR=2 crossings around the minimum-SWR sample (linear interp
  // between adjacent samples). Returns NaN bw if the curve never dips below 2.
  function bandwidthSWR2(samples) {
    if (samples.length === 0) return { lo: NaN, hi: NaN, bw: NaN, minSWR: NaN, minF: NaN };
    let mi = 0;
    for (let i = 1; i < samples.length; i++) {
      if (samples[i].swr < samples[mi].swr) mi = i;
    }
    const minSWR = samples[mi].swr;
    const minF = samples[mi].f;
    if (!isFinite(minSWR) || minSWR > SWR_THRESHOLD) {
      return { lo: NaN, hi: NaN, bw: NaN, minSWR: minSWR, minF: minF };
    }
    function cross(a, b) {
      const t = (SWR_THRESHOLD - a.swr) / (b.swr - a.swr);
      return a.f + t * (b.f - a.f);
    }
    let f_lo = samples[0].f;
    let leftClipped = samples[0].swr <= SWR_THRESHOLD;
    for (let i = mi; i > 0; i--) {
      if (samples[i].swr <= SWR_THRESHOLD && samples[i - 1].swr > SWR_THRESHOLD) {
        f_lo = cross(samples[i - 1], samples[i]);
        leftClipped = false;
        break;
      }
    }
    let f_hi = samples[samples.length - 1].f;
    let rightClipped = samples[samples.length - 1].swr <= SWR_THRESHOLD;
    for (let i = mi; i < samples.length - 1; i++) {
      if (samples[i].swr <= SWR_THRESHOLD && samples[i + 1].swr > SWR_THRESHOLD) {
        f_hi = cross(samples[i], samples[i + 1]);
        rightClipped = false;
        break;
      }
    }
    return {
      lo: f_lo,
      hi: f_hi,
      bw: f_hi - f_lo,
      minSWR: minSWR,
      minF: minF,
      leftClipped: leftClipped,
      rightClipped: rightClipped,
    };
  }

  // Build an inline SVG plot of SWR vs frequency. Returns markup string.
  function renderPlot(samples, fDesign, band, bw) {
    const W = 720, Hpx = 220;
    const m = { top: 18, right: 18, bottom: 34, left: 44 };
    const pw = W - m.left - m.right;
    const ph = Hpx - m.top - m.bottom;
    const ymax = 5, ymin = 1;
    const xs = (f) => m.left + ((f - band.lo) / (band.hi - band.lo)) * pw;
    const ys = (s) =>
      m.top + (1 - (Math.min(Math.max(s, ymin), ymax) - ymin) / (ymax - ymin)) * ph;

    let pts = "";
    for (let i = 0; i < samples.length; i++) {
      const s = samples[i];
      if (!isFinite(s.swr)) continue;
      pts += xs(s.f).toFixed(1) + "," + ys(s.swr).toFixed(1) + " ";
    }

    let bwRect = "";
    if (isFinite(bw.lo) && isFinite(bw.hi) && bw.hi > bw.lo) {
      const xLo = xs(bw.lo), xHi = xs(bw.hi);
      bwRect =
        '<rect x="' + xLo + '" y="' + m.top + '" width="' + (xHi - xLo) +
        '" height="' + ph + '" fill="rgba(255,176,0,0.13)" stroke="none"/>';
    }

    const yTickVals = [1, 1.5, 2, 3, 5];
    let yMarks = "";
    for (let i = 0; i < yTickVals.length; i++) {
      const v = yTickVals[i];
      const isThr = v === SWR_THRESHOLD;
      const stroke = isThr ? "#b08d57" : "#3a2e22";
      const dash = isThr ? "5 3" : "2 4";
      const sw = isThr ? "1.2" : "1";
      yMarks +=
        '<line x1="' + m.left + '" y1="' + ys(v) + '" x2="' + (m.left + pw) +
        '" y2="' + ys(v) + '" stroke="' + stroke + '" stroke-dasharray="' + dash +
        '" stroke-width="' + sw + '"/>';
      yMarks +=
        '<text x="' + (m.left - 7) + '" y="' + (ys(v) + 3.5) +
        '" text-anchor="end" font-size="10" fill="#a08658" font-family="ui-monospace, monospace">' +
        v + "</text>";
    }

    const xTickN = 5;
    let xMarks = "";
    for (let i = 0; i <= xTickN; i++) {
      const f = band.lo + ((band.hi - band.lo) * i) / xTickN;
      xMarks +=
        '<line x1="' + xs(f) + '" y1="' + (m.top + ph) + '" x2="' + xs(f) +
        '" y2="' + (m.top + ph + 4) + '" stroke="#7a5a32" stroke-width="1"/>';
      xMarks +=
        '<text x="' + xs(f) + '" y="' + (m.top + ph + 17) +
        '" text-anchor="middle" font-size="10" fill="#a08658" font-family="ui-monospace, monospace">' +
        f.toFixed(3) + "</text>";
    }

    let designLine = "";
    if (fDesign >= band.lo && fDesign <= band.hi) {
      const xD = xs(fDesign);
      designLine =
        '<line x1="' + xD + '" y1="' + m.top + '" x2="' + xD + '" y2="' + (m.top + ph) +
        '" stroke="#ffb000" stroke-width="1.2" opacity="0.85"/>' +
        '<text x="' + xD + '" y="' + (m.top - 5) +
        '" text-anchor="middle" font-size="10" fill="#ffb000" font-family="ui-monospace, monospace">f₀</text>';
    }

    return (
      '<svg viewBox="0 0 ' + W + " " + Hpx +
      '" preserveAspectRatio="xMidYMid meet" class="swr-plot" role="img" aria-label="SWR vs frequency, ' +
      band.name + '">' +
      '<rect x="' + m.left + '" y="' + m.top + '" width="' + pw + '" height="' + ph +
      '" fill="#1c1410" stroke="#7a5a32" stroke-width="1"/>' +
      bwRect + yMarks + xMarks + designLine +
      '<polyline points="' + pts.trim() +
      '" fill="none" stroke="#ffb000" stroke-width="1.6" stroke-linejoin="round" stroke-linecap="round" style="filter: drop-shadow(0 0 2px rgba(255,176,0,0.55));"/>' +
      '<text x="' + (m.left + pw / 2) + '" y="' + (Hpx - 4) +
      '" text-anchor="middle" font-size="10" fill="#5a5246" font-family="ui-sans-serif, system-ui, sans-serif" letter-spacing="0.18em">FREQUENCY (MHz)</text>' +
      '<text x="14" y="' + (m.top + ph / 2) +
      '" text-anchor="middle" font-size="10" fill="#5a5246" font-family="ui-sans-serif, system-ui, sans-serif" letter-spacing="0.18em" transform="rotate(-90, 14, ' +
      (m.top + ph / 2) + ')">SWR</text>' +
      "</svg>"
    );
  }

  // ---------- UI ----------
  const $ = (sel) => document.querySelector(sel);
  const $$ = (sel) => document.querySelectorAll(sel);

  const state = {
    units: "ft",
    target: "L",
    awg: 12,
    feedR: 36,
    band: null, // set in applyDefaults
    bandStore: {}, // bandName -> { H, B, L, f, target } in feet/µH/MHz
  };

  const inputs = {
    H: $("#in-H"),
    B: $("#in-B"),
    L: $("#in-L"),
    f: $("#in-f"),
  };
  const unitSpans = {
    H: $("#unit-H"),
    B: $("#unit-B"),
  };
  const errorEl = $("#error");
  const resultEl = $("#result");
  const resultUnitEl = $("#result-unit");
  const resultLabelEl = $("#result-label");
  const awgSelect = $("#awg");
  const feedREl = $("#in-R");
  const swrPlotEl = $("#swr-plot");
  const swrBwEl = $("#swr-bw");
  const swrBandEl = $("#swr-band");
  const swrMinEl = $("#swr-min");

  function readField(key) {
    const v = parseFloat(inputs[key].value);
    return isFinite(v) ? v : NaN;
  }

  // Read all four fields, converting H/B to feet for math. Returns {H,B,L,f} in
  // canonical units (feet, µH, MHz). NaN where fields are blank/invalid.
  function readVals() {
    const conv = state.units === "m" ? FT_PER_M : 1;
    return {
      H: readField("H") * conv,
      B: readField("B") * conv,
      L: readField("L"),
      f: readField("f"),
    };
  }

  function fmt(value, key) {
    if (!isFinite(value)) return "—";
    if (key === "L") return value.toFixed(2);
    if (key === "f") return value.toFixed(3);
    return value.toFixed(2);
  }

  function unitFor(key) {
    if (key === "L") return "µH";
    if (key === "f") return "MHz";
    return state.units; // H or B
  }

  function validate(vals, target) {
    const known = ["H", "B", "L", "f"].filter((k) => k !== target);
    for (const k of known) {
      if (!isFinite(vals[k]) || vals[k] <= 0) {
        return { ok: false, msg: "Enter positive values for all three known fields." };
      }
    }
    if (isFinite(vals.f) && (vals.f < HF_F_MIN || vals.f > HF_F_MAX)) {
      return { ok: false, msg: "Frequency outside supported range (1.8–54 MHz)." };
    }
    if (isFinite(vals.H) && isFinite(vals.B) && vals.B >= vals.H) {
      return { ok: false, msg: "Coil position B must be below the top of the antenna (B < H)." };
    }
    if (isFinite(vals.H) && isFinite(vals.f) && vals.H >= QUARTER_FT(vals.f)) {
      return {
        ok: false,
        msg: "Antenna already ≥ a full quarter-wave at this frequency — no loading coil needed.",
      };
    }
    return { ok: true };
  }

  function setError(msg) {
    errorEl.textContent = msg || "";
    resultEl.classList.toggle("invalid", !!msg);
  }

  function applyTargetStyling() {
    for (const key of ["H", "B", "L", "f"]) {
      const isTarget = key === state.target;
      inputs[key].readOnly = isTarget;
      inputs[key].classList.toggle("solved", isTarget);
    }
    resultLabelEl.textContent = labelFor(state.target);
    resultUnitEl.textContent = unitFor(state.target);
  }

  function labelFor(key) {
    return {
      H: "Total height",
      B: "Coil position",
      L: "Inductance",
      f: "Frequency",
    }[key];
  }

  function recalc() {
    const vals = readVals();
    const v = validate(vals, state.target);
    if (!v.ok) {
      setError(v.msg);
      resultEl.textContent = "—";
      inputs[state.target].value = "";
      clearSWR();
      return;
    }
    setError("");

    const D = AWG[state.awg];
    let answer;
    try {
      answer = solve(state.target, vals, D, state.band);
    } catch (e) {
      setError("Calculation produced a non-physical result; check inputs.");
      resultEl.textContent = "—";
      clearSWR();
      return;
    }
    if (!isFinite(answer)) {
      setError("No solution in the searched range. Try different inputs.");
      resultEl.textContent = "—";
      inputs[state.target].value = "";
      clearSWR();
      return;
    }

    // For lengths, convert internal-feet result back to user units before display.
    let displayValue = answer;
    if ((state.target === "H" || state.target === "B") && state.units === "m") {
      displayValue = answer / FT_PER_M;
    }
    inputs[state.target].value = fmt(displayValue, state.target);
    resultEl.textContent = fmt(displayValue, state.target);

    updateSWR();
  }

  function clearSWR(msg) {
    if (swrPlotEl) swrPlotEl.innerHTML = "";
    if (swrBwEl) swrBwEl.textContent = msg || "—";
    if (swrBandEl) swrBandEl.textContent = "—";
    if (swrMinEl) swrMinEl.textContent = "—";
  }

  function updateSWR() {
    if (!swrPlotEl) return;
    const v = readVals();
    const D = AWG[state.awg];
    const R = state.feedR;

    // Need all four parameters resolved + feedpoint R + a physical config.
    if (
      !(v.H > 0 && v.B > 0 && v.L > 0 && v.f > 0 && R > 0) ||
      v.B >= v.H ||
      v.H >= QUARTER_FT(v.f)
    ) {
      clearSWR();
      return;
    }

    const band = state.band || findBand(v.f);
    const samples = sweep(v.H, v.B, v.L, D, R, band);
    const bw = bandwidthSWR2(samples);

    swrPlotEl.innerHTML = renderPlot(samples, v.f, band, bw);
    const fInBand = v.f >= band.lo && v.f <= band.hi;
    swrBandEl.textContent =
      band.name +
      (band.custom ? " (no ham band match)" : fInBand ? "" : " (f outside band)");
    swrMinEl.textContent = isFinite(bw.minSWR)
      ? bw.minSWR.toFixed(2) + " @ " + bw.minF.toFixed(3) + " MHz"
      : "—";

    if (isFinite(bw.bw) && bw.bw > 0) {
      const kHz = (bw.bw * 1000).toFixed(1);
      let txt = kHz + " kHz  (" + bw.lo.toFixed(3) + " – " + bw.hi.toFixed(3) + " MHz)";
      if (bw.leftClipped || bw.rightClipped) txt += "  *extends past band edge";
      swrBwEl.textContent = txt;
    } else {
      swrBwEl.textContent = "SWR ≥ 2 across the band — adjust geometry or feedpoint R";
    }
  }

  // Capture the current parameter set in canonical units (feet, µH, MHz). NaNs
  // for blank fields are preserved so we don't pretend a missing target was set.
  function snapshotForBand() {
    const v = readVals();
    return { H: v.H, B: v.B, L: v.L, f: v.f, target: state.target };
  }

  function restoreFromSnapshot(snap) {
    const factor = state.units === "m" ? 1 / FT_PER_M : 1;
    inputs.H.value = isFinite(snap.H) ? (snap.H * factor).toFixed(2) : "";
    inputs.B.value = isFinite(snap.B) ? (snap.B * factor).toFixed(2) : "";
    inputs.L.value = isFinite(snap.L) ? snap.L.toFixed(2) : "";
    inputs.f.value = isFinite(snap.f) ? snap.f.toFixed(3) : "";
    if (snap.target) {
      state.target = snap.target;
      const r = document.querySelector(
        "input[name=solveFor][value=" + snap.target + "]"
      );
      if (r) r.checked = true;
    }
  }

  function highlightActiveBand() {
    const name = state.band ? state.band.name : null;
    document.querySelectorAll(".band-btn").forEach((b) => {
      b.classList.toggle("active", b.dataset.band === name);
    });
  }

  function selectBand(bandName) {
    const newBand = HAM_BANDS.find((b) => b.name === bandName);
    if (!newBand || newBand === state.band) return;
    if (state.band) {
      state.bandStore[state.band.name] = snapshotForBand();
    }
    state.band = newBand;
    const stored = state.bandStore[bandName];
    if (stored) {
      restoreFromSnapshot(stored);
    } else {
      // Fresh band: midpoint frequency, geometry from the per-band table,
      // solve for L so the user sees an immediate sanity-check answer.
      const d = BAND_DEFAULTS[bandName] || { H: 25, B: 12.5 };
      const fmid = 0.5 * (newBand.lo + newBand.hi);
      restoreFromSnapshot({ H: d.H, B: d.B, L: NaN, f: fmid, target: "L" });
    }
    highlightActiveBand();
    applyTargetStyling();
    recalc();
  }

  function onUnitsToggle(newUnits) {
    if (newUnits === state.units) return;
    // Convert displayed H and B in place so user doesn't lose the value.
    const factor = newUnits === "m" ? 1 / FT_PER_M : FT_PER_M;
    for (const key of ["H", "B"]) {
      const v = parseFloat(inputs[key].value);
      if (isFinite(v)) inputs[key].value = (v * factor).toFixed(2);
      unitSpans[key].textContent = newUnits;
    }
    state.units = newUnits;
    if (state.target === "H" || state.target === "B") {
      resultUnitEl.textContent = newUnits;
    }
    recalc();
  }

  function bind() {
    for (const key of ["H", "B", "L", "f"]) {
      inputs[key].addEventListener("input", recalc);
    }
    $$("input[name=solveFor]").forEach((r) =>
      r.addEventListener("change", (e) => {
        state.target = e.target.value;
        inputs[state.target].value = "";
        applyTargetStyling();
        recalc();
      })
    );
    awgSelect.addEventListener("change", (e) => {
      state.awg = parseInt(e.target.value, 10);
      recalc();
    });
    $$("input[name=units]").forEach((r) =>
      r.addEventListener("change", (e) => onUnitsToggle(e.target.value))
    );
    if (feedREl) {
      feedREl.addEventListener("input", () => {
        const v = parseFloat(feedREl.value);
        if (isFinite(v) && v > 0) state.feedR = v;
        updateSWR();
      });
    }
    document.querySelectorAll(".band-btn").forEach((btn) => {
      btn.addEventListener("click", () => selectBand(btn.dataset.band));
    });
  }

  function applyDefaults() {
    inputs.H.value = "25";
    inputs.B.value = "12.50";
    inputs.f.value = "3.750";
    inputs.L.value = "";
    awgSelect.value = "12";
    document.querySelector("input[name=solveFor][value=L]").checked = true;
    document.querySelector("input[name=units][value=ft]").checked = true;
    state.units = "ft";
    state.target = "L";
    state.awg = 12;
    state.feedR = 36;
    state.band = HAM_BANDS.find((b) => b.name === "80m");
    state.bandStore = {};
    if (feedREl) feedREl.value = "36";
    unitSpans.H.textContent = "ft";
    unitSpans.B.textContent = "ft";
    highlightActiveBand();
  }

  document.addEventListener("DOMContentLoaded", () => {
    applyDefaults();
    applyTargetStyling();
    bind();
    recalc();
  });

  // expose for console debugging
  window.CoilVert = {
    forwardL,
    bisect,
    solve,
    AWG,
    reactance,
    swrFromZ,
    findBand,
    sweep,
    bandwidthSWR2,
  };
})();
