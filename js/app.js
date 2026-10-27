/* app.js — interactive dashboard for the Exoplanet Transit Simulator. */
(function () {
  "use strict";
  const P = window.TransitPhysics;
  const COLORS = ["#2457d6", "#0e8a7e", "#b03a8c", "#b7791f"];
  const MAX_PLANETS = 4;
  const MAX_POINTS = 60000;
  const $ = (id) => document.getElementById(id);

  // ------------------------------------------------------------ state
  const state = {
    star: { mass: 1, radius: 1, teff: 5772, u1: 0.4, u2: 0.26 },
    ld: true,
    planets: [],
    sel: 0,
    obs: { window: 30, cadence: 2, noise: false, sigma: 200 },
    view: { unit: "auto", scale: "transits" },
    t: 0,
    playing: false,
    note: "",
  };
  let derived = null;      // per-planet derived quantities
  let lc = null;           // light-curve arrays
  let skyStarImage = null; // cached limb-darkened star

  const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();

  // ------------------------------------------------------- formatting
  const fmt = {
    n: (x, d = 2) => (Number.isFinite(x) ? x.toFixed(d) : "—"),
    sig: (x, s = 3) => (Number.isFinite(x) ? Number(x.toPrecision(s)).toLocaleString("en") : "—"),
    ppm: (x) => (x >= 0.01 ? `${(x * 100).toFixed(2)} %` : `${Math.round(x * 1e6).toLocaleString("en")} ppm`),
    hours: (d) => (d * 24 >= 1 ? `${(d * 24).toFixed(2)} h` : `${Math.round(d * 1440)} min`),
  };

  // ---------------------------------------------------- derivations
  function derive() {
    const s = state.star;
    const ld = { enabled: state.ld, u1: s.u1, u2: s.u2 };
    const lum = P.luminosity(s.radius, s.teff);
    const hz = P.habitableZone(s.teff, lum);
    const planets = state.planets.map((pl) => {
      const a = P.semiMajorAxis(pl.period, s.mass);
      const model = {
        period: pl.period, t0: pl.t0, inc: pl.inc, ecc: pl.ecc, omega: pl.omega,
        aRs: a / (s.radius * P.R_SUN_AU), p: (pl.radius * P.R_EARTH_RSUN) / s.radius,
      };
      const tp = P.transitParameters(model, ld);
      const teq = P.equilibriumTemperature(s.teff, s.radius, a, pl.albedo);
      return {
        ...pl, a, model, tp, teq,
        insolation: lum / (a * a),
        cls: P.planetClass(pl.radius, teq),
        hz: P.hzStatus(a, hz),
        apo: a * (1 + pl.ecc),
      };
    });
    derived = { lum, hz, ld, planets, spec: P.spectralType(s.teff) };
  }

  function computeLightCurve() {
    const W = state.obs.window;
    let cadence = state.obs.cadence;
    let n = Math.floor((W * 1440) / cadence) + 1;
    if (n > MAX_POINTS) { cadence = (W * 1440) / (MAX_POINTS - 1); n = MAX_POINTS; }
    const ts = new Float64Array(n);
    for (let i = 0; i < n; i++) ts[i] = (i * cadence) / 1440;
    const models = derived.planets.map((d) => d.model);
    const flux = P.lightCurve(ts, models, derived.ld);
    let data = null;
    if (state.obs.noise) {
      const noise = P.gaussianNoise(n, state.obs.sigma * 1e-6, 42);
      data = flux.map((f, i) => f + noise[i]);
    }
    lc = { ts, flux, data, cadence, capped: cadence !== state.obs.cadence };
  }

  // ----------------------------------------------------- controls UI
  function slider(container, o) {
    const wrap = document.createElement("div");
    wrap.className = "ctl";
    const id = `ctl-${o.key}`;
    const toPos = (v) => (o.log ? Math.log10(v) : v);
    const fromPos = (x) => (o.log ? Math.pow(10, x) : x);
    const hardMin = o.hardMin ?? o.min, hardMax = o.hardMax ?? o.max;
    wrap.innerHTML = `
      <div class="ctl-top"><label for="${id}">${o.label}</label>
        <span class="ctl-value"><input type="number" class="num" id="${id}-num" step="any"
          min="${hardMin}" max="${hardMax}" aria-label="${o.label}"><small>${o.unit || ""}</small></span></div>
      <input type="range" id="${id}" min="${toPos(o.min)}" max="${toPos(o.max)}" step="${o.log ? 0.001 : o.step}" aria-label="${o.label} slider">`;
    const input = wrap.querySelector('input[type="range"]');
    const num = wrap.querySelector(".num");
    const show = (v) => {
      if (document.activeElement !== num) num.value = Number(v.toPrecision(o.digits || 5));
      num.classList.remove("invalid");
    };
    // let the slider reach values typed beyond its normal range
    const fitRange = (v) => {
      const pos = toPos(v);
      if (pos < parseFloat(input.min)) input.min = pos;
      if (pos > parseFloat(input.max)) input.max = pos;
    };
    fitRange(o.get());
    input.value = toPos(o.get());
    show(o.get());
    input.addEventListener("input", () => {
      let v = fromPos(parseFloat(input.value));
      if (o.round) v = o.round(v);
      o.set(v);
      show(v);
      scheduleUpdate(o.light);
    });
    const commit = () => {
      const v = parseFloat(num.value);
      if (!Number.isFinite(v) || v < hardMin || v > hardMax) {
        num.classList.add("invalid");
        num.title = `Enter a value between ${hardMin} and ${hardMax}`;
        return;
      }
      num.classList.remove("invalid");
      o.set(o.round ? o.round(v) : v);
      fitRange(o.get());
      input.value = toPos(o.get());
      if (o.after) o.after();
      scheduleUpdate(o.light);
    };
    num.addEventListener("change", commit);
    num.addEventListener("keydown", (e) => { if (e.key === "Enter") { commit(); num.blur(); } });
    container.appendChild(wrap);
    return { refresh: () => { fitRange(o.get()); input.value = toPos(o.get()); show(o.get()); } };
  }

  let controls = [];
  function buildControls() {
    controls = [];
    const s = state.star;
    const sc = $("starControls"); sc.innerHTML = "";
    controls.push(
      slider(sc, { key: "teff", label: "Temperature", min: 2300, max: 10000, step: 10, unit: "K", hardMin: 1000, hardMax: 50000,
        get: () => s.teff, set: (v) => (s.teff = v), format: (v) => Math.round(v).toLocaleString("en"), round: Math.round }),
      slider(sc, { key: "srad", label: "Radius", min: 0.08, max: 3, log: true, unit: "R☉", hardMin: 0.01, hardMax: 100,
        get: () => s.radius, set: (v) => (s.radius = v), format: (v) => v.toFixed(3) }),
      slider(sc, { key: "smass", label: "Mass", min: 0.08, max: 3, log: true, unit: "M☉", hardMin: 0.01, hardMax: 100,
        get: () => s.mass, set: (v) => (s.mass = v), format: (v) => v.toFixed(3) }),
    );
    const ldc = $("ldControls"); ldc.innerHTML = "";
    controls.push(
      slider(ldc, { key: "u1", label: "Linear coefficient u₁", min: 0, max: 1, step: 0.01, hardMin: -1, hardMax: 2,
        get: () => s.u1, set: (v) => (s.u1 = v), format: (v) => v.toFixed(2) }),
      slider(ldc, { key: "u2", label: "Quadratic coefficient u₂", min: -0.3, max: 0.7, step: 0.01, hardMin: -1, hardMax: 2,
        get: () => s.u2, set: (v) => (s.u2 = v), format: (v) => v.toFixed(2) }),
    );
    const oc = $("obsControls"); oc.innerHTML = "";
    controls.push(
      slider(oc, { key: "window", label: "Time span", min: 1, max: 1500, log: true, unit: "days", hardMin: 0.1, hardMax: 5000,
        get: () => state.obs.window, set: (v) => (state.obs.window = v), format: (v) => fmt.sig(v, 3) }),
      slider(oc, { key: "cadence", label: "Cadence", min: 0.33, max: 30, log: true, unit: "min", hardMin: 0.05, hardMax: 120,
        get: () => state.obs.cadence, set: (v) => (state.obs.cadence = v), format: (v) => fmt.sig(v, 2) }),
    );
    const nc = $("noiseControls"); nc.innerHTML = "";
    controls.push(
      slider(nc, { key: "sigma", label: "Noise per cadence", min: 10, max: 5000, log: true, unit: "ppm", hardMin: 1, hardMax: 100000,
        get: () => state.obs.sigma, set: (v) => (state.obs.sigma = v), format: (v) => Math.round(v).toLocaleString("en") }),
    );
    buildPlanetControls();
  }

  function buildPlanetTabs() {
    const tabs = $("planetTabs");
    tabs.innerHTML = "";
    state.planets.forEach((pl, i) => {
      const b = document.createElement("button");
      b.type = "button";
      b.className = "planet-tab";
      b.setAttribute("role", "tab");
      b.setAttribute("aria-selected", String(i === state.sel));
      b.style.color = pl.color;
      b.innerHTML = `<span class="dot" style="background:${pl.color}"></span><span style="color:var(--ink)">${escapeHtml(pl.name)}</span>`;
      b.addEventListener("click", () => { state.sel = i; buildPlanetTabs(); buildPlanetControls(); drawAll(); });
      tabs.appendChild(b);
    });
    $("addPlanet").disabled = state.planets.length >= MAX_PLANETS;
    $("removePlanet").disabled = state.planets.length <= 1;
  }

  function buildPlanetControls() {
    const pc = $("planetControls");
    pc.innerHTML = "";
    const pl = state.planets[state.sel];
    if (!pl) return;
    const name = document.createElement("label");
    name.className = "name-field";
    name.innerHTML = `Name<input type="text" value="${escapeHtml(pl.name)}" maxlength="30">`;
    name.querySelector("input").addEventListener("input", (e) => {
      pl.name = e.target.value || "Planet";
      buildPlanetTabs(); scheduleUpdate(true);
    });
    pc.appendChild(name);
    let afterPeriod = () => {};
    const pctl = [
      slider(pc, { key: "prad", label: "Radius", min: 0.3, max: 25, log: true, unit: "R⊕", hardMin: 0.05, hardMax: 300,
        get: () => pl.radius, set: (v) => (pl.radius = v), format: (v) => v.toFixed(2) }),
      slider(pc, { key: "period", label: "Orbital period", min: 0.3, max: 1000, log: true, unit: "days", hardMin: 0.05, hardMax: 100000, digits: 7,
        get: () => pl.period, set: (v) => { pl.period = v; pl.t0 = Math.min(pl.t0, v); }, format: (v) => fmt.sig(v, 4),
        after: () => afterPeriod() }),
      slider(pc, { key: "inc", label: "Inclination", min: 80, max: 90, step: 0.01, unit: "°", hardMin: 0, hardMax: 90,
        get: () => pl.inc, set: (v) => (pl.inc = v), format: (v) => v.toFixed(2) }),
      slider(pc, { key: "ecc", label: "Eccentricity", min: 0, max: 0.9, step: 0.005, hardMax: 0.99,
        get: () => pl.ecc, set: (v) => (pl.ecc = v), format: (v) => v.toFixed(3) }),
      slider(pc, { key: "omega", label: "Argument of periastron", min: 0, max: 360, step: 1, unit: "°",
        get: () => pl.omega, set: (v) => (pl.omega = v), format: (v) => Math.round(v) }),
      slider(pc, { key: "t0", label: "First mid-transit", min: 0, max: 1000, step: 0.01, unit: "days", hardMax: 100000,
        get: () => pl.t0, set: (v) => (pl.t0 = ((v % pl.period) + pl.period) % pl.period), format: (v) => fmt.sig(v, 4) }),
      slider(pc, { key: "albedo", label: "Bond albedo", min: 0, max: 0.9, step: 0.01, light: true, hardMax: 0.99,
        get: () => pl.albedo, set: (v) => (pl.albedo = v), format: (v) => v.toFixed(2) }),
    ];
    // the mid-transit slider range follows the period
    const t0Input = pc.querySelector("#ctl-t0");
    const syncT0 = () => { t0Input.max = pl.period; t0Input.step = pl.period / 500; };
    syncT0();
    afterPeriod = () => { syncT0(); pctl[5].refresh(); };
    pc.querySelector("#ctl-period").addEventListener("input", afterPeriod);
  }

  function escapeHtml(s) {
    return String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
  }

  // --------------------------------------------------------- updates
  let pending = null;
  let pendingFull = false;
  function scheduleUpdate(lightOnly) {
    if (!lightOnly) pendingFull = true;
    if (pending) return;
    pending = requestAnimationFrame(() => {
      pending = null;
      const full = pendingFull;
      pendingFull = false;
      derive();
      if (full) { computeLightCurve(); skyStarImage = null; }
      updateTimeSlider();
      drawAll(full);
    });
  }

  function updateTimeSlider() {
    const ts = $("timeSlider");
    ts.max = state.obs.window;
    ts.step = state.obs.window / 4000;
    if (state.t > state.obs.window) state.t = 0;
    ts.value = state.t;
  }

  function drawAll(full = true) {
    renderStarReadout();
    drawSky();
    drawOrbit();
    if (full) { plotLightCurve(); plotZooms(); }
    renderTable();
    renderTime();
  }

  function renderStarReadout() {
    const hz = derived.hz;
    $("starReadout").innerHTML =
      `${derived.spec}-type star, <b>${fmt.sig(derived.lum, 3)} L☉</b>. ` +
      `Habitable zone <b>${fmt.sig(hz.conservativeInner, 3)}–${fmt.sig(hz.conservativeOuter, 3)} AU</b>` +
      (hz.valid ? "." : " (extrapolated outside 2600–7200 K).");
  }

  // ------------------------------------------------------- canvases
  function setupCanvas(cv) {
    const dpr = window.devicePixelRatio || 1;
    const w = cv.clientWidth, h = cv.clientHeight;
    if (cv.width !== Math.round(w * dpr) || cv.height !== Math.round(h * dpr)) {
      cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
      skyStarImage = null;
    }
    const ctx = cv.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h, dpr };
  }

  function starColor(teff) {
    // approximate blackbody tint, deliberately kept warm and readable on white
    const stops = [[2500, [255, 150, 80]], [3500, [255, 180, 105]], [5000, [255, 214, 150]],
      [5800, [255, 228, 170]], [7500, [235, 236, 255]], [10000, [205, 218, 255]]];
    if (teff <= stops[0][0]) return stops[0][1];
    for (let i = 1; i < stops.length; i++) {
      if (teff <= stops[i][0]) {
        const [t0, c0] = stops[i - 1], [t1, c1] = stops[i];
        const k = (teff - t0) / (t1 - t0);
        return c0.map((c, j) => Math.round(c + k * (c1[j] - c)));
      }
    }
    return stops[stops.length - 1][1];
  }

  function buildStarImage(Rpx, dpr) {
    const size = Math.ceil(2 * Rpx * dpr) + 2;
    const off = document.createElement("canvas");
    off.width = off.height = size;
    const g = off.getContext("2d");
    const img = g.createImageData(size, size);
    const [cr, cg, cb] = starColor(state.star.teff);
    const u1 = state.ld ? state.star.u1 : 0, u2 = state.ld ? state.star.u2 : 0;
    const c = size / 2, R = Rpx * dpr;
    for (let y = 0; y < size; y++) {
      for (let x = 0; x < size; x++) {
        const r = Math.hypot(x + 0.5 - c, y + 0.5 - c) / R;
        if (r > 1) continue;
        const I = Math.max(0, P.intensity(r, u1, u2));
        const edge = Math.min(1, (1 - r) * R);      // anti-aliased limb
        const k = (y * size + x) * 4;
        img.data[k] = cr * (0.25 + 0.75 * I);
        img.data[k + 1] = cg * (0.18 + 0.82 * I);
        img.data[k + 2] = cb * (0.12 + 0.88 * I);
        img.data[k + 3] = 255 * edge;
      }
    }
    g.putImageData(img, 0, 0);
    return off;
  }

  function drawSky() {
    const cv = $("skyCanvas");
    const { ctx, w, h, dpr } = setupCanvas(cv);
    ctx.clearRect(0, 0, w, h);
    const R = Math.min(h * 0.36, w * 0.3);
    const cx = w / 2, cy = h / 2 + 10;
    if (!skyStarImage || skyStarImage.R !== R) {
      skyStarImage = { img: buildStarImage(R, dpr), R };
    }
    // planets behind the star first, then the star, then planets in front
    const items = derived.planets.map((d) => ({ d, pos: P.position(state.t, d.model) }));
    const drawPlanet = ({ d, pos }, front) => {
      const x = cx + pos.X * R, y = cy - pos.Y * R;
      const rp = Math.max(2.2, d.model.p * R);
      if (x < -rp || x > w + rp) return;
      ctx.beginPath();
      ctx.arc(x, y, rp, 0, 2 * Math.PI);
      ctx.fillStyle = front ? "#0d1420" : hexAlpha(d.color, 0.35);
      ctx.fill();
      ctx.lineWidth = 1.6;
      ctx.strokeStyle = d.color;
      ctx.stroke();
    };
    items.filter((it) => it.pos.Z <= 0).forEach((it) => drawPlanet(it, false));
    ctx.drawImage(skyStarImage.img, cx - R - 1 / dpr, cy - R - 1 / dpr, (skyStarImage.img.width) / dpr, (skyStarImage.img.height) / dpr);
    // transit chords
    ctx.setLineDash([4, 5]);
    ctx.lineWidth = 1;
    derived.planets.forEach((d) => {
      const yc = cy + d.tp.b * R * Math.sign(Math.cos(d.inc * Math.PI / 180) || 1);
      if (!d.tp.transits) return;
      ctx.strokeStyle = hexAlpha(d.color, 0.7);
      ctx.beginPath(); ctx.moveTo(cx - 1.25 * R, yc); ctx.lineTo(cx + 1.25 * R, yc); ctx.stroke();
    });
    ctx.setLineDash([]);
    items.filter((it) => it.pos.Z > 0).forEach((it) => drawPlanet(it, true));
    // scale note
    ctx.fillStyle = css("--muted");
    ctx.font = "12px 'IBM Plex Sans', system-ui, sans-serif";
    ctx.textAlign = "right";
    const tiny = derived.planets.some((d) => d.model.p * R < 2.2);
    ctx.fillText(tiny ? "Planets drawn at true scale, with a minimum size so they stay visible" : "Planets drawn at true scale",
      w - 14, h - 12);
    ctx.textAlign = "left";
  }

  function hexAlpha(hex, a) {
    const n = parseInt(hex.slice(1), 16);
    return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${a})`;
  }

  function drawOrbit() {
    const cv = $("orbitCanvas");
    const { ctx, w, h } = setupCanvas(cv);
    ctx.clearRect(0, 0, w, h);
    const hz = derived.hz;
    const maxApo = Math.max(...derived.planets.map((d) => d.apo));
    // fit the orbits; show the habitable zone too when it is not much larger
    const extent = Math.max(maxApo, hz.conservativeOuter < 4 * maxApo ? hz.optimisticOuter : 0) * 1.12;
    const cx = w / 2, cy = h / 2 + 8;
    const scale = (Math.min(w, h - 40) / 2) / extent;
    // habitable zone rings
    const ring = (r1, r2, color) => {
      ctx.beginPath();
      ctx.arc(cx, cy, r2 * scale, 0, 2 * Math.PI);
      ctx.arc(cx, cy, r1 * scale, 0, 2 * Math.PI, true);
      ctx.fillStyle = color; ctx.fill();
    };
    ring(hz.optimisticInner, hz.optimisticOuter, css("--hz-softer"));
    ring(hz.conservativeInner, hz.conservativeOuter, css("--hz-soft"));
    // orbits and planets
    derived.planets.forEach((d) => {
      const Rs = state.star.radius * P.R_SUN_AU;
      ctx.beginPath();
      for (let k = 0; k <= 240; k++) {
        const pos = P.position(d.t0 + (k / 240) * d.period, d.model);
        const x = cx + pos.X * Rs * scale;
        const y = cy + pos.r * Math.sin(pos.angle) * Rs * scale;
        if (k === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      }
      ctx.strokeStyle = hexAlpha(d.color, 0.55);
      ctx.lineWidth = 1.2;
      ctx.stroke();
      const pos = P.position(state.t, d.model);
      ctx.beginPath();
      ctx.arc(cx + pos.X * Rs * scale, cy + pos.r * Math.sin(pos.angle) * Rs * scale, 5, 0, 2 * Math.PI);
      ctx.fillStyle = d.color; ctx.fill();
    });
    // star
    const [r, g, b] = starColor(state.star.teff);
    ctx.beginPath();
    ctx.arc(cx, cy, Math.max(4, state.star.radius * P.R_SUN_AU * scale), 0, 2 * Math.PI);
    ctx.fillStyle = `rgb(${r},${g},${b})`; ctx.fill();
    ctx.strokeStyle = css("--star"); ctx.lineWidth = 1.5; ctx.stroke();
    // observer direction and legend
    ctx.fillStyle = css("--muted");
    ctx.font = "12px 'IBM Plex Sans', system-ui, sans-serif";
    ctx.textAlign = "left";
    ctx.fillText("↓ towards Earth", 14, h - 12);
    ctx.textAlign = "right";
    ctx.fillText(`${fmt.sig(extent, 2)} AU to the edge`, w - 14, h - 12);
    ctx.fillStyle = css("--hz");
    ctx.fillText("Habitable zone", w - 14, 22);
    ctx.textAlign = "left";
  }

  // ------------------------------------------------------------ plots
  function plotTheme() {
    return {
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      font: { family: "IBM Plex Sans, system-ui, sans-serif", color: css("--ink"), size: 12 },
      grid: css("--rule"), muted: css("--muted"),
    };
  }

  function axis(title, th, extra) {
    return Object.assign({ title: { text: title, font: { size: 12, color: th.muted } }, gridcolor: th.grid,
      zerolinecolor: th.grid, linecolor: th.grid, tickfont: { color: th.muted } }, extra || {});
  }

  // ---- units and scaling shared by the light-curve plots
  function unitInfo() {
    let u = state.view.unit;
    if (u === "auto") u = Math.max(0, ...derived.planets.map((d) => d.tp.depth)) < 0.01 ? "ppm" : "pct";
    const conv = { ppm: (f) => (f - 1) * 1e6, pct: (f) => (f - 1) * 100, flux: (f) => f }[u];
    const label = { ppm: "Flux change (ppm)", pct: "Flux change (%)", flux: "Relative flux" }[u];
    const hover = { ppm: "%{y:,.0f} ppm", pct: "%{y:.3f} %", flux: "%{y:.6f}" }[u];
    const depthText = (depth) => ({ ppm: `−${Math.round(depth * 1e6).toLocaleString("en")} ppm`,
      pct: `−${(depth * 100).toFixed(3)} %`, flux: `${(1 - depth).toFixed(5)}` }[u]);
    return { u, conv, label, hover, depthText };
  }

  /** Keep only points where the model is not flat (plus run boundaries). */
  function compressModel(ts, flux) {
    const xs = [], ys = [];
    const n = ts.length;
    for (let i = 0; i < n; i++) {
      if (i === 0 || i === n - 1 || flux[i] !== 1 || flux[i - 1] !== 1 || flux[i + 1] !== 1) {
        xs.push(ts[i]); ys.push(flux[i]);
      }
    }
    return { xs, ys };
  }

  function binData(ts, data, width) {
    const xs = [], ys = [];
    let start = ts[0], sum = 0, cnt = 0, tsum = 0;
    for (let i = 0; i < ts.length; i++) {
      if (ts[i] - start >= width) {
        if (cnt) { xs.push(tsum / cnt); ys.push(sum / cnt); }
        start = start + Math.floor((ts[i] - start) / width) * width; sum = 0; cnt = 0; tsum = 0;
      }
      sum += data[i]; tsum += ts[i]; cnt++;
    }
    if (cnt) { xs.push(tsum / cnt); ys.push(sum / cnt); }
    return { xs, ys };
  }

  function binWidth() {
    const durs = derived.planets.filter((d) => d.tp.transits).map((d) => d.tp.t14);
    let w = durs.length ? Math.min(...durs) / 6 : state.obs.window / 300;
    w = Math.max(w, (lc.cadence * 3) / 1440);
    return Math.max(w, state.obs.window / 2500);   // at most ~2500 binned points
  }

  function plotLightCurve() {
    const th = plotTheme();
    const U = unitInfo();
    const base = U.conv(1);
    const model = compressModel(lc.ts, lc.flux);
    const minFlux = model.ys.reduce((m, v) => Math.min(m, v), 1);
    const transiting = derived.planets.filter((d) => d.tp.transits);
    const traces = [];
    if (lc.data) {
      traces.push({ type: "scattergl", mode: "markers", name: "Simulated data", x: lc.ts, y: Array.from(lc.data, U.conv),
        marker: { size: 2.5, color: "rgba(91,107,128,0.35)" }, hoverinfo: "skip" });
      const bw = binWidth();
      const bin = binData(lc.ts, lc.data, bw);
      traces.push({ type: "scattergl", mode: "markers", name: `Binned data (${fmt.hours(bw)})`, x: bin.xs, y: bin.ys.map(U.conv),
        marker: { size: 4.5, color: css("--ink"), opacity: 0.55 }, hoverinfo: "skip" });
    }
    // baseline then model filled down to it: each dip is a shaded notch
    traces.push({ type: "scatter", mode: "lines", x: [0, state.obs.window], y: [base, base], line: { width: 0 },
      hoverinfo: "skip", showlegend: false });
    traces.push({ type: "scatter", mode: "lines", name: "Model", x: model.xs, y: model.ys.map(U.conv),
      line: { color: css("--accent"), width: 2.4, simplify: false }, fill: "tonexty", fillcolor: hexAlpha(css("--accent").startsWith("#") ? css("--accent") : "#2457d6", 0.18),
      hovertemplate: `t = %{x:.4f} d<br>${U.hover}<extra></extra>` });

    // a marker at the bottom of every transit: always visible, even when the
    // dip is narrower than a pixel, and hover shows that transit's depth
    const dtDay = lc.cadence / 1440;
    transiting.forEach((d) => {
      const xs = [], ys = [], txt = [];
      for (let tc = d.t0; tc <= state.obs.window; tc += d.period) {
        const i0 = Math.max(0, Math.floor((tc - d.tp.t14 / 2) / dtDay)), i1 = Math.min(lc.ts.length - 1, Math.ceil((tc + d.tp.t14 / 2) / dtDay));
        let m = 1, tm = tc;
        for (let i = i0; i <= i1; i++) if (lc.flux[i] < m) { m = lc.flux[i]; tm = lc.ts[i]; }
        if (m < 1) { xs.push(tm); ys.push(U.conv(m)); txt.push(U.depthText(1 - m)); }
      }
      traces.push({ type: "scatter", mode: "markers", name: `${d.name} transits`, x: xs, y: ys, text: txt,
        marker: { size: 7, color: d.color, line: { color: css("--paper"), width: 1.2 } },
        hovertemplate: `${escapeHtml(d.name)}<br>mid-transit %{x:.3f} d<br>depth %{text}<extra></extra>`, showlegend: false });
    });

    const shapes = [cursorShape()];
    const annotations = [];
    // each transit gets a coloured band at least ~4 px wide, so short transits
    // stay visible even in light curves spanning years
    const minBand = state.obs.window / 220;
    transiting.forEach((d, k) => {
      const halfBand = Math.max(d.tp.t14, minBand) / 2;
      for (let tc = d.t0; tc <= state.obs.window; tc += d.period) {
        shapes.push({ type: "rect", xref: "x", yref: "paper", x0: tc - halfBand, x1: tc + halfBand, y0: 0, y1: 1,
          fillcolor: hexAlpha(d.color, 0.13), line: { width: 0 }, layer: "below" });
        shapes.push({ type: "line", xref: "x", yref: "paper", x0: tc, x1: tc, y0: 0.95, y1: 1, line: { color: d.color, width: 3 } });
      }
      const yd = U.conv(1 - d.tp.depth);
      shapes.push({ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: yd, y1: yd,
        line: { color: d.color, width: 1.2, dash: "dash" } });
      annotations.push({ xref: "paper", yref: "y", x: 1 - 0.27 * k, y: yd, xanchor: "right", yanchor: "bottom",
        showarrow: false, text: `${escapeHtml(d.name)}: ${U.depthText(d.tp.depth)}`, font: { size: 11, color: d.color },
        bgcolor: css("--paper"), borderpad: 2 });
    });
    if (lc.capped) annotations.push({ xref: "paper", yref: "paper", x: 0, y: 1.02, xanchor: "left", yanchor: "bottom", showarrow: false,
      text: `Cadence increased to ${fmt.sig(lc.cadence, 2)} min to keep the plot fast`, font: { size: 11, color: th.muted } });

    // y range: fit the transits (default) so the depth fills the plot
    let yaxisExtra = {};
    if (state.view.scale === "transits" && minFlux < 1) {
      const lo = U.conv(minFlux), span = base - lo;
      yaxisExtra = { range: [lo - 0.18 * span, base + 0.22 * span], autorange: false };
    } else {
      yaxisExtra = { autorange: true };
    }
    const layout = {
      ...th, margin: { l: 70, r: 16, t: 24, b: 44 }, hovermode: "closest",
      showlegend: !!lc.data, legend: { orientation: "h", x: 0, y: 1.02, yanchor: "bottom", font: { size: 11 } },
      xaxis: axis("Time (days)", th, { range: [0, state.obs.window] }),
      yaxis: axis(U.label, th, yaxisExtra),
      shapes, annotations,
    };
    Plotly.react("lcPlot", traces, layout, { responsive: true, displaylogo: false, modeBarButtonsToRemove: ["lasso2d", "select2d"] });
    const el = $("lcPlot");
    if (!el._clickBound) {
      // Jump in time when the plot area is clicked anywhere (not only on a data
      // point). A drag (zoom box) is ignored by comparing press and release.
      let down = null;
      el.addEventListener("mousedown", (e) => { down = [e.clientX, e.clientY]; });
      el.addEventListener("click", (e) => {
        const drag = el.querySelector(".nsewdrag");
        if (!drag || !down || Math.hypot(e.clientX - down[0], e.clientY - down[1]) > 4) return;
        const r = drag.getBoundingClientRect();
        if (e.clientX < r.left || e.clientX > r.right || e.clientY < r.top || e.clientY > r.bottom) return;
        const xa = el._fullLayout && el._fullLayout.xaxis;
        if (xa) setTime(xa.p2d(e.clientX - r.left));
      });
      el._clickBound = true;
    }
  }

  function cursorShape() {
    return { type: "line", xref: "x", yref: "paper", x0: state.t, x1: state.t, y0: 0, y1: 1,
      line: { color: css("--accent"), width: 1.5, dash: "dot" }, name: "cursor" };
  }

  let lastCursorUpdate = 0;
  function moveCursor(force) {
    const now = performance.now();
    if (!force && now - lastCursorUpdate < 100) return;
    lastCursorUpdate = now;
    const el = $("lcPlot");
    if (!el.layout) return;
    const shapes = (el.layout.shapes || []).slice();
    shapes[0] = cursorShape();
    Plotly.relayout(el, { shapes });
  }

  function plotZooms() {
    const grid = $("zoomGrid");
    const th = plotTheme();
    const U = unitInfo();
    if (grid.children.length !== derived.planets.length) {
      grid.innerHTML = "";
      derived.planets.forEach(() => {
        const cell = document.createElement("div");
        cell.className = "zoom-cell";
        grid.appendChild(cell);
      });
    }
    derived.planets.forEach((d, i) => {
      const cell = grid.children[i];
      if (!d.tp.transits) {
        if (cell._plot) { Plotly.purge(cell.firstChild); cell._plot = false; }
        cell.innerHTML = `<div class="zoom-empty"><div><b style="color:${d.color}">${escapeHtml(d.name)}</b><br>No transit from Earth’s direction: impact parameter ${d.tp.b.toFixed(2)} is larger than 1 + R<sub>p</sub>/R<sub>★</sub>.<br>Raise the inclination towards 90°.</div></div>`;
        return;
      }
      if (!cell._plot) { cell.innerHTML = '<div class="plot"></div>'; cell._plot = true; }
      const half = Math.max(d.tp.t14 * 0.85, 0.02);
      const n = 600;
      const ts = Array.from({ length: n }, (_, k) => d.t0 - half + (2 * half * k) / (n - 1));
      const model = P.lightCurve(ts, [d.model], derived.ld);
      const base = U.conv(1);
      const traces = [];
      if (lc.data) {
        const xs = [], ys = [];
        for (let k = 0; k < lc.ts.length; k++) {
          const dt = ((((lc.ts[k] - d.t0) % d.period) + 1.5 * d.period) % d.period) - 0.5 * d.period;
          if (Math.abs(dt) <= half) { xs.push(dt); ys.push(lc.data[k]); }
        }
        traces.push({ type: "scattergl", mode: "markers", x: xs.map((x) => x * 24), y: ys.map(U.conv),
          marker: { size: 3, color: "rgba(91,107,128,0.35)" }, hoverinfo: "skip" });
        // phase-binned data, ~20 bins across the transit
        const order = xs.map((x, k) => k).sort((a, b) => xs[a] - xs[b]);
        const bin = binData(order.map((k) => xs[k]), order.map((k) => ys[k]), d.tp.t14 / 10);
        traces.push({ type: "scatter", mode: "markers", x: bin.xs.map((x) => x * 24), y: bin.ys.map(U.conv),
          marker: { size: 6, color: css("--ink") }, hoverinfo: "skip" });
      }
      traces.push({ type: "scatter", mode: "lines", x: ts.map((t) => (t - d.t0) * 24), y: Array.from(model, U.conv),
        line: { color: d.color, width: 2.6 },
        hovertemplate: `%{x:.2f} h<br>${U.hover}<extra></extra>` });
      const yd = U.conv(1 - d.tp.depth), span = base - yd;
      const layout = {
        ...th, margin: { l: 62, r: 10, t: 30, b: 40 }, showlegend: false,
        title: { text: `${escapeHtml(d.name)}: depth ${U.depthText(d.tp.depth)}, ${fmt.hours(d.tp.t14)}`,
          font: { size: 13, color: d.color }, x: 0.02, xanchor: "left" },
        xaxis: axis("Hours from mid-transit", th),
        yaxis: axis(U.label.replace("Flux change ", ""), th, state.view.scale === "transits"
          ? { range: [yd - 0.3 * span, base + 0.3 * span], autorange: false } : { autorange: true }),
        shapes: [
          { type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: yd, y1: yd, line: { color: d.color, width: 1.2, dash: "dash" } },
          { type: "rect", xref: "x", yref: "paper", x0: -d.tp.t14 * 12, x1: d.tp.t14 * 12, y0: 0, y1: 1,
            fillcolor: hexAlpha(d.color, 0.07), line: { width: 0 }, layer: "below" },
        ],
      };
      Plotly.react(cell.firstChild, traces, layout, { responsive: true, displaylogo: false, displayModeBar: false });
    });
  }

  // ------------------------------------------------------------ table
  function renderTable() {
    const rows = derived.planets.map((d) => {
      const snr = P.transitSNR(d.tp.depth, d.tp.t14, d.tp.t23, d.period, state.obs.window, lc.cadence, state.obs.sigma);
      const hzTag = { conservative: ["conservative", "Conservative HZ"], optimistic: ["optimistic", "Optimistic HZ"],
        "too hot": ["hot", "Too hot"], "too cold": ["cold", "Too cold"] }[d.hz];
      const nTr = Math.floor(state.obs.window / d.period);
      return `<tr>
        <td><span class="swatch" style="background:${d.color}"></span>${escapeHtml(d.name)}</td>
        <td>${d.cls}</td>
        <td>${fmt.sig(d.a, 3)}</td>
        <td>${fmt.sig(d.model.aRs, 3)}</td>
        <td>${d.tp.b.toFixed(2)}${d.tp.grazing ? " (grazing)" : ""}</td>
        <td>${d.tp.transits ? fmt.ppm(d.tp.depth) : "no transit"}</td>
        <td>${d.tp.transits ? fmt.hours(d.tp.t14) : "—"}</td>
        <td>${Math.round(d.teq)} K</td>
        <td>${fmt.sig(d.insolation, 3)}</td>
        <td><span class="tag ${hzTag[0]}">${hzTag[1]}</span></td>
        <td>${(d.tp.probability * 100).toFixed(2)} %</td>
        <td>${state.obs.noise ? (d.tp.transits ? `${fmt.sig(snr, 3)} (${nTr} transit${nTr === 1 ? "" : "s"})` : "—") : "turn on noise"}</td>
      </tr>`;
    }).join("");
    $("resultsTable").innerHTML = `<thead><tr>
      <th>Planet</th><th>Type</th><th>a (AU)</th><th>a/R★</th><th>Impact b</th><th>Depth</th><th>Duration</th>
      <th>T<sub>eq</sub></th><th>Insolation (S⊕)</th><th>Habitable zone</th><th>Transit probability</th><th>SNR</th>
      </tr></thead><tbody>${rows}</tbody>`;
  }

  // ------------------------------------------------------------- time
  function currentFlux() {
    return P.lightCurve([state.t], derived.planets.map((d) => d.model), derived.ld)[0];
  }

  function renderTime() {
    const f = currentFlux();
    const inTransit = derived.planets.filter((d) => {
      const pos = P.position(state.t, d.model);
      return pos.Z > 0 && Math.hypot(pos.X, pos.Y) < 1 + d.model.p;
    }).map((d) => d.name);
    const ppm = Math.round((f - 1) * 1e6);
    const fluxText = f < 0.99 ? f.toFixed(4) : `${ppm < 0 ? "−" : ""}${Math.abs(ppm).toLocaleString("en")} ppm`;
    $("timeReadout").innerHTML = `Time <b>${state.t.toFixed(3)} d</b>, flux <b>${fluxText}</b>`;
    $("skyCaption").textContent = inTransit.length ? `${inTransit.join(" and ")} in transit` : "The star as seen from Earth";
  }

  function setTime(t) {
    state.t = Math.max(0, Math.min(state.obs.window, t));
    $("timeSlider").value = state.t;
    drawSky(); drawOrbit(); renderTime(); moveCursor(true);
  }

  let lastFrame = null;
  function tick(now) {
    if (!state.playing) { lastFrame = null; return; }
    if (lastFrame !== null) {
      const dt = (now - lastFrame) / 1000;
      // the selected planet's transit lasts ~4 s on screen, but one orbit is at most 25 s
      const sel = derived.planets[state.sel] || derived.planets[0];
      const speed = Math.max(sel.tp.transits ? sel.tp.t14 / 4 : 0, sel.period / 25, state.obs.window / 600);
      let t = state.t + dt * speed;
      if (t > state.obs.window) t = 0;
      state.t = t;
      $("timeSlider").value = t;
      drawSky(); drawOrbit(); renderTime(); moveCursor(false);
    }
    lastFrame = now;
    requestAnimationFrame(tick);
  }

  function setPlaying(on) {
    state.playing = on;
    const b = $("playBtn");
    b.textContent = on ? "❚❚" : "▶";
    b.setAttribute("aria-pressed", String(on));
    b.setAttribute("aria-label", on ? "Pause" : "Play");
    if (on) requestAnimationFrame(tick);
  }

  // ---------------------------------------------------------- presets
  function loadPreset(id) {
    const pr = window.PRESETS.find((x) => x.id === id) || window.PRESETS[0];
    state.star = { ...pr.star };
    state.planets = pr.planets.map((p, i) => ({ albedo: 0.3, ...p, color: COLORS[i % COLORS.length] }));
    state.sel = 0;
    state.note = pr.note;
    const maxP = Math.max(...state.planets.map((p) => p.period));
    const minT0 = Math.min(...state.planets.map((p) => p.t0));
    state.obs.window = Math.min(1500, Math.max(3.2 * maxP, minT0 + maxP));
    state.obs.cadence = state.obs.window > 200 ? 30 : state.obs.window > 40 ? 10 : 2;
    state.obs.sigma = state.planets[0].radius > 5 ? 2000 : 150;
    state.t = state.planets[0].t0;
    $("presetNote").textContent = pr.note;
    $("preset").value = pr.id;
    $("presetNote").textContent = pr.note;
    buildControls();
    buildPlanetTabs();
    derive(); computeLightCurve(); skyStarImage = null;
    updateTimeSlider();
    drawAll(true);
  }

  function addPlanet() {
    if (state.planets.length >= MAX_PLANETS) return;
    const last = state.planets[state.planets.length - 1];
    const used = new Set(state.planets.map((p) => p.color));
    state.planets.push({
      name: `Planet ${state.planets.length + 1}`, radius: 2, period: +(last.period * 1.7).toFixed(3),
      inc: 89.8, ecc: 0, omega: 90, t0: +(last.t0 + last.period * 0.37).toFixed(3) % (last.period * 1.7),
      albedo: 0.3, mass: null, color: COLORS.find((c) => !used.has(c)) || COLORS[0],
    });
    state.sel = state.planets.length - 1;
    $("presetNote").textContent = "Custom system.";
    buildPlanetTabs(); buildPlanetControls(); scheduleUpdate(false);
  }

  function removePlanet() {
    if (state.planets.length <= 1) return;
    state.planets.splice(state.sel, 1);
    state.sel = Math.max(0, state.sel - 1);
    buildPlanetTabs(); buildPlanetControls(); scheduleUpdate(false);
  }

  function downloadCsv() {
    const rows = ["time_days,model_flux" + (lc.data ? ",simulated_flux" : "")];
    for (let i = 0; i < lc.ts.length; i++) {
      rows.push(`${lc.ts[i].toFixed(6)},${lc.flux[i].toFixed(8)}${lc.data ? "," + lc.data[i].toFixed(8) : ""}`);
    }
    const blob = new Blob([rows.join("\n")], { type: "text/csv" });
    const a = document.createElement("a");
    a.href = URL.createObjectURL(blob);
    a.download = "transit_light_curve.csv";
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(a.href), 1000);
  }

  // ------------------------------------------------------ custom target
  const CUSTOM_KEY = "ets-custom-targets";
  function loadSavedTargets() {
    try { return JSON.parse(localStorage.getItem(CUSTOM_KEY) || "[]"); } catch (e) { return []; }
  }
  function saveTargets(list) {
    try { localStorage.setItem(CUSTOM_KEY, JSON.stringify(list.slice(0, 12))); } catch (e) { /* storage unavailable */ }
  }

  function fillPresetSelect(selectedId) {
    const sel = $("preset");
    sel.innerHTML = "";
    const add = (parent, p) => { const o = document.createElement("option"); o.value = p.id; o.textContent = p.name; parent.appendChild(o); };
    const builtIn = window.PRESETS.filter((p) => !p.custom);
    const mine = window.PRESETS.filter((p) => p.custom);
    const g1 = document.createElement("optgroup"); g1.label = "Known systems";
    builtIn.forEach((p) => add(g1, p)); sel.appendChild(g1);
    if (mine.length) {
      const g2 = document.createElement("optgroup"); g2.label = "My targets";
      mine.forEach((p) => add(g2, p)); sel.appendChild(g2);
    }
    if (selectedId) sel.value = selectedId;
  }

  function customPlanetRow(pl, i) {
    const fs = document.createElement("fieldset");
    fs.className = "cplanet";
    const color = COLORS[i % COLORS.length];
    fs.innerHTML = `
      <legend><span class="dot" style="background:${color}"></span>Planet ${i + 1}
        ${i > 0 ? '<button type="button" class="text-btn danger remove-planet">Remove</button>' : ""}</legend>
      <label class="field wide">Name <input name="pname" type="text" maxlength="30" value="${escapeHtml(pl.name)}"></label>
      <label class="field span2">Size
        <span class="pair"><input name="size" type="number" step="any" value="${+pl.radius.toPrecision(5)}">
          <select name="sizeMode" aria-label="Size unit">
            <option value="re">R⊕</option><option value="rj">R♃</option><option value="depth">depth (ppm)</option>
          </select></span></label>
      <label class="field">Orbital period (days) <input name="period" type="number" step="any" value="${+pl.period.toPrecision(8)}"></label>
      <label class="field span2">Orbit tilt
        <span class="pair"><input name="geom" type="number" step="any" value="${+pl.inc.toPrecision(6)}">
          <select name="geomMode" aria-label="Orbit tilt as">
            <option value="inc">inclination (°)</option><option value="b">impact parameter b</option>
          </select></span></label>
      <label class="field">Eccentricity <input name="ecc" type="number" step="any" min="0" max="0.99" value="${+pl.ecc.toPrecision(4)}"></label>
      <label class="field">Arg. of periastron (°) <input name="omega" type="number" step="any" value="${+pl.omega.toPrecision(4)}"></label>
      <label class="field">First mid-transit (days) <input name="t0" type="number" step="any" value="${+pl.t0.toPrecision(6)}"></label>
      <label class="field">Bond albedo <input name="albedo" type="number" step="any" min="0" max="0.99" value="${+pl.albedo.toPrecision(3)}"></label>`;
    const rm = fs.querySelector(".remove-planet");
    if (rm) rm.addEventListener("click", () => { fs.remove(); renumberCustomPlanets(); });
    return fs;
  }

  function renumberCustomPlanets() {
    const rows = [...$("customPlanets").children];
    rows.forEach((fs, i) => {
      const legend = fs.querySelector("legend");
      legend.firstChild.nextSibling.textContent = `Planet ${i + 1}\n        `;
      legend.querySelector(".dot").style.background = COLORS[i % COLORS.length];
    });
    $("customAddPlanet").disabled = rows.length >= MAX_PLANETS;
  }

  function readRowDefaults(fs) {
    const period = fs ? parseFloat(fs.querySelector('[name="period"]').value) : 10;
    return { period: Number.isFinite(period) && period > 0 ? period : 10 };
  }

  function openCustomDialog() {
    const f = $("customForm");
    const current = window.PRESETS.find((p) => p.id === $("preset").value);
    f.name.value = current && current.custom ? current.name : "My target";
    f.teff.value = Math.round(state.star.teff);
    f.srad.value = +state.star.radius.toPrecision(5);
    f.smass.value = +state.star.mass.toPrecision(5);
    f.u1.value = +state.star.u1.toPrecision(3);
    f.u2.value = +state.star.u2.toPrecision(3);
    const box = $("customPlanets");
    box.innerHTML = "";
    state.planets.forEach((pl, i) => box.appendChild(customPlanetRow(pl, i)));
    renumberCustomPlanets();
    $("customError").textContent = "";
    f.querySelectorAll(".invalid").forEach((el) => el.classList.remove("invalid"));
    setPlaying(false);
    $("customDialog").showModal();
  }

  /** Read and validate the form; returns a preset or throws with a message. */
  function readCustomForm() {
    const f = $("customForm");
    f.querySelectorAll(".invalid").forEach((el) => el.classList.remove("invalid"));
    const num = (input, label, lo, hi) => {
      const v = parseFloat(input.value);
      if (!Number.isFinite(v) || v < lo || v > hi) {
        input.classList.add("invalid");
        input.focus();
        throw new Error(`${label} must be a number between ${lo} and ${hi}.`);
      }
      return v;
    };
    const star = {
      teff: num(f.teff, "Star temperature", 1000, 50000),
      radius: num(f.srad, "Star radius", 0.01, 100),
      mass: num(f.smass, "Star mass", 0.01, 100),
      u1: num(f.u1, "Limb darkening u₁", -1, 2),
      u2: num(f.u2, "Limb darkening u₂", -1, 2),
    };
    const planets = [...$("customPlanets").children].map((fs, i) => {
      const q = (n) => fs.querySelector(`[name="${n}"]`);
      const label = `Planet ${i + 1}`;
      const sizeMode = q("sizeMode").value;
      const size = num(q("size"), `${label} size`, sizeMode === "depth" ? 1 : 0.01, sizeMode === "depth" ? 500000 : sizeMode === "rj" ? 30 : 300);
      const radius = sizeMode === "re" ? size : sizeMode === "rj" ? size * P.R_JUP_REARTH
        : (Math.sqrt(size * 1e-6) * star.radius) / P.R_EARTH_RSUN;          // geometric depth (Rp/R*)^2
      const period = num(q("period"), `${label} period`, 0.05, 100000);
      const ecc = num(q("ecc"), `${label} eccentricity`, 0, 0.99);
      const omega = num(q("omega"), `${label} argument of periastron`, -360, 720);
      let inc;
      if (q("geomMode").value === "inc") {
        inc = num(q("geom"), `${label} inclination`, 0, 90);
      } else {
        const b = num(q("geom"), `${label} impact parameter`, 0, 1000);
        const aRs = P.semiMajorAxis(period, star.mass) / (star.radius * P.R_SUN_AU);
        const w = (omega * Math.PI) / 180;
        const cosi = (b / aRs) * ((1 + ecc * Math.sin(w)) / (1 - ecc * ecc));
        if (cosi > 1) { q("geom").classList.add("invalid"); throw new Error(`${label}: impact parameter ${b} is larger than a/R★ = ${aRs.toFixed(1)}.`); }
        inc = (Math.acos(cosi) * 180) / Math.PI;
      }
      const t0raw = num(q("t0"), `${label} mid-transit time`, -1e6, 1e6);
      return {
        name: q("pname").value.trim() || label, radius, period, inc, ecc, omega: ((omega % 360) + 360) % 360,
        t0: ((t0raw % period) + period) % period, albedo: num(q("albedo"), `${label} albedo`, 0, 0.99), mass: null,
      };
    });
    if (!planets.length) throw new Error("Add at least one planet.");
    const name = f.name.value.trim() || "My target";
    return { id: `custom-${Date.now()}`, name, custom: true, star,
      note: `${name}: custom target (${planets.length} planet${planets.length > 1 ? "s" : ""}).`, planets };
  }

  function submitCustom(e) {
    e.preventDefault();
    let preset;
    try { preset = readCustomForm(); } catch (err) { $("customError").textContent = err.message; return; }
    // replace a saved target of the same name, otherwise add it
    const saved = loadSavedTargets().filter((p) => p.name !== preset.name);
    saved.unshift(preset);
    saveTargets(saved);
    window.PRESETS = window.PRESETS.filter((p) => !(p.custom && p.name === preset.name));
    window.PRESETS.push(preset);
    fillPresetSelect(preset.id);
    $("customDialog").close();
    loadPreset(preset.id);
  }

  // ------------------------------------------------------------ theme
  function applyTheme(theme) {
    document.documentElement.setAttribute("data-theme", theme);
    const btn = $("themeToggle");
    btn.textContent = theme === "dark" ? "☀" : "☾";
    btn.setAttribute("aria-label", theme === "dark" ? "Switch to light theme" : "Switch to dark theme");
    try { localStorage.setItem("ets-theme", theme); } catch (e) { /* storage unavailable */ }
    if (derived) { skyStarImage = null; drawAll(true); }
  }

  // ------------------------------------------------------------- init
  function init() {
    const sel = $("preset");
    loadSavedTargets().forEach((p) => { if (!window.PRESETS.some((q) => q.id === p.id)) window.PRESETS.push({ ...p, custom: true }); });
    fillPresetSelect();
    sel.addEventListener("change", () => { setPlaying(false); loadPreset(sel.value); });
    $("customBtn").addEventListener("click", openCustomDialog);
    $("customCancel").addEventListener("click", () => $("customDialog").close());
    $("customForm").addEventListener("submit", submitCustom);
    $("customAddPlanet").addEventListener("click", () => {
      const box = $("customPlanets");
      if (box.children.length >= MAX_PLANETS) return;
      const last = readRowDefaults(box.lastElementChild);
      box.appendChild(customPlanetRow({ name: `Planet ${box.children.length + 1}`, radius: 2, period: +(last.period * 1.8).toPrecision(5),
        inc: 89.8, ecc: 0, omega: 90, t0: +(last.period * 0.4).toPrecision(4), albedo: 0.3 }, box.children.length));
      renumberCustomPlanets();
    });
    $("unitSelect").addEventListener("change", (e) => { state.view.unit = e.target.value; plotLightCurve(); plotZooms(); });
    $("scaleSelect").addEventListener("change", (e) => { state.view.scale = e.target.value; plotLightCurve(); plotZooms(); });
    $("ldEnabled").addEventListener("change", (e) => { state.ld = e.target.checked; $("ldControls").style.opacity = state.ld ? 1 : 0.45; scheduleUpdate(false); });
    $("noiseEnabled").addEventListener("change", (e) => { state.obs.noise = e.target.checked; $("noiseControls").style.display = state.obs.noise ? "" : "none"; scheduleUpdate(false); });
    $("noiseControls").style.display = "none";
    $("addPlanet").addEventListener("click", addPlanet);
    $("removePlanet").addEventListener("click", removePlanet);
    $("playBtn").addEventListener("click", () => setPlaying(!state.playing));
    $("timeSlider").addEventListener("input", (e) => setTime(parseFloat(e.target.value)));
    $("downloadCsv").addEventListener("click", downloadCsv);
    $("themeToggle").addEventListener("click", () =>
      applyTheme(document.documentElement.getAttribute("data-theme") === "dark" ? "light" : "dark"));
    let saved = null;
    try { saved = localStorage.getItem("ets-theme"); } catch (e) { saved = null; }
    if (saved === "dark") applyTheme("dark");
    window.addEventListener("resize", () => { skyStarImage = null; drawSky(); drawOrbit(); });
    document.addEventListener("keydown", (e) => {
      if (e.code === "Space" && !/INPUT|SELECT|TEXTAREA|BUTTON/.test(document.activeElement.tagName)) {
        e.preventDefault(); setPlaying(!state.playing);
      }
    });
    loadPreset(window.PRESETS[0].id);
  }

  window.TransitApp = { state, get derived() { return derived; }, get lc() { return lc; }, loadPreset, setTime, setPlaying };
  document.addEventListener("DOMContentLoaded", init);
})();
