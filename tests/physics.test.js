// Run with:  node --test tests/
const test = require("node:test");
const assert = require("node:assert/strict");
const P = require("../js/physics.js");

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg || ""} expected ${b} ± ${tol}, got ${a}`);

function planet(o) {
  const star = o.star || { mass: 1, radius: 1, teff: 5772 };
  const a = P.semiMajorAxis(o.period, star.mass);
  return { period: o.period, t0: o.t0 || 0, inc: o.inc ?? 90, ecc: o.ecc || 0, omega: o.omega ?? 90,
           aRs: a / (star.radius * P.R_SUN_AU), p: (o.radius * P.R_EARTH_RSUN) / star.radius, a };
}

test("Kepler's third law: Earth orbit is 1 AU", () => {
  close(P.semiMajorAxis(365.25, 1), 1, 1e-9);
});

test("Earth equilibrium temperature is 255 K (original code gave ~3700 K)", () => {
  close(P.equilibriumTemperature(5772, 1, 1, 0.3), 255, 1);
});

test("Solar habitable zone matches Kopparapu et al. 2014", () => {
  const hz = P.habitableZone(5780, 1);
  close(hz.conservativeInner, 0.95, 0.005);
  close(hz.conservativeOuter, 1.676, 0.005);
  close(hz.optimisticInner, 0.75, 0.005);
  close(hz.optimisticOuter, 1.77, 0.005);
});

test("Earth transit: 84 ppm deep and ~13 h long", () => {
  const pl = planet({ period: 365.25, radius: 1 });
  const tp = P.transitParameters(pl);
  close(tp.depth * 1e6, 83.9, 0.5);
  close(tp.t14 * 24, 13.0, 0.2);
});

test("HD 209458 b: ~1.4 % deep, ~3 h long (with radius in Earth units)", () => {
  const star = { mass: 1.148, radius: 1.203, teff: 6065 };
  const pl = planet({ star, period: 3.52474859, radius: 1.38 * P.R_JUP_REARTH, inc: 86.71 });
  const tp = P.transitParameters(pl);
  close(tp.geometricDepth, 0.0139, 0.0005);
  close(tp.t14 * 24, 3.0, 0.25);
  close(tp.b, 0.51, 0.05);
});

test("exactly one transit per orbit: no dip when the planet is behind the star", () => {
  const pl = planet({ period: 10, radius: 3, t0: 2 });
  const ts = Array.from({ length: 20001 }, (_, i) => i * 0.001);   // 20 days, two orbits
  const f = P.lightCurve(ts, [pl], { enabled: false });
  const inTransit = ts.filter((t, i) => f[i] < 1);
  assert.ok(inTransit.length > 0);
  for (const t of inTransit) {
    const ph = (((t - 2) % 10) + 10) % 10;
    assert.ok(ph < 0.5 || ph > 9.5, `dip at t = ${t} is not at a transit`);
  }
  // occultation half an orbit later shows no dip
  const occ = P.lightCurve([7], [pl], { enabled: false })[0];
  assert.equal(occ, 1);
});

test("mid-transit happens at t0 even for eccentric orbits", () => {
  for (const [ecc, omega] of [[0, 90], [0.3, 30], [0.5, 200]]) {
    const pl = planet({ period: 8, radius: 2, t0: 3, ecc, omega });
    const pos = P.position(3, pl);
    assert.ok(pos.Z > 0);
    close(pos.X, 0, 1e-9, `ecc ${ecc}`);
  }
});

test("minimum flux equals 1 - depth (uniform disk)", () => {
  const pl = planet({ period: 5, radius: 5, t0: 1 });
  const ts = Array.from({ length: 4001 }, (_, i) => 0.9 + i * 0.00005);
  const f = P.lightCurve(ts, [pl], { enabled: false });
  close(1 - Math.min(...f), pl.p * pl.p, 1e-9);
});

test("overlap fraction: limits and continuity", () => {
  assert.equal(P.overlapFraction(0.1, 2), 0);
  close(P.overlapFraction(0.1, 0.5), 0.01, 1e-12);
  close(P.overlapFraction(0.1, 0.9 + 1e-9), 0.01, 1e-6);    // continuous at second contact
  close(P.overlapFraction(0.1, 1.1 - 1e-9), 0, 1e-6);       // continuous at first contact
  close(P.overlapFraction(0.1, 1.0), 0.0049, 0.0003);       // roughly half covered at the limb
});

test("limb darkening: deeper at centre, shallower near the limb, flux conserved", () => {
  const centre = P.fluxDeficit(0.05, 0, 0.4, 0.26);
  const limb = P.fluxDeficit(0.05, 0.9, 0.4, 0.26);
  assert.ok(centre > 0.0025 && limb < 0.0025);
  // integrating intensity over the disk gives the normalisation
  let s = 0; const n = 20000;
  for (let i = 0; i < n; i++) { const r = (i + 0.5) / n; s += P.intensity(r, 0.4, 0.26) * 2 * r / n; }
  close(s, 1 - 0.4 / 3 - 0.26 / 6, 1e-6);
});

test("non-transiting inclination gives no transit", () => {
  const pl = planet({ period: 365.25, radius: 1, inc: 89 });
  assert.equal(P.transitParameters(pl).transits, false);
});

test("SNR grows with the square root of the number of transits", () => {
  const one = P.transitSNR(1e-3, 0.1, 0.09, 10, 10, 2, 1000);
  const four = P.transitSNR(1e-3, 0.1, 0.09, 10, 40, 2, 1000);
  close(four / one, 2, 1e-9);
});

test("noise is reproducible and has the requested scatter", () => {
  const a = P.gaussianNoise(20000, 1e-3, 7), b = P.gaussianNoise(20000, 1e-3, 7);
  assert.deepEqual(a, b);
  const sd = Math.sqrt(a.reduce((s, x) => s + x * x, 0) / a.length);
  close(sd, 1e-3, 3e-5);
});

test("planet classes", () => {
  assert.equal(P.planetClass(1.0, 255), "Earth-sized");
  assert.equal(P.planetClass(15, 1450), "Hot Jupiter");
  assert.equal(P.planetClass(2.4, 260), "Sub-Neptune");
});
