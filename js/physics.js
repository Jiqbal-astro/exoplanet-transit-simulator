/*
 * physics.js — transit and planet physics for the Exoplanet Transit Simulator.
 *
 * Pure functions with no DOM access, usable in the browser (window.TransitPhysics)
 * and in Node (module.exports) for testing.
 *
 * Conventions follow Winn (2010), "Transits and Occultations":
 *   sky-plane coordinates in units of the stellar radius,
 *     X = -r cos(w + f),  Y = -r sin(w + f) cos i,  Z = r sin(w + f) sin i,
 *   with Z > 0 towards the observer. Mid-transit happens at true anomaly
 *   f = pi/2 - w. A planet only blocks starlight when Z > 0 (in front of the star).
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.TransitPhysics = api;
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  const R_SUN_AU = 0.00465047;      // solar radius in AU
  const R_EARTH_RSUN = 0.0091577;   // Earth radius in solar radii (IAU nominal)
  const R_JUP_REARTH = 11.209;      // Jupiter equatorial radius in Earth radii
  const M_JUP_MEARTH = 317.83;
  const T_SUN = 5772;               // K (IAU 2015 nominal)
  const DAY_PER_YEAR = 365.25;

  // Kopparapu et al. (2014): S_eff = S_sun + a T + b T^2 + c T^3 + d T^4, T = Teff - 5780 K
  const HZ_COEFFS = {
    recentVenus:       [1.776, 2.136e-4, 2.533e-8, -1.332e-11, -3.097e-15],
    runawayGreenhouse: [1.107, 1.332e-4, 1.580e-8, -8.308e-12, -1.931e-15],
    maximumGreenhouse: [0.356, 6.171e-5, 1.698e-9, -3.198e-12, -5.575e-16],
    earlyMars:         [0.320, 5.547e-5, 1.526e-9, -2.874e-12, -5.011e-16],
  };

  const deg = Math.PI / 180;

  // ---------------------------------------------------------------- star
  function luminosity(radius, teff) {
    return radius * radius * Math.pow(teff / T_SUN, 4);
  }

  function spectralType(teff) {
    if (teff >= 30000) return "O";
    if (teff >= 10000) return "B";
    if (teff >= 7500) return "A";
    if (teff >= 6000) return "F";
    if (teff >= 5200) return "G";
    if (teff >= 3700) return "K";
    return "M";
  }

  function habitableZone(teff, lum) {
    const t = teff - 5780;
    const d = (c) => Math.sqrt(lum / (c[0] + c[1] * t + c[2] * t * t + c[3] * t ** 3 + c[4] * t ** 4));
    return {
      optimisticInner: d(HZ_COEFFS.recentVenus),
      conservativeInner: d(HZ_COEFFS.runawayGreenhouse),
      conservativeOuter: d(HZ_COEFFS.maximumGreenhouse),
      optimisticOuter: d(HZ_COEFFS.earlyMars),
      valid: teff >= 2600 && teff <= 7200,
    };
  }

  // ----------------------------------------------------------- orbits
  function semiMajorAxis(periodDays, starMass) {
    return Math.cbrt(starMass * Math.pow(periodDays / DAY_PER_YEAR, 2));
  }

  function solveKepler(M, e) {
    M = M % (2 * Math.PI);
    if (M < 0) M += 2 * Math.PI;
    let E = e < 0.8 ? M : Math.PI;
    for (let k = 0; k < 60; k++) {
      const dE = (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
      E -= dE;
      if (Math.abs(dE) < 1e-12) break;
    }
    return E;
  }

  function trueToMean(f, e) {
    const E = 2 * Math.atan2(Math.sqrt(1 - e) * Math.sin(f / 2), Math.sqrt(1 + e) * Math.cos(f / 2));
    return E - e * Math.sin(E);
  }

  /** Planet position at time t (days). Returns sky coords X, Y, Z in stellar radii.
   *  `pl` needs period, t0 (mid-transit time), ecc, omega (deg), inc (deg) and
   *  aRs (= a / R_star). */
  function position(t, pl) {
    const e = pl.ecc || 0;
    const w = (pl.omega || 0) * deg;
    const Mt = trueToMean(Math.PI / 2 - w, e);                 // mean anomaly at mid-transit
    const M = Mt + (2 * Math.PI * (t - pl.t0)) / pl.period;
    const E = solveKepler(M, e);
    const f = 2 * Math.atan2(Math.sqrt(1 + e) * Math.sin(E / 2), Math.sqrt(1 - e) * Math.cos(E / 2));
    const r = pl.aRs * (1 - e * Math.cos(E));
    const inc = pl.inc * deg;
    return {
      X: -r * Math.cos(w + f),
      Y: -r * Math.sin(w + f) * Math.cos(inc),
      Z: r * Math.sin(w + f) * Math.sin(inc),
      r,
      angle: w + f,
    };
  }

  // ---------------------------------------------------------- transits
  /** Fraction of the stellar disk area covered by a planet of radius p at
   *  projected separation z, both in stellar radii (exact circle overlap). */
  function overlapFraction(p, z) {
    if (z >= 1 + p) return 0;
    if (p >= 1 && z <= p - 1) return 1;
    if (z <= 1 - p) return p * p;
    const k0 = Math.acos(Math.min(1, Math.max(-1, (p * p + z * z - 1) / (2 * p * z))));
    const k1 = Math.acos(Math.min(1, Math.max(-1, (1 - p * p + z * z) / (2 * z))));
    const s = Math.max(0, 4 * z * z - Math.pow(1 + z * z - p * p, 2));
    return (p * p * k0 + k1 - 0.5 * Math.sqrt(s)) / Math.PI;
  }

  function intensity(r, u1, u2) {
    const mu = Math.sqrt(Math.max(0, 1 - r * r));
    return 1 - u1 * (1 - mu) - u2 * (1 - mu) * (1 - mu);
  }

  /** Fractional flux drop for planet radius p at separation z with quadratic
   *  limb darkening: small-planet approximation (Mandel & Agol 2002) combined
   *  with the exact overlap area, so ingress and egress have the right shape. */
  function fluxDeficit(p, z, u1, u2) {
    const area = overlapFraction(p, z);
    if (area <= 0) return 0;
    const norm = 1 - u1 / 3 - u2 / 6;
    const rEff = z <= 1 - p ? z : Math.min(1, (1 + Math.max(0, z - p)) / 2);
    return (area * intensity(rEff, u1, u2)) / norm;
  }

  /** Relative stellar flux at times `ts` for a list of planets. */
  function lightCurve(ts, planets, ld) {
    const u1 = ld && ld.enabled ? ld.u1 : 0;
    const u2 = ld && ld.enabled ? ld.u2 : 0;
    const flux = new Float64Array(ts.length);
    for (let i = 0; i < ts.length; i++) {
      let f = 1;
      for (const pl of planets) {
        const pos = position(ts[i], pl);
        if (pos.Z > 0) {                       // only when in front of the star
          const z = Math.hypot(pos.X, pos.Y);
          if (z < 1 + pl.p) f -= fluxDeficit(pl.p, z, u1, u2);
        }
      }
      flux[i] = f;
    }
    return flux;
  }

  /** Impact parameter, durations, depth and transit probability (Winn 2010). */
  function transitParameters(pl, ld) {
    const e = pl.ecc || 0;
    const w = (pl.omega || 0) * deg;
    const inc = pl.inc * deg;
    const eccFactor = (1 - e * e) / (1 + e * Math.sin(w));
    const b = pl.aRs * Math.cos(inc) * eccFactor;
    const transits = b < 1 + pl.p;
    const dur = (k) => {
      if (!(k > b)) return 0;
      const arg = Math.sqrt(k * k - b * b) / (pl.aRs * Math.sin(inc));
      return (pl.period / Math.PI) * Math.asin(Math.min(1, arg)) * (Math.sqrt(1 - e * e) / (1 + e * Math.sin(w)));
    };
    const u1 = ld && ld.enabled ? ld.u1 : 0;
    const u2 = ld && ld.enabled ? ld.u2 : 0;
    return {
      b,
      transits,
      grazing: transits && b > 1 - pl.p,
      t14: transits ? dur(1 + pl.p) : 0,
      t23: transits ? dur(1 - pl.p) : 0,
      depth: transits ? fluxDeficit(pl.p, Math.abs(b), u1, u2) : 0,
      geometricDepth: pl.p * pl.p,
      probability: Math.min(1, (1 + e * Math.sin(w)) / (1 - e * e) / pl.aRs),
    };
  }

  // ---------------------------------------------------- planet properties
  /** Equilibrium temperature [K] with full heat redistribution. Stellar radius in
   *  solar radii and a in AU (converted to the same unit). */
  function equilibriumTemperature(teff, starRadius, aAU, albedo) {
    return teff * Math.sqrt((starRadius * R_SUN_AU) / (2 * aAU)) * Math.pow(1 - albedo, 0.25);
  }

  function planetClass(radiusEarth, teq) {
    if (radiusEarth >= 8) return teq > 1000 ? "Hot Jupiter" : "Gas giant";
    if (radiusEarth >= 4) return "Neptune-like";
    if (radiusEarth >= 2) return "Sub-Neptune";
    if (radiusEarth >= 1.25) return "Super-Earth";
    if (radiusEarth >= 0.8) return "Earth-sized";
    return "Sub-Earth";
  }

  function hzStatus(aAU, hz) {
    if (aAU >= hz.conservativeInner && aAU <= hz.conservativeOuter) return "conservative";
    if (aAU >= hz.optimisticInner && aAU <= hz.optimisticOuter) return "optimistic";
    return aAU < hz.optimisticInner ? "too hot" : "too cold";
  }

  /** Expected transit signal-to-noise for white noise `sigmaPpm` per cadence. */
  function transitSNR(depth, t14Days, t23Days, periodDays, windowDays, cadenceMin, sigmaPpm) {
    if (!(depth > 0)) return 0;
    if (!(sigmaPpm > 0)) return Infinity;
    const nTransits = Math.floor(windowDays / periodDays);
    const nIn = (nTransits * ((t14Days + t23Days) / 2) * 1440) / cadenceMin;
    return (depth / (sigmaPpm * 1e-6)) * Math.sqrt(Math.max(0, nIn));
  }

  // ------------------------------------------------------------ noise
  function rng(seed) {
    let a = seed >>> 0;
    return function () {
      a = (a + 0x6d2b79f5) >>> 0;
      let t = a;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  /** Reproducible Gaussian noise (Box–Muller on a seeded generator). */
  function gaussianNoise(n, sigma, seed) {
    const rand = rng(seed || 12345);
    const out = new Float64Array(n);
    for (let i = 0; i < n; i += 2) {
      const u = Math.max(rand(), 1e-12), v = rand();
      const m = Math.sqrt(-2 * Math.log(u));
      out[i] = sigma * m * Math.cos(2 * Math.PI * v);
      if (i + 1 < n) out[i + 1] = sigma * m * Math.sin(2 * Math.PI * v);
    }
    return out;
  }

  return {
    R_SUN_AU, R_EARTH_RSUN, R_JUP_REARTH, M_JUP_MEARTH, T_SUN,
    luminosity, spectralType, habitableZone, semiMajorAxis, solveKepler, position,
    overlapFraction, intensity, fluxDeficit, lightCurve, transitParameters,
    equilibriumTemperature, planetClass, hzStatus, transitSNR, gaussianNoise,
  };
});
