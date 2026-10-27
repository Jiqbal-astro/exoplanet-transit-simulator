/*
 * Planet presets. Units: stellar mass and radius in solar units, planet radius in
 * Earth radii, planet mass in Earth masses, period in days, inclination and
 * argument of periastron in degrees. Semi-major axes are derived from the
 * period and stellar mass (Kepler's third law), so they always agree.
 */
window.PRESETS = [
  {
    id: "trappist1",
    name: "TRAPPIST-1 e, f, g",
    note: "Three temperate Earth-sized planets transiting an ultracool dwarf (Agol et al. 2021). Transit times are illustrative.",
    star: { mass: 0.0898, radius: 0.1192, teff: 2566, u1: 0.168, u2: 0.245 },
    planets: [
      { name: "TRAPPIST-1 e", radius: 0.920, period: 6.099615, inc: 89.793, ecc: 0, omega: 90, t0: 1.0, mass: 0.692 },
      { name: "TRAPPIST-1 f", radius: 1.045, period: 9.206690, inc: 89.740, ecc: 0, omega: 90, t0: 3.0, mass: 1.039 },
      { name: "TRAPPIST-1 g", radius: 1.129, period: 12.35294, inc: 89.742, ecc: 0, omega: 90, t0: 5.0, mass: 1.321 },
    ],
  },
  {
    id: "hd209458",
    name: "HD 209458 b",
    note: "The first exoplanet seen to transit (1999): a hot Jupiter, 1.38 Jupiter radii (Torres et al. 2008; Knutson et al. 2007).",
    star: { mass: 1.148, radius: 1.203, teff: 6065, u1: 0.38, u2: 0.27 },
    planets: [
      { name: "HD 209458 b", radius: 1.38 * 11.209, period: 3.52474859, inc: 86.71, ecc: 0, omega: 90, t0: 1.0, mass: 0.73 * 317.83 },
    ],
  },
  {
    id: "kepler186",
    name: "Kepler-186 f",
    note: "The first Earth-sized planet found in a habitable zone (Quintana et al. 2014).",
    star: { mass: 0.54, radius: 0.52, teff: 3755, u1: 0.45, u2: 0.30 },
    planets: [
      { name: "Kepler-186 f", radius: 1.11, period: 129.944, inc: 89.9, ecc: 0.04, omega: 90, t0: 20, mass: null },
    ],
  },
  {
    id: "kepler22",
    name: "Kepler-22 b",
    note: "The first Kepler planet confirmed in a habitable zone (Borucki et al. 2012).",
    star: { mass: 0.97, radius: 0.979, teff: 5518, u1: 0.42, u2: 0.25 },
    planets: [
      { name: "Kepler-22 b", radius: 2.38, period: 289.862, inc: 89.764, ecc: 0, omega: 90, t0: 40, mass: null },
    ],
  },
  {
    id: "proxima",
    name: "Proxima Centauri b (hypothetical transit)",
    note: "Proxima b was found with radial velocities and has never been seen to transit. Its inclination is unknown; here the orbit is set edge-on and the radius assumed Earth-like to show what a transit would look like (Anglada-Escudé et al. 2016).",
    star: { mass: 0.122, radius: 0.154, teff: 3050, u1: 0.6, u2: 0.25 },
    planets: [
      { name: "Proxima b", radius: 1.07, period: 11.186, inc: 90, ecc: 0, omega: 90, t0: 2.0, mass: 1.07 },
    ],
  },
  {
    id: "earth",
    name: "Earth (seen from afar)",
    note: "What Earth’s transit would look like to a distant observer: an 84 ppm, 13-hour dip once a year.",
    star: { mass: 1.0, radius: 1.0, teff: 5772, u1: 0.40, u2: 0.26 },
    planets: [
      { name: "Earth", radius: 1.0, period: 365.256, inc: 90, ecc: 0.0167, omega: 90, t0: 100, mass: 1.0 },
    ],
  },
];
