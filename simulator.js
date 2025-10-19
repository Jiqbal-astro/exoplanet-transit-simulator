// Global variables
let exoplanetDatabase = [];
let lightCurveChart;

// Initialize on load
document.addEventListener('DOMContentLoaded', function() {
  updatePlanetInputs();
  loadExoplanetDatabase();
  initializeChart();
});

// Chart initialization
function initializeChart() {
  try { if (typeof ChartZoom !== 'undefined') { Chart.register(ChartZoom); } } catch(e) { console.warn('Zoom plugin not registered', e); }
  const ctx = document.getElementById('lightCurveChart').getContext('2d');
  lightCurveChart = new Chart(ctx, {
    type: 'line',
    data: {
      datasets: [{
        label: 'Relative Flux',
        borderColor: '#4CAF50',
        backgroundColor: 'rgba(76, 175, 80, 0.1)',
        borderWidth: 2,
        fill: true,
        tension: 0.4
      }]
    },
    options: {
      responsive: true,
      maintainAspectRatio: false,
      scales: {
        x: {
          title: { display: true, text: 'Time (days)' },
          grid: { color: 'rgba(255, 255, 255, 0.1)' }
        },
        y: {
          title: { display: true, text: 'Relative Flux' },
          suggestedMin: 0.998,
          suggestedMax: 1.0015,
          grid: { color: 'rgba(255, 255, 255, 0.1)' }
        }
      },
      plugins: {
        legend: { labels: { color: 'white' } },
        tooltip: {
          mode: 'index',
          intersect: false,
          callbacks: {
            afterBody: (items) => {
              if (!items || !items.length) return;
              const i = items[0].dataIndex;
              const tp = currentLightCurveData[i]?.transitingPlanets || [];
              return tp.length ? `Transiting: ${tp.join(', ')}` : 'No transit';
            }
          }
        },
        zoom: {
          limits: { y: { min: 0.8, max: 1.05 } },
          pan: { enabled: true, mode: 'xy' },
          zoom: { wheel: { enabled: true }, pinch: { enabled: true }, mode: 'xy' }
        }
      }
    }
  });
}

// Real Exoplanet Database
async function loadExoplanetDatabase() {
  try {
    const response = await fetch('exoplanet_database.json');
    exoplanetDatabase = await response.json();
    populateExoplanetSelector();
    console.log(`Loaded ${exoplanetDatabase.length} exoplanets`);
  } catch (error) {
    console.log('Using default exoplanet data');
    exoplanetDatabase = getDefaultExoplanetData();
    populateExoplanetSelector();
  }
}

function getDefaultExoplanetData() {
  return [
    {
      name: "TRAPPIST-1e",
      discovery_year: 2017,
      star_mass: 0.089,
      star_radius: 0.119,
      star_temp: 2550,
      planet_radius: 0.92,
      orbital_period: 6.10,
      semi_major_axis: 0.029,
      eccentricity: 0.005,
      inclination: 89.86,
      mass: 0.692,
      description: "Potentially habitable Earth-sized exoplanet in TRAPPIST-1 system"
    },
    {
      name: "Kepler-186f",
      discovery_year: 2014,
      star_mass: 0.54,
      star_radius: 0.52,
      star_temp: 3755,
      planet_radius: 1.11,
      orbital_period: 129.9,
      semi_major_axis: 0.432,
      eccentricity: 0.04,
      inclination: 89.9,
      mass: null,
      description: "First Earth-sized planet discovered in habitable zone"
    },
    {
      name: "Proxima Centauri b",
      discovery_year: 2016,
      star_mass: 0.12,
      star_radius: 0.14,
      star_temp: 3042,
      planet_radius: 1.07,
      orbital_period: 11.186,
      semi_major_axis: 0.0485,
      eccentricity: 0.109,
      inclination: 89.5,
      mass: 1.07,
      description: "Closest known exoplanet to Earth"
    },
    {
      name: "HD 209458 b",
      discovery_year: 1999,
      star_mass: 1.148,
      star_radius: 1.203,
      star_temp: 6065,
      planet_radius: 1.38,
      orbital_period: 3.5247,
      semi_major_axis: 0.04707,
      eccentricity: 0.007,
      inclination: 86.71,
      mass: 0.73,
      description: "First exoplanet discovered by transit method"
    },
    {
      name: "Kepler-22b",
      discovery_year: 2011,
      star_mass: 0.97,
      star_radius: 0.98,
      star_temp: 5518,
      planet_radius: 2.38,
      orbital_period: 289.9,
      semi_major_axis: 0.849,
      eccentricity: 0.0,
      inclination: 89.7,
      mass: null,
      description: "First Kepler planet found in the habitable zone"
    }
  ];
}

function populateExoplanetSelector() {
  const selector = document.getElementById('exoplanetSelector');
  // Clear existing options except first two
  while (selector.children.length > 2) {
    selector.removeChild(selector.lastChild);
  }
  
  exoplanetDatabase.forEach(planet => {
    const option = document.createElement('option');
    option.value = planet.name;
    option.textContent = `${planet.name} (${planet.discovery_year})`;
    selector.appendChild(option);
  });
}

function loadExoplanetData() {
  const selector = document.getElementById('exoplanetSelector');
  const selectedPlanetName = selector.value;
  
  if (selectedPlanetName === 'custom' || !selectedPlanetName) {
    return;
  }
  
  const planet = exoplanetDatabase.find(p => p.name === selectedPlanetName);
  if (planet) {
    // Set star properties
    document.getElementById('starMass').value = planet.star_mass;
    document.getElementById('starRadius').value = planet.star_radius;
    document.getElementById('starTemp').value = planet.star_temp;
    
    // Set single planet
    document.getElementById('planetCount').value = '1';
    updatePlanetInputs();
    
    // Set planet properties
    document.getElementById('planetRadius_0').value = planet.planet_radius;
    document.getElementById('period_0').value = planet.orbital_period;
    document.getElementById('semiMajorAxis_0').value = planet.semi_major_axis;
    document.getElementById('inclination_0').value = planet.inclination;
    document.getElementById('eccentricity_0').value = planet.eccentricity;
    document.getElementById('argPeriastron_0').value = 0;
    
    simulateTransit();
  }
}

// Multiple Planet Systems
function updatePlanetInputs() {
  const planetCount = parseInt(document.getElementById('planetCount').value);
  const container = document.getElementById('planetInputsContainer');
  
  container.innerHTML = '<h3> Planet Properties</h3>';
  
  for (let i = 0; i < planetCount; i++) {
    const planetGroup = document.createElement('div');
    planetGroup.className = 'planet-group';
    planetGroup.innerHTML = `
      <h4>Planet ${i + 1}</h4>
      <div class="input-grid">
        <label>Radius (R⊕): <input type="number" id="planetRadius_${i}" value="${1.0 + i * 0.2}" step="0.1" min="0.1" max="20"></label>
        <label>Orbital Period (days): <input type="number" id="period_${i}" value="${365 / (i + 1)}" step="1" min="0.5" max="10000"></label>
        <label>Semi-major Axis (AU): <input type="number" id="semiMajorAxis_${i}" value="${1.0 + i * 0.3}" step="0.1" min="0.01" max="100"></label>
        <label>Inclination (°): <input type="number" id="inclination_${i}" value="${90 - i * 5}" step="1" min="0" max="180"></label>
        <label>Eccentricity: <input type="number" id="eccentricity_${i}" value="0.0" min="0" max="0.99" step="0.01"></label>
        <label>Argument of Periastron (°): <input type="number" id="argPeriastron_${i}" value="0" step="1" min="0" max="360"></label>
        <label>Bond Albedo A: <input type="number" id="albedo_${i}" value="0.30" step="0.01" min="0" max="0.95"></label>
        <label>Greenhouse ε: <input type="number" id="emissivity_${i}" value="1.00" step="0.05" min="0.1" max="1.0"></label>
      </div>
    `;
    container.appendChild(planetGroup);
  }
}


// Limb darkening (quadratic) small-planet approximation
function limbDarkenedDepthSmallPlanet(p, z, u1, u2) {
  // p = Rp/Rs, z = sep/Rs
  if (z > 1 + p) return 0;
  const mu = Math.max(0, Math.sqrt(Math.max(0, 1 - z*z)));
  const I = 1 - u1*(1 - mu) - u2*Math.pow(1 - mu, 2);
  const Iavg = 1 - u1/3 - u2/6;
  return (p*p) * (I / Iavg);
}

// Robust Kepler solver
function solveEccentricAnomaly(M, e, tol = 1e-10, maxIter = 50) {
  M = M % (2*Math.PI);
  if (M < 0) M += 2*Math.PI;
  const s = Math.sin(M), c = Math.cos(M);
  let E = (e < 0.8) ? M + e*s/(1 - e*c) : Math.PI;
  for (let i = 0; i < maxIter; i++) {
    const f = E - e*Math.sin(E) - M;
    const fp = 1 - e*Math.cos(E);
    const fpp = e*Math.sin(E);
    const fppp = e*Math.cos(E);
    const d1 = -f/fp;
    const d2 = -f/(fp + 0.5*d1*fpp);
    const d3 = -f/(fp + 0.5*d2*fpp + (1/6)*d2*d2*fppp);
    E += d3;
    if (Math.abs(d3) < tol) break;
  }
  return E;
}

function calculateTransitLightCurveWithEccentricity(starRadius, planetRadius, semiMajorAxis, inclination, eccentricity, argPeriastron, time, period) {
  const R_star = starRadius * 6.957e8;
  const R_planet = planetRadius * 6.371e6;
  const a = semiMajorAxis * 1.496e11;
  
  const meanAnomaly = ((2 * Math.PI * time) / period) % (2*Math.PI);
  const E = solveEccentricAnomaly(meanAnomaly, eccentricity);
  const trueAnomaly = 2 * Math.atan2(
    Math.sqrt(1 + eccentricity) * Math.sin(E / 2),
    Math.sqrt(1 - eccentricity) * Math.cos(E / 2)
  );
  const distance = a * (1 - eccentricity * eccentricity) / (1 + eccentricity * Math.cos(trueAnomaly));
  
  const argPeriRad = argPeriastron * Math.PI / 180;
  const inclinationRad = inclination * Math.PI / 180;
  
  const x = distance * (Math.cos(argPeriRad) * Math.cos(trueAnomaly) - Math.sin(argPeriRad) * Math.sin(trueAnomaly));
  const y = distance * (Math.sin(argPeriRad) * Math.cos(trueAnomaly) + Math.cos(argPeriRad) * Math.sin(trueAnomaly)) * Math.cos(inclinationRad);
  
  const impactParameter = Math.sqrt(x * x + y * y) / R_star;
  const isTransiting = impactParameter <= (1 + R_planet / R_star);
  const p = R_planet / R_star;
  const depthUniform = p*p;
  let fluxHere = 1.0;
  if (isTransiting) {
    const ldEnabled = document.getElementById('ldToggle')?.checked;
    const u1 = parseFloat(document.getElementById('ldU1')?.value || '0');
    const u2 = parseFloat(document.getElementById('ldU2')?.value || '0');
    if (ldEnabled) {
      const d = limbDarkenedDepthSmallPlanet(p, impactParameter, u1, u2);
      fluxHere = 1 - d;
    } else {
      fluxHere = 1 - depthUniform;
    }
  }
  const depth = depthUniform;
  
  return {
    flux: fluxHere,
    impactParameter: impactParameter,
    distance: distance,
    isTransiting: isTransiting,
    depth: depth
  };
}

// Main Simulation Function
let currentLightCurveData = [];
function simulateTransit() {
  const starMass = parseFloat(document.getElementById('starMass').value);
  const starRadius = parseFloat(document.getElementById('starRadius').value);
  const starTemp = parseFloat(document.getElementById('starTemp').value);
  const planetCount = parseInt(document.getElementById('planetCount').value);
  
  if (!validateInputs(starMass, starRadius, starTemp)) return;
  
  const planets = [];
  for (let i = 0; i < planetCount; i++) {
    const planet = {
      radius: parseFloat(document.getElementById(`planetRadius_${i}`).value),
      period: parseFloat(document.getElementById(`period_${i}`).value),
      semiMajorAxis: parseFloat(document.getElementById(`semiMajorAxis_${i}`).value),
      inclination: parseFloat(document.getElementById(`inclination_${i}`).value),
      eccentricity: parseFloat(document.getElementById(`eccentricity_${i}`).value),
      argPeriastron: parseFloat(document.getElementById(`argPeriastron_${i}`).value),
      albedo: parseFloat(document.getElementById(`albedo_${i}`).value),
      emissivity: parseFloat(document.getElementById(`emissivity_${i}`).value),
      name: `Planet ${i + 1}`
    };
    if (!validatePlanetParameters(planet)) {
      showError(`Invalid parameters for Planet ${i + 1}`);
      return;
    }
    planets.push(planet);
  }
  
  const starLuminosity = Math.pow(starRadius, 2) * Math.pow(starTemp / 5778, 4);
  const hz = calculateHabitableZone(starTemp, starLuminosity);
  const lightCurveData = generateMultiPlanetLightCurve(starRadius, planets);
  currentLightCurveData = lightCurveData;
  
  updateLightCurveChart(lightCurveData);
  updateSystemVisualization(starRadius, planets);
  displayResults(hz, planets, starLuminosity, getSelectedExoplanet());
}

function generateMultiPlanetLightCurve(starRadius, planets) {
  const data = [];
  const points = 200;
  const maxPeriod = Math.max(...planets.map(p => p.period));
  const observationTime = maxPeriod * 1.2;
  
  for (let i = 0; i <= points; i++) {
    const time = (i / points) * observationTime;
    let totalFlux = 1.0;
    const transitingPlanets = [];
    
    planets.forEach((planet, index) => {
      const transit = calculateTransitLightCurveWithEccentricity(
        starRadius, planet.radius, planet.semiMajorAxis, 
        planet.inclination, planet.eccentricity, planet.argPeriastron, 
        time, planet.period
      );
      totalFlux *= transit.flux;
      if (transit.isTransiting) transitingPlanets.push(index + 1);
    });
    
    data.push({ time: time, flux: totalFlux, transitingPlanets: [...transitingPlanets] });
  }
  
  return data;
}

// Validation Functions
function validateInputs(starMass, starRadius, starTemp) {
  if (starMass <= 0 || starRadius <= 0 || starTemp <= 0) {
    showError("Star parameters must be positive values");
    return false;
  }
  return true;
}

function validatePlanetParameters(planet) {
  if (planet.radius <= 0 || planet.period <= 0 || planet.semiMajorAxis <= 0) return false;
  if (planet.eccentricity < 0 || planet.eccentricity >= 1) return false;
  if (planet.inclination < 0 || planet.inclination > 180) return false;
  return true;
}

function showError(message) { alert("Error: " + message); }

function getSelectedExoplanet() {
  const selector = document.getElementById('exoplanetSelector');
  const selectedName = selector.value;
  if (selectedName && selectedName !== 'custom') {
    return exoplanetDatabase.find(p => p.name === selectedName);
  }
  return null;
}

// Habitable Zone Calculations (Kopparapu et al. 2013-inspired coefficients)
function calculateHabitableZone(starTemp, starLuminosity) {
  const T_star = starTemp - 5780;
  const seff_inner = 1.776 + (1.037e-4 * T_star) + (1.235e-8 * T_star * T_star) - (3.764e-12 * T_star * T_star * T_star);
  const innerHZ = Math.sqrt(starLuminosity / seff_inner);
  const seff_outer = 0.3207 + (5.5471e-5 * T_star) + (1.5265e-9 * T_star * T_star) - (2.874e-12 * T_star * T_star * T_star);
  const outerHZ = Math.sqrt(starLuminosity / seff_outer);
  return { inner: innerHZ, outer: outerHZ };
}

function calculateExtendedHabitableZone(starTemp, starLuminosity) {
  const T_star = starTemp - 5780;
  const seff_inner_optimistic = 2.136 + (1.332e-4 * T_star) + (2.343e-8 * T_star * T_star) - (8.844e-12 * T_star * T_star * T_star);
  const innerHZ_optimistic = Math.sqrt(starLuminosity / seff_inner_optimistic);
  const seff_outer_optimistic = 0.2486 + (1.113e-4 * T_star) + (1.399e-8 * T_star * T_star) - (4.895e-12 * T_star * T_star * T_star);
  const outerHZ_optimistic = Math.sqrt(starLuminosity / seff_outer_optimistic);
  return { inner: innerHZ_optimistic, outer: outerHZ_optimistic };
}

// Results Display
function displayResults(hz, planets, starLuminosity, selectedExoplanet = null) {
  const hzResult = document.getElementById('habitableZoneResult');
  
  let hzHTML = `
    <p><strong>Star Type:</strong> ${classifyStar(document.getElementById('starTemp').value, starLuminosity)}</p>
    <p><strong>Star Luminosity:</strong> ${starLuminosity.toFixed(3)} L</p>
    <p><strong>Conservative HZ:</strong> ${hz.inner.toFixed(3)} - ${hz.outer.toFixed(3)} AU</p>
  `;
  
  const extendedHZ = calculateExtendedHabitableZone(document.getElementById('starTemp').value, starLuminosity);
  hzHTML += `<p><strong>Optimistic HZ:</strong> ${extendedHZ.inner.toFixed(3)} - ${extendedHZ.outer.toFixed(3)} AU</p>`;
  
  planets.forEach((planet, index) => {
    const inConservativeHZ = planet.semiMajorAxis >= hz.inner && planet.semiMajorAxis <= hz.outer;
    const inOptimisticHZ = planet.semiMajorAxis >= extendedHZ.inner && planet.semiMajorAxis <= extendedHZ.outer;
    
    let hzStatus = ' Not in habitable zone';
    let hzClass = 'not-habitable';
    
    if (inConservativeHZ) { hzStatus = ' <strong>IN CONSERVATIVE HABITABLE ZONE</strong>'; hzClass = 'habitable'; }
    else if (inOptimisticHZ) { hzStatus = ' In optimistic habitable zone'; hzClass = 'marginally-habitable'; }
    
    const Aplanet = isFinite(planet.albedo) ? planet.albedo : (document.getElementById('albedo') ? parseFloat(document.getElementById('albedo').value) : 0.3);
    const epsPlanet = isFinite(planet.emissivity) ? planet.emissivity : 1.0;
    const planetType = classifyPlanet(planet.radius, planet.semiMajorAxis, document.getElementById('starTemp').value, Aplanet, epsPlanet);
    const eqTemp = calculateEquilibriumTemp(document.getElementById('starTemp').value, document.getElementById('starRadius').value, planet.semiMajorAxis, Aplanet, epsPlanet);
    
    hzHTML += `
      <div class="planet-hz-analysis">
        <h4>${selectedExoplanet ? selectedExoplanet.name : `Planet ${index + 1}`} - ${planetType}</h4>
        <p><strong>Orbital Distance:</strong> ${planet.semiMajorAxis.toFixed(3)} AU</p>
        <p><strong>Equilibrium Temperature:</strong> ${eqTemp.toFixed(0)} K</p>
        <p><strong>A (Bond Albedo):</strong> ${Aplanet.toFixed(2)} | <strong>ε:</strong> ${epsPlanet.toFixed(2)}</p>
        <p><strong>HZ Status:</strong> <span class="${hzClass}">${hzStatus}</span></p>
        ${selectedExoplanet && selectedExoplanet.description ? 
         `<p><em>${selectedExoplanet.description}</em></p>` : 
         `<p><em>Custom planetary simulation</em></p>`}
      </div>
    `;
  });
  
  hzResult.innerHTML = hzHTML;
  displayTransitDetails(planets, selectedExoplanet);
}

function displayTransitDetails(planets, selectedExoplanet) {
  const transitDiv = document.getElementById('transitDetails');
  let transitHTML = '';
  
  if (selectedExoplanet) {
    transitHTML += `
      <div class="real-exoplanet-info">
        <h4>${selectedExoplanet.name} Real Data</h4>
        <p><strong>Discovery Year:</strong> ${selectedExoplanet.discovery_year}</p>
        ${selectedExoplanet.mass ? `<p><strong>Planet Mass:</strong> ${selectedExoplanet.mass} M⊕</p>` : ''}
        <p><strong>Orbital Eccentricity:</strong> ${selectedExoplanet.eccentricity}</p>
      </div>
    `;
  }
  
  planets.forEach((planet, index) => {
    const starRadius = parseFloat(document.getElementById('starRadius').value);
    const RE_OVER_RS = 0.009157; // Earth radii to solar radii
    const depth = Math.pow((planet.radius * RE_OVER_RS) / starRadius, 2);
    
    transitHTML += `
      <div class="planet-transit-info">
        <h4>${selectedExoplanet ? selectedExoplanet.name : `Planet ${index + 1}`} Transit</h4>
        <p><strong>Transit Depth:</strong> ${(depth * 1e6).toFixed(1)} ppm</p>
        <p><strong>Orbital Period:</strong> ${planet.period.toFixed(2)} days</p>
        <p><strong>Semi-major Axis:</strong> ${planet.semiMajorAxis.toFixed(4)} AU</p>
        <p><strong>Planet Radius:</strong> ${planet.radius.toFixed(2)} R⊕</p>
        <p><strong>Eccentricity:</strong> ${planet.eccentricity.toFixed(3)}</p>
      </div>
    `;
  });
  
  transitDiv.innerHTML = transitHTML;
}

// Helper Functions
function classifyStar(temperature, luminosity) {
  if (temperature >= 30000) return 'O-type';
  if (temperature >= 10000) return 'B-type';
  if (temperature >= 7500) return 'A-type';
  if (temperature >= 6000) return 'F-type';
  if (temperature >= 5200) return 'G-type';
  if (temperature >= 3700) return 'K-type';
  return 'M-type';
}

function classifyPlanet(radius, distance, starTemp, albedo, emissivity) {
  const A = (typeof albedo !== 'undefined') ? parseFloat(albedo) : (document.getElementById('albedo') ? parseFloat(document.getElementById('albedo').value) : 0.3);
  const eps = (document.getElementById('emissivity_global') ? parseFloat(document.getElementById('emissivity_global').value) : 1.0);
  const eqTemp = calculateEquilibriumTemp(starTemp, 1, distance, A, eps);
  if (radius < 1) return 'Sub-Earth';
  if (radius < 1.5) return eqTemp > 400 ? 'Hot Earth' : 'Earth-like';
  if (radius < 3) return 'Super-Earth';
  if (radius < 6) return 'Mini-Neptune';
  if (radius < 12) return eqTemp > 1000 ? 'Hot Jupiter' : 'Gas Giant';
  return 'Brown Dwarf';
}

function calculateEquilibriumTemp(starTemp, starRadius, distance, albedo, emissivity) {
  const A = Math.max(0, Math.min(0.99, parseFloat(albedo)));
  const eps = Math.max(0.1, Math.min(1.0, parseFloat(emissivity)) || 1.0);
  // Simple greenhouse via effective emissivity: lower eps -> higher Teq
  return starTemp * Math.sqrt(starRadius / (2 * distance)) * Math.pow((1 - A) / eps, 0.25);
}


// p5.js Orbit View
let p5Instance = null;
function initP5SystemView(planets, observationTimeDays) {
  const container = document.getElementById('systemView');
  container.innerHTML = ''; // clear previous
  const AU_TO_PX = 100;
  const starColor = '#ffd27f';
  const orbitColor = 'rgba(255,255,255,0.3)';
  const planetColor = '#7fc8ff';
  const W = container.clientWidth || 600;
  const H = 300;
  let t = 0;
  let last = null;
  const maxA = Math.max(...planets.map(p => p.semiMajorAxis));
  const scale = Math.min(AU_TO_PX, (Math.min(W, H)/2 - 20) / (maxA * (1 + Math.max(...planets.map(p => p.eccentricity || 0)))));

  p5Instance = new p5(p => {
    p.setup = () => {
      const cnv = p.createCanvas(W, H);
      cnv.parent('systemView');
    };
    p.draw = () => {
      const now = performance.now();
      if (last === null) last = now;
      const dt = (now - last)/1000.0; // seconds
      last = now;
      // advance time so that a full observation loops in ~12 seconds
      const speed = observationTimeDays / 12.0; // days per second
      t = (t + dt * speed) % observationTimeDays;

      p.background(0, 0, 0, 200);
      p.push();
      p.translate(p.width/2, p.height/2);

      // star
      p.noStroke();
      p.fill(starColor);
      p.circle(0, 0, 12);

      // orbits + planets
      planets.forEach(pl => {
        const a = pl.semiMajorAxis * scale;
        const e = pl.eccentricity;
        const b = a * Math.sqrt(1 - e*e);
        // orbit ellipse (center offset by ae to the left)
        p.noFill();
        p.stroke(orbitColor);
        p.ellipse(-a*e, 0, 2*a, 2*b);

        // mean anomaly at time t
        const M = ((2*Math.PI * t) / pl.period) % (2*Math.PI);
        const E = solveEccentricAnomaly(M, e);
        const v = 2 * Math.atan2(Math.sqrt(1+e)*Math.sin(E/2), Math.sqrt(1-e)*Math.cos(E/2));
        // distance r = a(1-e^2)/(1+e cos v)
        const r = a * (1 - e*e) / (1 + e*Math.cos(v));

        // position in ellipse coords (periapsis to the right), rotate by argPeriastron
        const w = (pl.argPeriastron || 0) * Math.PI/180;
        const x = r * Math.cos(v + w) - a*e; // center at star
        const y = r * Math.sin(v + w);

        p.noStroke();
        p.fill(planetColor);
        p.circle(x, y, Math.max(4, 6 - 3*e));
      });

      p.pop();
    };
  });
}

// Visualization Updates
function updateLightCurveChart(lightCurveData) {
  const labels = lightCurveData.map(d => d.time.toFixed(1));
  const data = lightCurveData.map(d => d.flux);
  const minFlux = Math.min(...data);
  
  lightCurveChart.data.labels = labels;
  lightCurveChart.data.datasets[0].data = data;
  lightCurveChart.options.scales.y.suggestedMin = Math.min(0.98, minFlux - 0.0005);
  lightCurveChart.options.scales.y.suggestedMax = 1.0015;
  lightCurveChart.update();
}

function updateSystemVisualization(starRadius, planets) {
  // Simple system visualization placeholder
  const points = 200;
  const maxPeriod = Math.max(...planets.map(p => p.period));
  const observationTime = maxPeriod * 1.2;
  initP5SystemView(planets, observationTime);
}

function resetZoom(){ if(lightCurveChart && lightCurveChart.resetZoom){ lightCurveChart.resetZoom(); } }
