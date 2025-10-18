# Exoplanet Transit Simulator (Fixed Edition)

A free, open-source web application that simulates exoplanet transits, models light curves, and analyzes habitable zones. Perfect for education, research, and astronomy enthusiasts!

## Live Demo
After you publish with GitHub Pages, your app will be at:
`https://Jiqbal1921/exoplanet-transit-simulator`

## Features
- Real exoplanet presets (sample of 5 included; add more in `exoplanet_database.json`).
- Multiple-planet systems (up to 4).
- Eccentric orbits via a robust Kepler solver.
- Real-time light-curve modeling.
- Conservative & optimistic habitable-zone estimates (Kopparapu et al. 2013-style).
- Planet classification and equilibrium temperature estimates.

## Quick Start
```bash
git clone https://github.com/Jiqbal1921/exoplanet-transit-simulator.git
cd exoplanet-transit-simulator
open index.html  # or just double-click index.html
```

## Deploy on GitHub Pages
```bash
git init
git add .
git commit -m "Initial commit: Exoplanet Transit Simulator (fixed)"
git branch -M main
git remote add origin https://github.com/Jiqbal1921/exoplanet-transit-simulator.git
git push -u origin main
```
Then enable **Settings → Pages → Build from main**.

## License
MIT — see `LICENSE`.
