"""End-to-end browser test: operates every control of the dashboard in headless
Chromium and checks the results.   Run:  python tests/test_browser.py
(needs `pip install playwright` and a local copy of plotly.min.js, see tools/browser.py)."""
import sys, pathlib, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "tools"))
from browser import serve, open_page
from playwright.sync_api import sync_playwright

checks = []
def check(name, cond, detail=""):
    checks.append((name, bool(cond), detail))
    print(("PASS " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else ""))

def planets(page):
    return page.evaluate("""TransitApp.derived.planets.map(d => ({name: d.name, depth: d.tp.depth, geom: d.tp.geometricDepth,
        t14h: d.tp.t14 * 24, transits: d.tp.transits, teq: d.teq, a: d.a, hz: d.hz}))""")

def set_slider(page, sel, value):
    page.eval_on_selector(sel, "(el, v) => { el.value = v; el.dispatchEvent(new Event('input', {bubbles: true})); }", value)
    page.wait_for_timeout(250)

serve()
with sync_playwright() as p:
    b, page, errors = open_page(p, "http://127.0.0.1:8765/index.html")

    # ---- presets
    # geometric depth (Rp/R*)^2 and HZ status against literature values
    expected = {"hd209458": ("geom", 0.0139, 0.0006), "earth": ("geom", 83.9e-6, 1e-6),
                "trappist1": ("t14h", 0.93, 0.05), "kepler186": ("hz", "conservative", None)}
    for pid in page.evaluate("PRESETS.map(p => p.id)"):
        page.select_option("#preset", pid)
        page.wait_for_timeout(400)
        pl = planets(page)
        n_dips = page.evaluate("(() => { const f = TransitApp.lc.flux; let n = 0, inn = false;"
                               "for (const x of f) { if (x < 1 && !inn) { n++; inn = true; } if (x === 1) inn = false; } return n; })()")
        check(f"preset {pid}: loads, transits visible in light curve", n_dips >= 1 and all(x["transits"] for x in pl), f"{n_dips} dips")
        if pid in expected:
            key, val, tol = expected[pid]
            got = pl[0][key]
            ok = got == val if tol is None else abs(got - val) <= tol
            check(f"preset {pid}: {key} correct", ok, f"got {got}, expected {val}")
        page.wait_for_timeout(100)

    # duration measured from the simulated light curve matches the analytic T14
    page.select_option("#preset", "kepler22"); page.wait_for_timeout(300)
    measured = page.evaluate("""(() => { const {ts, flux, cadence} = TransitApp.lc; const t0 = TransitApp.state.planets[0].t0;
        let a = null, b = null; for (let i = 0; i < ts.length; i++) if (flux[i] < 1 && Math.abs(ts[i] - t0) < 2) { if (a === null) a = ts[i]; b = ts[i]; }
        return (b - a) * 24 + cadence / 60; })()""")
    analytic = planets(page)[0]["t14h"]
    check("Kepler-22 b: light-curve duration matches analytic T14", abs(measured - analytic) < 0.6, f"{measured:.2f} h vs {analytic:.2f} h")

    # one dip per orbit: HD 209458 b over its window
    page.select_option("#preset", "hd209458"); page.wait_for_timeout(300)
    window, period = page.evaluate("[TransitApp.state.obs.window, TransitApp.state.planets[0].period]")
    n = page.evaluate("(() => { const f = TransitApp.lc.flux; let n = 0, inn = false; for (const x of f) { if (x < 1 && !inn) { n++; inn = true; } if (x === 1) inn = false; } return n; })()")
    t0 = page.evaluate("TransitApp.state.planets[0].t0")
    expected_n = int((window - t0) // period) + 1
    check("exactly one dip per orbit (no false occultation dips)", n == expected_n, f"{n} dips, expected {expected_n}")

    # ---- star sliders change derived values
    page.select_option("#preset", "earth"); page.wait_for_timeout(300)
    teq0 = planets(page)[0]["teq"]
    check("Earth T_eq is 255 K", abs(teq0 - 255) < 2, f"{teq0:.1f} K")
    set_slider(page, "#ctl-teff", 6500)
    check("star temperature slider raises T_eq", planets(page)[0]["teq"] > teq0 + 20)
    d0 = planets(page)[0]["depth"]
    set_slider(page, "#ctl-srad", 0)   # log10(1) = 0 -> 1 R_sun; then shrink
    set_slider(page, "#ctl-srad", -0.3)
    check("star radius slider deepens transit", planets(page)[0]["depth"] > d0 * 3)
    set_slider(page, "#ctl-smass", 0.3)
    check("star mass slider changes semi-major axis (Kepler III)", planets(page)[0]["a"] > 1.2)

    # ---- planet sliders
    page.select_option("#preset", "hd209458"); page.wait_for_timeout(300)
    base = planets(page)[0]
    set_slider(page, "#ctl-prad", 0.5)   # 10^0.5 = 3.2 R_earth
    check("planet radius slider changes depth", planets(page)[0]["depth"] < base["depth"] / 10)
    set_slider(page, "#ctl-period", 1.0)  # 10 days
    check("period slider changes a and duration", planets(page)[0]["a"] > base["a"] * 1.5)
    set_slider(page, "#ctl-inc", 84)
    zoom_text = page.inner_text("#zoomGrid")
    check("low inclination: no transit and an explanation", not planets(page)[0]["transits"] and "No transit" in zoom_text)
    set_slider(page, "#ctl-inc", 90)
    set_slider(page, "#ctl-ecc", 0.4)
    set_slider(page, "#ctl-omega", 30)
    check("eccentricity and periastron sliders accepted, still transiting", planets(page)[0]["transits"])
    set_slider(page, "#ctl-t0", 2)
    check("mid-transit slider moves the transit", abs(page.evaluate("TransitApp.state.planets[0].t0") - 2) < 0.05)
    teq_a = planets(page)[0]["teq"]
    set_slider(page, "#ctl-albedo", 0.8)
    check("albedo slider lowers T_eq", planets(page)[0]["teq"] < teq_a)
    page.fill(".name-field input", "My planet"); page.wait_for_timeout(200)
    check("planet name field renames the tab", "My planet" in page.inner_text("#planetTabs"))

    # ---- add / remove planets
    page.click("#addPlanet"); page.wait_for_timeout(300)
    page.click("#addPlanet"); page.wait_for_timeout(300)
    page.click("#addPlanet"); page.wait_for_timeout(300)
    check("Add planet: up to 4 planets", page.evaluate("TransitApp.state.planets.length") == 4)
    check("Add planet disabled at 4", page.is_disabled("#addPlanet"))
    page.click(".planet-tab >> nth=1"); page.wait_for_timeout(200)
    check("planet tabs switch the selected planet", page.evaluate("TransitApp.state.sel") == 1)
    page.click("#removePlanet"); page.wait_for_timeout(300)
    check("Remove deletes a planet", page.evaluate("TransitApp.state.planets.length") == 3)

    # ---- limb darkening
    page.select_option("#preset", "earth"); page.wait_for_timeout(300)
    with_ld = planets(page)[0]["depth"]
    page.click("#ldEnabled"); page.wait_for_timeout(300)
    without_ld = planets(page)[0]["depth"]
    check("limb-darkening toggle changes central depth", with_ld > without_ld * 1.05, f"{with_ld*1e6:.1f} vs {without_ld*1e6:.1f} ppm")
    page.click("#ldEnabled"); page.wait_for_timeout(200)
    set_slider(page, "#ctl-u1", 0.8)
    check("u1 slider changes depth", planets(page)[0]["depth"] != with_ld)

    # ---- observation and noise
    set_slider(page, "#ctl-cadence", 1.0)
    check("cadence slider sets 10-min cadence", abs(page.evaluate("TransitApp.state.obs.cadence") - 10) < 0.05
          and page.evaluate("TransitApp.lc.cadence") >= page.evaluate("TransitApp.state.obs.cadence") - 1e-9)
    set_slider(page, "#ctl-window", 3)   # 1000 days
    check("time-span slider extends the light curve", page.evaluate("TransitApp.lc.ts[TransitApp.lc.ts.length-1]") > 900)
    check("very long spans coarsen the cadence to stay fast", page.evaluate("TransitApp.lc.capped") and "Cadence increased" in page.inner_text("#lcPlot"))
    page.click("#noiseEnabled"); page.wait_for_timeout(400)
    check("noise toggle adds simulated data", page.evaluate("TransitApp.lc.data !== null"))
    snr_cell = page.inner_text("#resultsTable tbody tr td:last-child")
    check("SNR shown when noise is on", any(c.isdigit() for c in snr_cell), snr_cell)
    set_slider(page, "#ctl-sigma", 3)    # 1000 ppm
    check("noise slider sets 1000 ppm", abs(page.evaluate("TransitApp.state.obs.sigma") - 1000) < 1)

    # ---- time controls
    page.select_option("#preset", "trappist1"); page.wait_for_timeout(400)
    t_start = page.evaluate("TransitApp.state.t")
    page.click("#playBtn"); page.wait_for_timeout(1500)
    t_mid = page.evaluate("TransitApp.state.t")
    page.click("#playBtn"); page.wait_for_timeout(100)
    check("Play animates time", t_mid > t_start, f"{t_start:.3f} -> {t_mid:.3f}")
    t_paused = page.evaluate("TransitApp.state.t"); page.wait_for_timeout(500)
    check("Pause stops time", page.evaluate("TransitApp.state.t") == t_paused)
    set_slider(page, "#timeSlider", 12.5)
    check("time slider sets time", abs(page.evaluate("TransitApp.state.t") - 12.5) < 0.02)
    cursor = page.evaluate("document.getElementById('lcPlot').layout.shapes[0].x0")
    check("light-curve cursor follows time", abs(cursor - page.evaluate("TransitApp.state.t")) < 1e-9)
    box = page.locator("#lcPlot .nsewdrag").first.bounding_box()
    page.mouse.click(box["x"] + box["width"] * 0.25, box["y"] + box["height"] * 0.5)
    page.wait_for_timeout(400)
    window = page.evaluate("TransitApp.state.obs.window")
    check("clicking empty plot area jumps to that time", abs(page.evaluate("TransitApp.state.t") - 0.25 * window) < 0.03 * window,
          f"t = {page.evaluate('TransitApp.state.t'):.2f}, expected ~{0.25*window:.2f}")
    page.mouse.move(box["x"] + box["width"] * 0.6, box["y"] + 40); page.mouse.down()
    page.mouse.move(box["x"] + box["width"] * 0.8, box["y"] + 120, steps=5); page.mouse.up(); page.wait_for_timeout(400)
    xr = page.evaluate("document.getElementById('lcPlot')._fullLayout.xaxis.range")
    check("dragging still zooms (and does not jump)", xr[1] - xr[0] < 0.5 * window and abs(page.evaluate("TransitApp.state.t") - 0.25 * window) < 0.03 * window)
    page.focus("body"); page.keyboard.press("Space"); page.wait_for_timeout(300)
    check("space bar toggles play", page.evaluate("TransitApp.state.playing"))
    page.keyboard.press("Space")

    # ---- download and theme
    with page.expect_download() as dl:
        page.click("#downloadCsv")
    path = dl.value.path()
    lines = open(path).read().splitlines()
    check("Download CSV produces the light curve", lines[0].startswith("time_days") and len(lines) > 1000, f"{len(lines)} rows")
    page.click("#themeToggle"); page.wait_for_timeout(300)
    check("theme toggle switches to dark", page.evaluate("document.documentElement.dataset.theme") == "dark")
    page.click("#themeToggle"); page.wait_for_timeout(300)
    check("theme toggle back to white", page.evaluate("document.documentElement.dataset.theme") == "light")


    # ---- typed values in the number boxes
    page.select_option("#preset", "earth"); page.wait_for_timeout(300)
    page.fill("#ctl-period-num", "200"); page.press("#ctl-period-num", "Enter"); page.wait_for_timeout(300)
    check("typing a period sets it exactly", page.evaluate("TransitApp.state.planets[0].period") == 200)
    page.fill("#ctl-period-num", "5000"); page.press("#ctl-period-num", "Enter"); page.wait_for_timeout(300)
    check("typed values beyond the slider range are accepted", page.evaluate("TransitApp.state.planets[0].period") == 5000)
    page.fill("#ctl-ecc-num", "1.5"); page.press("#ctl-ecc-num", "Enter"); page.wait_for_timeout(200)
    check("invalid typed value is rejected and flagged",
          page.evaluate("TransitApp.state.planets[0].ecc") < 1 and "invalid" in page.get_attribute("#ctl-ecc-num", "class"))

    # ---- custom target: radius in Earth radii, inclination
    page.click("#customBtn"); page.wait_for_timeout(200)
    check("Custom target opens a form prefilled from the current system", page.is_visible("#customDialog")
          and page.input_value("[name=teff]") == "5772")
    page.fill("[name=name]", "Test star A")
    page.fill("[name=teff]", "4500"); page.fill("[name=srad]", "0.7"); page.fill("[name=smass]", "0.75")
    row = "#customPlanets fieldset:nth-child(1) "
    page.fill(row + "[name=pname]", "A b"); page.fill(row + "[name=size]", "2.5"); page.select_option(row + "[name=sizeMode]", "re")
    page.fill(row + "[name=period]", "12.3"); page.fill(row + "[name=geom]", "89.5"); page.fill(row + "[name=t0]", "3")
    page.click("#customAddPlanet"); page.wait_for_timeout(100)
    row2 = "#customPlanets fieldset:nth-child(2) "
    page.fill(row2 + "[name=size]", "1.1"); page.select_option(row2 + "[name=sizeMode]", "rj")
    page.fill(row2 + "[name=period]", "30"); page.fill(row2 + "[name=geom]", "0.3"); page.select_option(row2 + "[name=geomMode]", "b")
    page.click("#customDialog button[type=submit]"); page.wait_for_timeout(500)
    st = page.evaluate("JSON.parse(JSON.stringify({s: TransitApp.state.star, p: TransitApp.state.planets}))")
    pl = planets(page)
    check("custom target is modelled with its star", not page.is_visible("#customDialog") and st["s"]["teff"] == 4500 and st["s"]["radius"] == 0.7)
    check("custom planet 1: radius in Earth radii and inclination", abs(st["p"][0]["radius"] - 2.5) < 1e-9 and st["p"][0]["inc"] == 89.5)
    check("custom planet 2: Jupiter radii converted", abs(st["p"][1]["radius"] - 1.1 * 11.209) < 1e-6)
    b2 = page.evaluate("TransitApp.derived.planets[1].tp.b")
    check("custom planet 2: impact parameter converted to inclination", abs(b2 - 0.3) < 1e-6, f"b = {b2:.4f}")
    check("custom target appears under My targets", "Test star A" in page.inner_text("#preset optgroup[label='My targets']"))
    check("custom target transits are drawn", all(x["transits"] for x in pl))

    # ---- custom target: size from transit depth, and validation
    page.click("#customBtn"); page.wait_for_timeout(200)
    page.fill("[name=name]", "Depth target")
    page.select_option(row + "[name=sizeMode]", "depth"); page.fill(row + "[name=size]", "10000")
    page.click("#customDialog button[type=submit]"); page.wait_for_timeout(500)
    check("size given as transit depth reproduces that depth", abs(page.evaluate("TransitApp.derived.planets[0].tp.geometricDepth") - 0.01) < 1e-9)
    page.click("#customBtn"); page.wait_for_timeout(200)
    page.fill(row + "[name=period]", "-4")
    page.click("#customDialog button[type=submit]"); page.wait_for_timeout(200)
    check("invalid custom input shows a clear error and keeps the form open",
          page.is_visible("#customDialog") and "period" in page.inner_text("#customError"), page.inner_text("#customError"))
    page.click("#customCancel"); page.wait_for_timeout(200)
    check("Cancel closes the form", not page.is_visible("#customDialog"))

    # ---- depth clearly visible
    for pid in ("earth", "trappist1", "hd209458", "kepler22"):
        page.select_option("#preset", pid); page.wait_for_timeout(400)
        frac = page.evaluate("""(() => { const el = document.getElementById('lcPlot'); const r = el._fullLayout.yaxis.range;
            const minModel = Math.min(...TransitApp.lc.flux); const u = el._fullLayout.yaxis.title.text;
            const c = u.includes('ppm') ? (f) => (f - 1) * 1e6 : u.includes('%') ? (f) => (f - 1) * 100 : (f) => f;
            return (c(1) - c(minModel)) / (r[1] - r[0]); })()""")
        lines = page.evaluate("document.getElementById('lcPlot').layout.shapes.filter(s => s.line && s.line.dash === 'dash').length")
        check(f"{pid}: transit depth fills most of the light-curve height, with depth line", frac > 0.6 and lines >= 1, f"{frac:.0%} of axis")
    page.check("#noiseEnabled"); page.wait_for_timeout(400)
    names = page.evaluate("document.getElementById('lcPlot').data.map(t => t.name)")
    check("noise on: binned data drawn over the raw points", any(n and n.startswith("Binned") for n in names))
    page.select_option("#unitSelect", "pct"); page.wait_for_timeout(300)
    check("units switch to %", "%" in page.evaluate("document.getElementById('lcPlot')._fullLayout.yaxis.title.text"))
    page.select_option("#unitSelect", "flux"); page.wait_for_timeout(300)
    check("units switch to relative flux", page.evaluate("document.getElementById('lcPlot')._fullLayout.yaxis.title.text") == "Relative flux")
    page.select_option("#scaleSelect", "data"); page.wait_for_timeout(300)
    check("scale switch to fit all data", page.evaluate("document.getElementById('lcPlot').layout.yaxis.autorange") is True)
    page.select_option("#scaleSelect", "transits"); page.select_option("#unitSelect", "auto"); page.uncheck("#noiseEnabled"); page.wait_for_timeout(300)

    check("no JavaScript errors during the whole session", not errors, "; ".join(errors[:3]))
    b.close()

failed = [c for c in checks if not c[1]]
print(f"\n{len(checks) - len(failed)} / {len(checks)} checks passed")
sys.exit(1 if failed else 0)
