"""Regenerate docs/ screenshots:  python tools/screenshots.py"""
import sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from browser import serve, open_page, ROOT
from playwright.sync_api import sync_playwright

OUT = ROOT / "docs"
serve(port=8767)
with sync_playwright() as p:
    b, page, errors = open_page(p, "http://127.0.0.1:8767/index.html", 1440, 1000)
    page.click("#noiseEnabled"); page.wait_for_timeout(600)
    page.screenshot(path=str(OUT / "screenshot.png"))
    page.screenshot(path=str(OUT / "screenshot_full.png"), full_page=True)
    page.select_option("#preset", "hd209458"); page.wait_for_timeout(600)
    page.evaluate("window.scrollTo(0, 0)"); page.wait_for_timeout(200)
    page.screenshot(path=str(OUT / "screenshot_hd209458.png"))
    b.close()
    b, page, errors = open_page(p, "http://127.0.0.1:8767/index.html", 390, 844)
    page.screenshot(path=str(OUT / "screenshot_mobile.png"))
    b.close()
print("done", errors)
