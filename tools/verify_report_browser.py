#!/usr/bin/env python3
"""Offline report layout probes; requires Playwright and its Chromium/Firefox runtimes.

Usage: python tools/verify_report_browser.py EVIDENCE_DIR RESULTS_DIR [RESULTS_DIR ...]
Screenshots and Chromium PDFs require visual review; this is not a biological audit.
"""

from pathlib import Path
import json
import sys
from urllib.parse import unquote, urlparse
from playwright.sync_api import sync_playwright

if len(sys.argv) < 3:
    raise SystemExit(__doc__)
out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=True)
results = []
with sync_playwright() as p:
    for name in ("chromium", "firefox"):
        browser = getattr(p, name).launch(headless=True)
        for root in map(Path, sys.argv[2:]):
            for path in [root / "report.html", *sorted((root / "reports").glob("*.html"))]:
                for mode, width in [("desktop", 1440), ("narrow", 375), ("print", 794)]:
                    page = browser.new_page(viewport={"width": width, "height": 950}, offline=True)
                    errors = []
                    page.on("pageerror", lambda e: errors.append(str(e)))
                    page.goto(path.resolve().as_uri())
                    page.emulate_media(media="print" if mode == "print" else "screen")
                    page.locator("img").evaluate_all("(imgs)=>Promise.all(imgs.map(i=>i.decode()))")
                    state = page.evaluate(
                        """()=>({width:innerWidth, scroll:document.documentElement.scrollWidth,brokenImages:[...document.images].filter(i=>!i.complete||!i.naturalWidth).length,external:[...document.querySelectorAll('script[src],link[rel="stylesheet"],img')].map(e=>e.src||e.href).filter(s=>!s.startsWith('data:')),text:document.body.innerText.length})"""
                    )
                    broken_links = []
                    for href in page.locator("a[href]").evaluate_all("els=>els.map(e=>e.href)"):
                        target = urlparse(href)
                        if target.scheme == "file" and not Path(unquote(target.path)).is_file():
                            broken_links.append(href)
                    state["broken_links"] = broken_links
                    state.update(
                        browser=name,
                        version=browser.version,
                        root=str(root),
                        page=path.name,
                        mode=mode,
                        errors=errors,
                    )
                    results.append(state)
                    page.screenshot(
                        path=str(out / f"{name}-{root.name}-{path.stem}-{mode}.png"), full_page=True
                    )
                    if name == "chromium" and mode == "print":
                        page.pdf(
                            path=str(out / f"{root.name}-{path.stem}.pdf"),
                            format="A4",
                            print_background=True,
                        )
                    page.close()
        browser.close()
(out / "results.json").write_text(json.dumps(results, indent=2))
print(
    json.dumps(
        {
            "pages": len(results),
            "problems": [
                r
                for r in results
                if r["scroll"] > r["width"]
                or r["brokenImages"]
                or r["external"]
                or r["errors"]
                or r["broken_links"]
            ],
        },
        indent=2,
    )
)

if any(
    r["scroll"] > r["width"]
    or r["brokenImages"]
    or r["external"]
    or r["errors"]
    or r["broken_links"]
    for r in results
):
    raise SystemExit(1)
