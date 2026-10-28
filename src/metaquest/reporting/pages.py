"""Small offline HTML reports with escaped tables and embedded figures."""

import base64
import html
from pathlib import Path


STYLE = """
:root{color-scheme:light;--ink:#193247;--muted:#526779;--line:#dce5ea;--accent:#116b78}
*{box-sizing:border-box}body{margin:0;background:#f3f6f8;color:var(--ink);font:15px/1.6 system-ui,sans-serif}
main{max-width:1220px;margin:auto;padding:36px 28px}header{padding:25px 0;border-bottom:3px solid var(--accent)}
.eyebrow{font-size:12px;letter-spacing:.12em;text-transform:uppercase;color:var(--accent);font-weight:700}
h1{font-size:clamp(26px,4vw,40px);line-height:1.2;margin:12px 0}h2{font-size:22px;margin:0 0 16px}
nav{display:flex;flex-wrap:wrap;gap:8px;margin:20px 0}a{color:#096779}nav a{padding:6px 12px;background:white;border:1px solid var(--line);border-radius:20px;text-decoration:none}
section,figure{background:white;border:1px solid var(--line);border-radius:12px;padding:24px;margin:20px 0}
.note{padding:14px 18px;border-left:4px solid #c79334;background:#fff9e9}.scroll{overflow:auto}
table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:left;border-bottom:1px solid var(--line);padding:9px 12px;vertical-align:top;overflow-wrap:anywhere}
th{background:#eff5f7}img{max-width:100%;height:auto;display:block;margin:auto}figcaption{color:var(--muted);margin-top:12px}
pre{white-space:pre-wrap;overflow-wrap:anywhere}footer{color:var(--muted);padding:20px 0;font-size:12px}
@media print{body{background:white}main{padding:0}section,figure{break-inside:avoid}nav{display:none}}
"""


def figure_html(path: Path, caption: str) -> str:
    mime = "image/svg+xml" if path.suffix == ".svg" else "image/png"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    label = html.escape(caption)
    return f'<figure><img src="data:{mime};base64,{encoded}" alt="{label}"><figcaption>{label}</figcaption></figure>'


def write_page(path: Path, title: str, sections: list[tuple[str, str]], *, note: str = "") -> None:
    navigation = "".join(f'<a href="#section-{i}">{html.escape(name)}</a>' for i, (name, _) in enumerate(sections))
    contents = "".join(f'<section id="section-{i}"><h2>{html.escape(name)}</h2>{body}</section>'
                       for i, (name, body) in enumerate(sections))
    path.write_text(
        '<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{html.escape(title)}</title><style>{STYLE}</style></head><body><main>'
        f'<header><div class="eyebrow">MetaQuest / Research evidence</div><h1>{html.escape(title)}</h1></header>'
        f'<p class="note">{html.escape(note)}</p><nav>{navigation}</nav>{contents}'
        '<footer>Generated from recorded analysis outputs. Inspect source tables and provenance before interpretation.</footer>'
        '</main></body></html>', encoding="utf-8")
