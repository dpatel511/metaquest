"""Per-stage evidence pages, including explicit skipped/not-run states."""

import html
import json
from pathlib import Path

import pandas as pd

from metaquest.reporting.pages import figure_html, write_page


STEPS = (
    ("validation", "Input validation", ("validation_",)),
    ("preprocessing", "Read preprocessing", ("qc_",)),
    ("taxonomy", "Taxonomic classification", ("taxonomy_",)),
    ("assembly", "Metagenomic assembly", ("assembly_",)),
    ("gene_prediction", "Gene prediction", ("gene_prediction_",)),
    ("annotation", "Functional annotation", ("functional_",)),
    ("abundance", "Gene and functional abundance", ("abundance_",)),
    ("screening", "Research screening", ("screening_",)),
)


def _metric_rows(value, prefix=""):
    if isinstance(value, dict):
        for key, item in value.items():
            yield from _metric_rows(item, f"{prefix} / {key.replace('_', ' ')}".strip(" /"))
    elif isinstance(value, list):
        for number, item in enumerate(value, 1):
            yield from _metric_rows(item, f"{prefix} / {number}")
    else:
        yield {"metric": prefix, "value": "not available" if value is None else value}


def write_step_reports(output: Path, summary: dict) -> list[dict]:
    directory = output / "reports"
    directory.mkdir(exist_ok=True)
    index = []
    for position, (key, label, prefixes) in enumerate(STEPS):
        data = summary.get(key)
        status = "completed" if data else "not_run"
        if data:
            status = data.get("status", status)
            if key == "annotation":
                status = data.get("functional_status", status)
        payload = {"stage": label, "status": status, "metrics": data or {}}
        stem = f"{position:02d}_{key}"
        (directory / f"{stem}.json").write_text(
            json.dumps(payload, indent=2, default=lambda value: value.item()) + "\n", encoding="utf-8")
        frame = pd.DataFrame(_metric_rows(data or {}), columns=["metric", "value"])
        frame.to_csv(directory / f"{stem}.tsv", sep="\t", index=False)
        sections = [("Stage status", f"<p>{html.escape(status)}</p><p><a href='../analysis_metadata.json'>Run provenance and settings</a> · <a href='../report.html'>Consolidated report</a></p>"),
                    ("Recorded evidence", '<div class="scroll">' + frame.head(150).to_html(index=False, escape=True) + '</div>'
                     + '<p>Showing up to 150 metrics; the companion TSV and JSON contain the full record.</p>')]
        if key == "screening" and data:
            for title, path in (("Taxonomic candidates", "taxon_table"), ("Protein evidence", "protein_table")):
                table = pd.read_csv(output / data[path], sep="\t", keep_default_na=False)
                sections.append((title, '<div class="scroll">' + table.head(100).to_html(index=False, escape=True) + '</div>'
                                 + f"<p>Showing up to 100 rows; full table: {html.escape(data[path])}</p>"))
        selected = [output / relative for relative in summary.get("figures", [])
                    if Path(relative).stem.startswith(prefixes) and Path(relative).suffix in (".svg", ".png")]
        selected = [path for path in selected if path.suffix == ".svg" or not path.with_suffix(".svg").is_file()]
        if selected:
            sections.append(("Figures", "".join(figure_html(path, path.stem.replace("_", " ")) for path in selected)))
        write_page(directory / f"{stem}.html", f"MetaQuest · {label}", sections,
                   note="Research evidence. Not run or skipped is not a negative finding. Counts, identities and source provenance should be inspected before biological interpretation.")
        index.append({"stage": label, "status": status, "html": f"reports/{stem}.html",
                      "json": f"reports/{stem}.json", "table": f"reports/{stem}.tsv"})
    return index
