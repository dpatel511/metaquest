"""Human-readable optional host-filter evidence shared by both report views."""
import html


def host_html(evidence, prefix=""):
    state = evidence["state"]
    fraction = evidence["retained_fraction"]
    rows = {
        "Count unit": evidence["count_unit"], "Input to host filtering": evidence["input"],
        "Retained for downstream analysis": evidence["retained"], "Removed": evidence["removed"],
        "Retained fraction": f"{fraction:.2%}" if fraction is not None else "not available",
        "Bowtie2 version": state["tool_versions"]["bowtie2"],
        "Samtools version": state["tool_versions"]["samtools"],
        "Alignment preset": state["preset"], "Reference index prefix": state["index"]["prefix"],
        "Paired-read policy": "Remove the whole pair if either mate aligns",
    }
    table = "".join(f"<tr><th>{html.escape(key)}</th><td>{html.escape(str(value))}</td></tr>"
                    for key, value in rows.items())
    notes = "".join(f"<li>{html.escape(note)}</li>" for note in evidence["interpretation"])
    return (f"<table>{table}</table><p>Taxonomy, assembly and gene abundance use retained reads. "
            f"<a href='{prefix}preprocessing/host_removal/summary.json'>Host-filter methods, counts and source hashes JSON</a></p>"
            f"<ul>{notes}</ul>")
