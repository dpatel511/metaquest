import json
from pathlib import Path
from dataclasses import replace
from types import SimpleNamespace

from metaquest.core.preprocessing import _summary_metrics
from metaquest.io.output_formatter import OutputFormatter
from metaquest.pipeline.context import PipelineContext
from metaquest.pipeline.runner import PipelineRunner
from metaquest.reporting.stable_reporter import generate_stable_reports
from metaquest.settings import load_config


def test_fastp_nested_version_and_paired_mean_length():
    result = _summary_metrics(
        {
            "summary": {
                "fastp_version": "0.24.0",
                "before_filtering": {
                    "total_reads": 200,
                    "read1_mean_length": 100,
                    "read2_mean_length": 98,
                },
                "after_filtering": {
                    "total_reads": 150,
                    "read1_mean_length": 95,
                    "read2_mean_length": 93,
                },
            }
        }
    )

    assert result["fastp_version"] == "0.24.0"
    assert result["after"]["mean_length"] == 94


def test_relative_output_and_embedded_html_figures(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "relative-results"
    output.mkdir()
    context = SimpleNamespace(
        output_dir=Path("relative-results"),
        completed_stages=[],
        config=load_config(db_dir=tmp_path / "db"),
        read_mode="single",
        classification=None,
        assembly=None,
        annotation=None,
        preprocessing={
            "before": {"reads": 200, "q20_rate": 0.9, "q30_rate": 0.8},
            "after": {
                "reads": 150,
                "mean_length": 94,
                "q20_rate": 0.98,
                "q30_rate": 0.95,
            },
            "fastp_version": "0.24.0",
            "filtering_result": {"low_quality_reads": 50},
        },
    )

    generate_stable_reports(context)

    report = (output / "report.html").read_text()
    assert "data:image/svg+xml;base64," in report
    assert 'src="plots/' not in report
    assert (output / "analysis_summary.json").is_file()


def test_resume_preserves_log_and_attempt_history(tmp_path):
    log_path = tmp_path / "metaquest.log"
    initial = OutputFormatter(log_file=log_path)
    initial._log("original run")
    initial.log_handle.close()
    resumed = OutputFormatter(log_file=log_path, append_log=True)
    resumed._log("resumed run")
    resumed.log_handle.close()

    log = log_path.read_text()
    assert "original run" in log
    assert "=== MetaQuest resumed run ===" in log
    assert "resumed run" in log

    previous = {
        "status": "completed",
        "started_at": "2026-01-01T00:00:00",
        "timestamp": "2026-01-01T00:10:00",
        "stage_timings_seconds": {"Reporting": 10.0},
        "completed_stages": ["Reporting"],
    }
    (tmp_path / "analysis_metadata.json").write_text(json.dumps(previous))
    config = load_config(db_dir=tmp_path / "databases")
    context = PipelineContext(
        config=config,
        input_files=[],
        output_dir=tmp_path,
        resume=True,
    )
    runner = PipelineRunner(config)
    runner.add_stage("Reporting", lambda ctx: ctx)

    runner.run(context)

    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    assert metadata["started_at"] == "2026-01-01T00:00:00"
    assert metadata["original_stage_timings_seconds"]["Reporting"] == 10.0
    assert metadata["stage_timings_seconds"]["Reporting"] == 10.0
    assert metadata["resume_stage_timings_seconds"]["Reporting"] >= 0
    assert [item["status"] for item in metadata["attempt_history"]] == [
        "completed",
        "completed",
    ]


def test_stage_reports_keep_skipped_annotation_distinct_from_zero_hits(tmp_path):
    config = load_config(db_dir=tmp_path / "db")
    config = replace(config, reporting=replace(config.reporting, plot_formats=("png",)))
    ctx = PipelineContext(config, [], tmp_path, skip_functional=True)
    from metaquest.pipeline.context import AnnotationResult

    prediction = tmp_path / "gene_prediction"
    prediction.mkdir()
    (prediction / "genes.gff3").write_text("c\tfixture\tCDS\t1\t90\t.\t+\t0\tID=g;partial=00\n")
    ctx.annotation = AnnotationResult(prediction, gene_count=1)
    ctx.metadata["validation"] = {"status": "skipped", "bypassed": True}
    generate_stable_reports(ctx)
    summary = json.loads((tmp_path / "analysis_summary.json").read_text())
    status = {item["stage"]: item["status"] for item in summary["stage_reports"]}
    assert status["Input validation"] == "skipped"
    assert status["Gene prediction"] == "completed"
    assert status["Functional annotation"] == "skipped"
    assert status["Research screening"] == "not_run"
    assert summary["annotation"]["annotated_genes"] is None
    assert summary["annotation"]["unannotated_genes"] is None
    assert "data:image/png;base64," in (tmp_path / "reports/04_gene_prediction.html").read_text()
    assert "not assessed" in (tmp_path / "report.html").read_text()


def test_identical_contig_lengths_can_be_plotted(tmp_path):
    from metaquest.pipeline.context import AssemblyResult

    ctx = PipelineContext(load_config(db_dir=tmp_path / "db"), [], tmp_path)
    ctx.assembly = AssemblyResult(tmp_path / "contigs.fa", lengths=[200, 200, 200])
    generate_stable_reports(ctx)
    assert (tmp_path / "plots/assembly_contig_lengths.svg").is_file()


def test_validation_records_full_md5_checksum(tmp_path):
    import hashlib
    from metaquest.io.file_validator import FileValidator

    reads = tmp_path / "reads.fastq"
    reads.write_bytes(b"@one\nACGT\n+\nIIII\n")
    stats = FileValidator()._generate_statistics(str(reads))
    assert stats["md5_checksum"] == hashlib.md5(reads.read_bytes()).hexdigest()
    ctx = PipelineContext(load_config(db_dir=tmp_path / "db"), [], tmp_path)
    ctx.metadata["validation"] = {"status": "passed", "reports": [stats]}
    generate_stable_reports(ctx)
    import pandas as pd

    quality = pd.read_csv(tmp_path / "plots/plot_data/validation_quality.tsv", sep="\t")
    assert quality["percent"].tolist() == [100, 100]
    assert "data:image/svg+xml;base64," in (tmp_path / "reports/00_validation.html").read_text()
