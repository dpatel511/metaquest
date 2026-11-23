"""Read preprocessing pipeline stage."""

from __future__ import annotations

import json
from pathlib import Path

from metaquest.pipeline.context import PipelineContext


def run_preprocessing_stage(ctx: PipelineContext) -> PipelineContext:
    from metaquest.core.preprocessing import _summary_metrics, run_fastp
    from metaquest.io.output_formatter import get_formatter

    config = ctx.config.preprocessing
    if not config.enabled:
        ctx.analysis_input_files = list(ctx.input_files)
        ctx.preprocessing = {"status": "disabled"}
        return _remove_host(ctx)

    json_path = ctx.output_dir / "preprocessing" / "fastp.json"
    cleaned = [ctx.output_dir / "preprocessing" / "cleaned_R1.fastq.gz"]
    if ctx.read_mode != "single":
        cleaned.append(ctx.output_dir / "preprocessing" / "cleaned_R2.fastq.gz")
    if ctx.resume and json_path.is_file() and all(path.is_file() for path in cleaned):
        payload = json.loads(json_path.read_text(encoding="utf-8"))
        ctx.preprocessing = _summary_metrics(payload)
        ctx.preprocessing.update({"json_report": str(json_path), "html_report": str(json_path.with_suffix('.html'))})
        get_formatter().info("Reusing completed fastp preprocessing")
    else:
        if ctx.read_mode == "interleaved":
            from metaquest.io.utils import split_interleaved
            source = [Path(p) for p in split_interleaved(str(ctx.input_files[0]), ctx.output_dir, get_formatter())]
        else:
            source = ctx.input_files
        cleaned, ctx.preprocessing = run_fastp(
            source, ctx.output_dir, paired=ctx.read_mode != "single",
            quality=config.qualified_quality_phred,
            min_length=config.length_required,
            threads=config.threads,
        )
    ctx.analysis_input_files = cleaned
    return _remove_host(ctx)


def _remove_host(ctx):
    config = ctx.config.preprocessing
    if config.host_index is None:
        return ctx
    from metaquest.core.host_removal import run_host_removal
    from metaquest.exceptions import MetaQuestError
    from metaquest.io.output_formatter import get_formatter

    reads = ctx.analysis_input_files
    if ctx.read_mode == "interleaved" and len(reads) == 1:
        from metaquest.io.utils import split_interleaved
        split = [ctx.output_dir / f"split_R{mate}.fastq.gz" for mate in (1, 2)]
        reads = split if ctx.resume and all(path.is_file() for path in split) else [
            Path(p) for p in split_interleaved(str(reads[0]), ctx.output_dir, get_formatter())]
    get_formatter().info("Filtering reads against the supplied host index")
    retained, evidence = run_host_removal(reads, ctx.output_dir, config.host_index,
                                        paired=ctx.read_mode != "single", threads=config.threads,
                                        reuse=ctx.resume)
    ctx.preprocessing.update(status="completed", fastp_status="completed" if config.enabled else "disabled",
                             host_removal=evidence)
    ctx.analysis_input_files = retained
    if evidence["retained"] == 0:
        raise MetaQuestError("No reads retained after host removal; inspect preprocessing/host_removal/summary.json")
    return ctx
