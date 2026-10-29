"""Failure injection uses tiny files and disposable installations only."""
from dataclasses import asdict, replace
import gzip
import hashlib
import io
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import sqlite3
import tarfile
import threading
import urllib.error
from types import SimpleNamespace

import pytest

from metaquest import database_manager as db
from metaquest.core.analysis import _jsonable, _sha256, _validate_resume
from metaquest.exceptions import ConfigError, PipelineStageError
from metaquest.pipeline.context import PipelineContext
from metaquest.pipeline.runner import PipelineRunner
from metaquest.settings import load_config


class Response(io.BytesIO):
    def __init__(self, payload, status, headers):
        super().__init__(payload)
        self.status, self.headers = status, headers


def test_completed_unknown_size_download_handles_http_416(tmp_path, monkeypatch):
    partial = tmp_path / "data.partial"
    partial.write_bytes(b"complete")

    def response(*_, **__):
        raise urllib.error.HTTPError("https://fixture", 416, "range", {"Content-Range": "bytes */8"}, None)

    monkeypatch.setattr(db.urllib.request, "urlopen", response)
    db._download("https://fixture", partial, max_requests=1)
    assert partial.read_bytes() == b"complete"


def test_range_body_cannot_exceed_advertised_range(tmp_path, monkeypatch):
    partial = tmp_path / "data.partial"
    partial.write_bytes(b"abc")
    monkeypatch.setattr(db.urllib.request, "urlopen", lambda *_, **__: Response(
        b"defghij", 206, {"Content-Range": "bytes 3-4/10"}))
    with pytest.raises(db.DatabaseSetupError, match="range"):
        db._download("https://fixture", partial, expected_bytes=10, max_requests=1)
    assert partial.read_bytes() == b"abc"


def test_download_uses_a_network_timeout(tmp_path, monkeypatch):
    def response(*_, **kwargs):
        assert 0 < kwargs.get("timeout", 0) <= 60
        return Response(b"abc", 200, {"Content-Length": "3"})

    monkeypatch.setattr(db.urllib.request, "urlopen", response)
    db._download("https://fixture", tmp_path / "data.partial", expected_bytes=3)


def test_http_disconnect_resumes_against_local_server(tmp_path):
    ranges = []

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            requested = self.headers.get("Range")
            ranges.append(requested)
            self.send_response(206 if requested else 200)
            self.send_header("Content-Length", "3" if requested else "6")
            if requested:
                self.send_header("Content-Range", "bytes 3-5/6")
            self.end_headers()
            # Initial response ends early; the next request supplies the suffix.
            self.wfile.write(b"def" if requested else b"abc")
            self.close_connection = True

        def log_message(self, *_):
            pass

    with ThreadingHTTPServer(("127.0.0.1", 0), Handler) as server:
        thread = threading.Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            path = tmp_path / "data.partial"
            db._download(f"http://127.0.0.1:{server.server_port}/fixture", path, expected_bytes=6)
        finally:
            server.shutdown()
            thread.join(timeout=3)
    assert path.read_bytes() == b"abcdef"
    assert ranges == [None, "bytes=3-5"]


def toy_database(tmp_path, monkeypatch, key):
    target = tmp_path / key
    target.mkdir()
    (target / "preserve.txt").write_text("original installation")
    monkeypatch.setattr(db.shutil, "disk_usage", lambda _: SimpleNamespace(free=100 * 1024**3))
    payloads = {}
    if key == "taxonomy":
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w:gz") as bundle:
            for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
                item = tarfile.TarInfo(name)
                item.size = 5
                bundle.addfile(item, io.BytesIO(b"index"))
        payloads[db.DATABASES[key].url] = buffer.getvalue()
        archive_name = Path(db.DATABASES[key].url).name
        payloads[db.DATABASES[key].checksum_url] = (
            hashlib.md5(buffer.getvalue()).hexdigest() + "  " + archive_name + "\n"
        ).encode()
    else:
        source = tmp_path / "toy.sqlite"
        with sqlite3.connect(source) as connection:
            connection.execute("CREATE TABLE fixture (id INTEGER)")
        artifacts = []
        for name in ("eggnog.db", "eggnog.taxa.db", "eggnog_proteins.dmnd"):
            content = gzip.compress(source.read_bytes() if name.endswith(".db") else b"fixture")
            artifact = {"archive": name + ".gz", "output": name,
                        "format": "gzip", "compressed_bytes": len(content)}
            artifacts.append(artifact)
            payloads[f"{db.DATABASES[key].url}/{artifact['archive']}"] = content
        monkeypatch.setattr(db, "FUNCTIONAL_ARTIFACTS", artifacts)
        monkeypatch.setattr(db, "_validate_functional_tools", lambda *_: None)

    def download(url, path, **_):
        path.write_bytes(payloads[url])

    monkeypatch.setattr(db, "_download", download)
    return target


@pytest.mark.parametrize("key", ["taxonomy", "functional"])
@pytest.mark.parametrize("failure", [OSError, KeyboardInterrupt])
def test_failed_database_promotion_preserves_previous_installation(tmp_path, monkeypatch, key, failure):
    target = toy_database(tmp_path, monkeypatch, key)
    original = Path.replace

    def fail_promotion(path, destination):
        if Path(destination) == target and path.name in ("extracted", "functional"):
            raise failure("injected promotion failure")
        return original(path, destination)

    monkeypatch.setattr(Path, "replace", fail_promotion)
    with pytest.raises(failure):
        db.install_database(key, tmp_path, force=True)
    assert (target / "preserve.txt").read_text() == "original installation"


@pytest.mark.parametrize("key", ["taxonomy", "functional"])
def test_low_space_preserves_database_without_network(tmp_path, monkeypatch, key):
    target = toy_database(tmp_path, monkeypatch, key)
    monkeypatch.setattr(db.shutil, "disk_usage", lambda _: SimpleNamespace(free=0))
    monkeypatch.setattr(db, "_download", lambda *_args, **_kwargs: pytest.fail("network before disk check"))
    with pytest.raises(db.DatabaseSetupError, match="Insufficient"):
        db.install_database(key, tmp_path, force=True)
    assert (target / "preserve.txt").is_file()


def test_recovers_database_backup_left_by_hard_termination(tmp_path, monkeypatch):
    backup = tmp_path / ".taxonomy.previous"
    backup.mkdir()
    (backup / "metaquest-db.json").write_text(json.dumps(asdict(db.DATABASES["taxonomy"])))
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (backup / name).write_text("original")
    monkeypatch.setattr(db, "_download", lambda *_args, **_kwargs: pytest.fail("lost recoverable database"))
    target = db.install_database("taxonomy", tmp_path)
    assert (target / "hash.k2d").read_text() == "original"
    assert not backup.exists()


@pytest.mark.parametrize("manifest", ["invalid JSON", "[]", '{"key":"wrong"}', '{"key":"taxonomy"}'])
def test_invalid_manifest_is_not_an_installed_database(tmp_path, manifest):
    target = tmp_path / "taxonomy"
    target.mkdir()
    (target / "metaquest-db.json").write_text(manifest)
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (target / name).write_text("index")
    assert not db.is_installed("taxonomy", tmp_path)


def test_empty_index_is_not_installed(tmp_path):
    target = tmp_path / "taxonomy"
    target.mkdir()
    (target / "metaquest-db.json").write_text(json.dumps(asdict(db.DATABASES["taxonomy"])))
    for name in ("hash.k2d", "opts.k2d", "taxo.k2d"):
        (target / name).touch()
    assert not db.is_installed("taxonomy", tmp_path)


def test_archive_rejects_named_pipe_before_extraction(tmp_path):
    archive = tmp_path / "pipe.tar.gz"
    with tarfile.open(archive, "w:gz") as bundle:
        member = tarfile.TarInfo("named-pipe")
        member.type = tarfile.FIFOTYPE
        bundle.addfile(member)
    destination = tmp_path / "extracted"
    destination.mkdir()
    with pytest.raises(db.DatabaseSetupError, match="Unsafe"):
        db._safe_extract(archive, destination)


def test_corrupt_gzip_is_discarded_so_retry_can_redownload(tmp_path, monkeypatch):
    target = toy_database(tmp_path, monkeypatch, "functional")
    original = db._download
    bad = db.FUNCTIONAL_ARTIFACTS[0]

    def corrupt_download(url, path, **kwargs):
        original(url, path, **kwargs)
        if path.name == bad["archive"] + ".partial":
            path.write_bytes(b"X" * int(bad["compressed_bytes"]))

    monkeypatch.setattr(db, "_download", corrupt_download)
    with pytest.raises(db.DatabaseSetupError, match="gzip"):
        db.install_database("functional", tmp_path, force=True)
    assert (target / "preserve.txt").is_file()
    assert not (tmp_path / ".downloads/functional" / (bad["archive"] + ".partial")).exists()
    monkeypatch.setattr(db, "_download", original)
    db.install_database("functional", tmp_path, force=True)
    assert db.is_installed("functional", tmp_path)


def test_disk_full_during_extraction_preserves_old_database_and_archives(tmp_path, monkeypatch):
    target = toy_database(tmp_path, monkeypatch, "functional")
    original = Path.open

    def no_space(path, *args, **kwargs):
        if path.name == "eggnog.db" and args and args[0] == "wb":
            raise OSError(28, "No space left on device")
        return original(path, *args, **kwargs)

    monkeypatch.setattr(Path, "open", no_space)
    with pytest.raises(OSError, match="No space"):
        db.install_database("functional", tmp_path, force=True)
    assert (target / "preserve.txt").is_file()
    assert len(list((tmp_path / ".downloads/functional").glob("*.partial"))) == 3


def context(tmp_path, monkeypatch, *, resume=False):
    from metaquest.pipeline import runner
    monkeypatch.setattr(runner, "_tool_version", lambda _: "fixture-v1")
    return PipelineContext(load_config(db_dir=tmp_path / "db"), [], tmp_path,
                           read_mode="single", resume=resume)


def test_interrupt_is_recorded_and_not_completed(tmp_path, monkeypatch):
    ctx = context(tmp_path, monkeypatch)
    pipeline = PipelineRunner(ctx.config)

    def interrupt(_):
        raise KeyboardInterrupt

    pipeline.add_stage("fixture", interrupt)
    with pytest.raises(KeyboardInterrupt):
        pipeline.run(ctx)
    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    assert metadata["status"] == metadata["attempt_history"][-1]["status"] == "interrupted"
    assert metadata["failed_stage"] == "fixture"
    assert "fixture" not in metadata["completed_stages"]


def test_failed_then_resumed_stage_clears_stale_failure_and_preserves_history(tmp_path, monkeypatch):
    ctx = context(tmp_path, monkeypatch)
    pipeline = PipelineRunner(ctx.config)

    def fail(_):
        raise RuntimeError("injected failure")

    pipeline.add_stage("fixture", fail)
    with pytest.raises(PipelineStageError):
        pipeline.run(ctx)
    resumed = context(tmp_path, monkeypatch, resume=True)
    pipeline = PipelineRunner(ctx.config).add_stage("fixture", lambda c: c)
    pipeline.run(resumed)
    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    assert metadata["status"] == "completed"
    assert "failed_stage" not in metadata and "failure_reason" not in metadata
    assert [attempt["status"] for attempt in metadata["attempt_history"]] == ["failed", "completed"]


@pytest.mark.parametrize("damage", [None, "modify", "remove"])
def test_runner_reuses_only_intact_completed_stages(tmp_path, monkeypatch, damage):
    seen = []

    def stage(name, relative):
        def run(ctx):
            seen.append((name, ctx.resume))
            path = ctx.output_dir / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
            return ctx
        return run

    # Real stage names/required artifacts, but no native tools or biological data.
    files = ["megahit_assembly/final.contigs.fa"]
    pipeline = PipelineRunner(load_config(db_dir=tmp_path / "db"))
    pipeline.add_stage("Metagenomic Assembly", stage("assembly", files[0]))
    pipeline.add_stage("downstream", stage("downstream", "downstream.txt"))
    pipeline.run(context(tmp_path, monkeypatch))
    if damage == "modify":
        (tmp_path / files[0]).write_text("damaged")
    elif damage == "remove":
        (tmp_path / files[0]).unlink()
    seen.clear()
    pipeline.run(context(tmp_path, monkeypatch, resume=True))
    assert seen == [("assembly", damage is None), ("downstream", damage is None)]
    assert (tmp_path / files[0]).read_text() == "assembly"


def test_changed_interleaved_split_invalidates_classification_cache(tmp_path, monkeypatch):
    seen = []

    def run(ctx):
        seen.append(ctx.resume)
        for name in ("kraken_report.txt", "bracken_report.tsv", "split_R1.fastq.gz", "split_R2.fastq.gz"):
            (tmp_path / name).write_text("fixture")
        return ctx

    def configured(resume):
        ctx = context(tmp_path, monkeypatch, resume=resume)
        ctx.read_mode = "interleaved"
        ctx.config = replace(ctx.config, preprocessing=replace(ctx.config.preprocessing, enabled=False))
        return ctx

    pipeline = PipelineRunner(load_config(db_dir=tmp_path / "db"))
    pipeline.add_stage("Taxonomic Classification", run)
    pipeline.run(configured(False))
    (tmp_path / "split_R2.fastq.gz").write_text("damaged")
    pipeline.run(configured(True))
    assert seen == [False, False]


def test_resume_rejects_changed_database_manifest(tmp_path, monkeypatch):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@fixture\nACGT\n+\nIIII\n")
    ctx = context(tmp_path, monkeypatch)
    ctx.input_files = [reads]
    database = ctx.config.databases.kraken_db
    database.mkdir(parents=True)
    manifest = database / "metaquest-db.json"
    manifest.write_text('{"release": "original"}')
    PipelineRunner(ctx.config).run(ctx)
    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    manifest.write_text('{"release": "replaced"}')
    with pytest.raises(ConfigError, match="database"):
        _validate_resume(tmp_path, [reads], ctx.config, {
            key: value for key, value in metadata["workflow"].items() if key != "resume"
        })


def test_resume_rejects_changed_inputs_and_config(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("original")
    config = load_config(db_dir=tmp_path / "db")
    metadata = {"input_files": [{"path": str(reads.resolve()), "size_bytes": reads.stat().st_size,
                                  "sha256": _sha256(reads)}],
                "effective_config": _jsonable(asdict(config)), "workflow": {}}
    (tmp_path / "analysis_metadata.json").write_text(json.dumps(metadata))
    reads.write_text("modified")
    with pytest.raises(ConfigError, match="input"):
        _validate_resume(tmp_path, [reads], config, {})
    reads.write_text("original")
    changed = replace(config, assembly=replace(config.assembly, threads=1))
    with pytest.raises(ConfigError, match="configuration"):
        _validate_resume(tmp_path, [reads], changed, {})


def test_resume_rejects_changed_tool_version(tmp_path, monkeypatch):
    from metaquest.pipeline import runner
    ctx = context(tmp_path, monkeypatch)
    PipelineRunner(ctx.config).run(ctx)
    metadata = json.loads((tmp_path / "analysis_metadata.json").read_text())
    monkeypatch.setattr(runner, "_tool_version", lambda _: "fixture-v2")
    with pytest.raises(ConfigError, match="tool version"):
        _validate_resume(tmp_path, [], ctx.config, {
            key: value for key, value in metadata["workflow"].items() if key != "resume"
        })


def test_previous_running_attempt_is_preserved_as_interrupted(tmp_path, monkeypatch):
    metadata = {"status": "running", "attempt_history": [{"status": "running"}]}
    (tmp_path / "analysis_metadata.json").write_text(json.dumps(metadata))
    ctx = context(tmp_path, monkeypatch, resume=True)
    PipelineRunner(ctx.config).run(ctx)
    saved = json.loads((tmp_path / "analysis_metadata.json").read_text())
    assert [attempt["status"] for attempt in saved["attempt_history"]] == ["interrupted", "completed"]


def test_metadata_promotion_failure_leaves_last_complete_json(tmp_path, monkeypatch):
    ctx = context(tmp_path, monkeypatch)
    pipeline = PipelineRunner(ctx.config)
    pipeline._save_metadata(ctx)
    metadata = tmp_path / "analysis_metadata.json"
    before = metadata.read_bytes()
    original = Path.replace

    def fail_metadata(path, destination):
        if Path(destination) == metadata:
            raise OSError(28, "injected full filesystem")
        return original(path, destination)

    monkeypatch.setattr(Path, "replace", fail_metadata)
    ctx.metadata["status"] = "changed"
    with pytest.raises(OSError, match="filesystem"):
        pipeline._save_metadata(ctx)
    assert metadata.read_bytes() == before
    assert isinstance(json.loads(metadata.read_text()), dict)
