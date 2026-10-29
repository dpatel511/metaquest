import gzip
import http.client
import io
import sqlite3
import sys
import tarfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from metaquest import database_manager
from metaquest.core import functional_analysis
from metaquest.core.functional_analysis import run_functional_annotation, run_gene_prediction
from metaquest.database_manager import (
    DATABASES,
    DatabaseSetupError,
    _download,
    _extract_functional_artifact,
    _expected_md5,
    _safe_extract,
    database_rows,
)


def test_database_catalog_distinguishes_available_and_planned(tmp_path):
    rows = {row["Database"]: row for row in database_rows(tmp_path)}

    assert rows["taxonomy"]["Status"] == "available"
    assert rows["taxonomy"]["Release"] == "2026-06-26"
    assert rows["functional"]["Status"] == "available"
    assert rows["functional"]["Release"] == "5.0.2"
    assert rows["amr"]["Status"] == "available"
    assert rows["virulence"]["Status"] == "planned"


def test_available_databases_record_upstream_provenance():
    for key in ("taxonomy", "functional", "amr"):
        spec = DATABASES[key]
        assert spec.release != "pending"
        assert spec.url
        assert spec.provenance_url


def test_database_subdirectory_is_diagnosed_instead_of_doubled(tmp_path):
    from metaquest.settings import load_config
    from metaquest.exceptions import ConfigError

    taxonomy = tmp_path / "databases/taxonomy"
    taxonomy.mkdir(parents=True)
    (taxonomy / "hash.k2d").write_text("index")
    with pytest.raises(ConfigError, match="root"):
        load_config(db_dir=taxonomy)


def test_amr_install_validates_before_replacing_existing_database(tmp_path, monkeypatch):
    from metaquest.core.screening import AMR_FILES, resolve_amr_database

    target = tmp_path / "amr"
    target.mkdir()
    (target / "previous").write_text("preserve")
    monkeypatch.setattr(database_manager.shutil, "which", lambda command: command)
    monkeypatch.setattr(database_manager.shutil, "disk_usage", lambda path: SimpleNamespace(free=10 * 1024**3))

    def failing_update(command, **kwargs):
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(database_manager, "run_command", failing_update)
    with pytest.raises(DatabaseSetupError, match="setup failed"):
        database_manager.install_database("amr", tmp_path, force=True)
    assert (target / "previous").read_text() == "preserve"

    def successful_update(command, **kwargs):
        if command[0] == "amrfinder_update":
            staged = Path(command[-1])
            version = staged / "test-v1"
            version.mkdir()
            for name in AMR_FILES:
                (version / name).write_text("test-v1\n")
            (staged / "latest").symlink_to("test-v1", target_is_directory=True)
        return SimpleNamespace(returncode=0, stdout="validated", stderr="")

    monkeypatch.setattr(database_manager, "run_command", successful_update)
    installed = database_manager.install_database("amr", tmp_path, force=True)
    assert resolve_amr_database(installed).name == "test-v1"
    assert database_manager.is_installed("amr", tmp_path)


def test_safe_extract_rejects_path_traversal(tmp_path):
    archive = tmp_path / "unsafe.tar.gz"
    payload = b"not allowed"
    with tarfile.open(archive, "w:gz") as bundle:
        member = tarfile.TarInfo("../outside.txt")
        member.size = len(payload)
        bundle.addfile(member, io.BytesIO(payload))

    destination = tmp_path / "extract"
    destination.mkdir()
    with pytest.raises(DatabaseSetupError, match="Unsafe path"):
        _safe_extract(archive, destination)

    assert not (tmp_path / "outside.txt").exists()


def test_archive_checksum_is_selected_by_filename(tmp_path):
    checksums = tmp_path / "checksums.md5"
    checksums.write_text(
        "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa  hash.k2d\n"
        "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb  database.tar.gz\n",
        encoding="utf-8",
    )

    assert _expected_md5(checksums, "database.tar.gz") == "b" * 32


def test_functional_gzip_extraction_streams_to_expected_name(tmp_path):
    archive = tmp_path / "eggnog.db.gz.partial"
    with gzip.open(archive, "wb") as handle:
        handle.write(b"SQLite format 3\x00test payload")
    destination = tmp_path / "output"
    destination.mkdir()

    _extract_functional_artifact(
        archive,
        destination,
        {"format": "gzip", "output": "eggnog.db"},
    )

    assert (destination / "eggnog.db").read_bytes() == b"SQLite format 3\x00test payload"


def test_functional_taxonomy_archive_uses_safe_extraction(tmp_path):
    source = tmp_path / "source.db"
    with sqlite3.connect(source) as connection:
        connection.execute("CREATE TABLE taxa (id INTEGER)")
    archive = tmp_path / "eggnog.taxa.tar.gz.partial"
    with tarfile.open(archive, "w:gz") as bundle:
        bundle.add(source, arcname="eggnog.taxa.db")
    destination = tmp_path / "output"
    destination.mkdir()

    _extract_functional_artifact(
        archive,
        destination,
        {"format": "tar", "output": "eggnog.taxa.db"},
    )

    assert (destination / "eggnog.taxa.db").is_file()


def test_pyrodigal_outputs_are_streamed_per_contig(tmp_path, monkeypatch):
    class FakePredictions:
        def __len__(self):
            return 1

        def write_translations(self, handle, sequence_id, **_kwargs):
            handle.write(f">{sequence_id}_1\nMPEPTIDE\n")

        def write_genes(self, handle, sequence_id, **_kwargs):
            handle.write(f">{sequence_id}_1\nATGCCCTAA\n")

        def write_gff(self, handle, sequence_id, **_kwargs):
            handle.write(
                f"{sequence_id}\tPyrodigal\tCDS\t1\t9\t.\t+\t0\tID={sequence_id}_1\n"
            )

    class FakeGeneFinder:
        def __init__(self, *, meta):
            assert meta is True

        def find_genes(self, sequence):
            assert isinstance(sequence, bytes)
            return FakePredictions()

    monkeypatch.setitem(
        sys.modules,
        "pyrodigal",
        SimpleNamespace(GeneFinder=FakeGeneFinder, __version__="test"),
    )

    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(
        ">long_contig\n" + "ATG" * 100 + "\n>short_contig\nATG\n",
        encoding="utf-8",
    )

    output_dir, count = run_gene_prediction(
        contigs,
        tmp_path / "results",
        min_contig_length=200,
    )

    assert count == 1
    assert ">contig_" in (output_dir / "genes.faa").read_text()
    assert "short_contig" not in (output_dir / "genes.faa").read_text()
    assert (output_dir / "genes.fna").exists()
    assert (output_dir / "genes.gff3").read_text().startswith("##gff-version 3")
    assert (output_dir / "contigs.stable.fasta").is_file()
    assert (output_dir / "summary.json").exists()
    assert "\tlong_contig\t" in (output_dir / "contig_id_map.tsv").read_text()


def test_eggnog_outputs_include_every_gene_and_are_restart_safe(tmp_path, monkeypatch):
    proteins = tmp_path / "genes.faa"
    proteins.write_text(
        ">gene_1\nMPEPTIDE\n>gene_2\nMSECOND\n>gene_3\nMTHIRD\n",
        encoding="utf-8",
    )
    database = tmp_path / "database"
    database.mkdir()
    for name in ("eggnog.db", "eggnog.taxa.db", "eggnog_proteins.dmnd"):
        (database / name).write_text("fixture\n", encoding="utf-8")

    emapper_runs = []

    run_options = []

    def fake_run(command, **kwargs):
        if command[-1] == "--version":
            return SimpleNamespace(returncode=0, stdout="emapper-2.1.15", stderr="")
        emapper_runs.append(command)
        run_options.append(kwargs)
        raw = tmp_path / "results" / "functional_annotation" / "metaquest.emapper.annotations"
        raw.write_text(
            "#query\tseed_ortholog\tevalue\tscore\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\tEC\tKEGG_ko\tPFAMs\n"
            "gene_1\t1.A\t1e-20\t100\tBacteria\tJ\tRibosomal protein\trplA\tGO:1,GO:2\t1.1.1.1\tko:K00001\tPF001\n"
            "gene_2\t2.B\t1e-10\t80\tBacteria\tK,L\tRegulator\tregA\t-\t-\tko:K00002\t-\n",
            encoding="utf-8",
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(functional_analysis.subprocess, "run", fake_run)
    monkeypatch.setattr(functional_analysis, "run_command", fake_run)

    table, categories, annotated, reused = run_functional_annotation(
        proteins,
        tmp_path / "results",
        database,
        threads=2,
    )

    rows = table.read_text(encoding="utf-8").splitlines()
    assert len(rows) == 4
    assert "gene_3\tunannotated" in rows[3]
    assert annotated == 2
    assert reused is False
    assert "COG\tJ\t1" in categories.read_text(encoding="utf-8")
    assert "KO\tko:K00001\t1" in categories.read_text(encoding="utf-8")
    assert len(emapper_runs) == 1
    assert Path(run_options[0]["cwd"]).parent == table.parent

    _, _, annotated_again, reused = run_functional_annotation(
        proteins,
        tmp_path / "results",
        database,
        threads=2,
    )

    assert annotated_again == 2
    assert reused is True
    assert len(emapper_runs) == 1


def test_functional_cache_hash_ignores_fasta_record_order(tmp_path):
    first = tmp_path / "first.faa"
    second = tmp_path / "second.faa"
    first.write_text(">gene_a\nMAAA\n>gene_b\nMBBB\n", encoding="utf-8")
    second.write_text(">gene_b\nMBBB\n>gene_a\nMAAA\n", encoding="utf-8")

    assert functional_analysis._canonical_fasta_sha256(
        first
    ) == functional_analysis._canonical_fasta_sha256(second)


class _FakeResponse(io.BytesIO):
    def __init__(self, payload, *, status, headers):
        super().__init__(payload)
        self.status = status
        self.headers = headers


class _InterruptedResponse(_FakeResponse):
    def __init__(self, payload, *, status, headers):
        super().__init__(payload, status=status, headers=headers)
        self._interrupted = False

    def read(self, size=-1):
        if self.tell() < len(self.getvalue()):
            return super().read(size)
        if not self._interrupted:
            self._interrupted = True
            raise http.client.IncompleteRead(b"")
        return b""


def test_download_continues_across_bounded_range_responses(tmp_path, monkeypatch):
    destination = tmp_path / "database.partial"
    responses = iter(
        [
            _FakeResponse(b"abc", status=200, headers={"Content-Length": "10"}),
            _FakeResponse(
                b"defg",
                status=206,
                headers={"Content-Range": "bytes 3-6/10"},
            ),
            _FakeResponse(
                b"hij",
                status=206,
                headers={"Content-Range": "bytes 7-9/10"},
            ),
        ]
    )
    ranges = []

    def fake_urlopen(request, **_kwargs):
        ranges.append(request.get_header("Range"))
        return next(responses)

    monkeypatch.setattr(database_manager.urllib.request, "urlopen", fake_urlopen)

    _download("https://example.test/database", destination, expected_bytes=10)

    assert destination.read_bytes() == b"abcdefghij"
    assert ranges == [None, "bytes=3-9", "bytes=7-9"]


def test_download_retries_after_interrupted_response(tmp_path, monkeypatch):
    destination = tmp_path / "database.partial"
    responses = iter(
        [
            _InterruptedResponse(
                b"abc",
                status=200,
                headers={"Content-Length": "6"},
            ),
            _FakeResponse(
                b"def",
                status=206,
                headers={"Content-Range": "bytes 3-5/6"},
            ),
        ]
    )
    monkeypatch.setattr(
        database_manager.urllib.request,
        "urlopen",
        lambda _request, **_kwargs: next(responses),
    )

    _download("https://example.test/database", destination, expected_bytes=6)

    assert destination.read_bytes() == b"abcdef"


def test_download_rejects_invalid_resume_range(tmp_path, monkeypatch):
    destination = tmp_path / "database.partial"
    destination.write_bytes(b"abc")
    response = _FakeResponse(
        b"def",
        status=206,
        headers={"Content-Range": "bytes 0-2/6"},
    )
    monkeypatch.setattr(
        database_manager.urllib.request,
        "urlopen",
        lambda _request, **_kwargs: response,
    )

    with pytest.raises(DatabaseSetupError, match="invalid byte range"):
        _download("https://example.test/database", destination, expected_bytes=6)

    assert destination.read_bytes() == b"abc"


def test_download_rejects_oversized_partial_without_request(tmp_path, monkeypatch):
    destination = tmp_path / "database.partial"
    destination.write_bytes(b"abcd")
    monkeypatch.setattr(
        database_manager.urllib.request,
        "urlopen",
        lambda _request, **_kwargs: pytest.fail("network request should not occur"),
    )

    with pytest.raises(DatabaseSetupError, match="larger than expected"):
        _download("https://example.test/database", destination, expected_bytes=3)


def test_download_stall_is_bounded_and_preserves_partial(tmp_path, monkeypatch):
    destination = tmp_path / "database.partial"

    def fake_urlopen(_request, **_kwargs):
        return _FakeResponse(
            b"",
            status=206,
            headers={"Content-Range": "bytes 0-2/3"},
        )

    monkeypatch.setattr(database_manager.urllib.request, "urlopen", fake_urlopen)

    with pytest.raises(DatabaseSetupError, match="did not complete after 2 requests"):
        _download(
            "https://example.test/database",
            destination,
            expected_bytes=3,
            max_requests=2,
        )

    assert destination.read_bytes() == b""
