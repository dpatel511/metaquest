"""Versioned reference-database installation outside the source tree."""

from __future__ import annotations

import gzip
import hashlib
import http.client
import json
import re
import shutil
import sqlite3
import tarfile
import tempfile
import urllib.error
import urllib.request
import zlib
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

from metaquest.exceptions import MetaQuestError
from metaquest.io.processes import run_command


class DatabaseSetupError(MetaQuestError):
    """A reference database could not be downloaded or validated."""


@dataclass(frozen=True)
class DatabaseSpec:
    key: str
    name: str
    release: str
    status: str
    download_gb: float | None
    installed_gb: float | None
    url: str | None = None
    checksum_url: str | None = None
    provenance_url: str | None = None
    terms_url: str | None = None


DATABASES = {
    "taxonomy": DatabaseSpec(
        key="taxonomy",
        name="Kraken2 Standard-8 with Bracken",
        release="2026-06-26",
        status="available",
        download_gb=5.5,
        installed_gb=7.5,
        url=(
            "https://genome-idx.s3.amazonaws.com/kraken/"
            "k2_standard_08_GB_20260626.tar.gz"
        ),
        checksum_url=(
            "https://genome-idx.s3.amazonaws.com/kraken/"
            "standard_08_GB_20260626/standard_08_GB.md5"
        ),
        provenance_url="https://benlangmead.github.io/aws-indexes/k2",
        terms_url="https://www.ncbi.nlm.nih.gov/home/about/policies/",
    ),
    "functional": DatabaseSpec(
        key="functional",
        name="eggNOG-mapper core functional database",
        release="5.0.2",
        status="available",
        download_gb=11.3,
        installed_gb=53.0,
        url="http://eggnog5.embl.de/download/emapperdb-5.0.2",
        provenance_url="http://eggnog5.embl.de/",
    ),
    "amr": DatabaseSpec(
        "amr", "AMRFinderPlus core + selected virulence/stress data",
        "tool-compatible release recorded at installation", "available", None, None,
        url="https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/",
        provenance_url="https://github.com/ncbi/amr/wiki/AMRFinderPlus-database",
    ),
    "virulence": DatabaseSpec(
        "virulence", "VFDB Set A", "pending", "planned", None, None,
    ),
}


FUNCTIONAL_ARTIFACTS = (
    {
        "archive": "eggnog.db.gz",
        "output": "eggnog.db",
        "compressed_bytes": 6_776_977_123,
        "format": "gzip",
    },
    {
        "archive": "eggnog.taxa.tar.gz",
        "output": "eggnog.taxa.db",
        "compressed_bytes": 72_797_584,
        "format": "tar",
    },
    {
        "archive": "eggnog_proteins.dmnd.gz",
        "output": "eggnog_proteins.dmnd",
        "compressed_bytes": 5_208_806_170,
        "format": "gzip",
    },
)


def database_rows(db_root: Path) -> list[dict[str, str]]:
    """Return display-ready database status rows."""
    rows = []
    for spec in DATABASES.values():
        installed = is_installed(spec.key, db_root)
        size = (
            f"{spec.download_gb:g} GB download / {spec.installed_gb:g} GB installed"
            if spec.download_gb is not None
            else "upstream release size" if spec.status == "available" else "not finalized"
        )
        rows.append(
            {
                "Database": spec.key,
                "Description": spec.name,
                "Release": spec.release,
                "Size": size,
                "Status": "installed" if installed else spec.status,
            }
        )
    return rows


def is_installed(key: str, db_root: Path) -> bool:
    """Check the manifest and essential files for an installed database."""
    target = Path(db_root) / key
    manifest = target / "metaquest-db.json"
    if not manifest.exists():
        return False
    try:
        recorded = json.loads(manifest.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    if (not isinstance(recorded, dict) or recorded.get("key") != key
            or not isinstance(recorded.get("release"), str) or not recorded["release"].strip()):
        return False
    if key == "taxonomy":
        return all((target / name).is_file() and (target / name).stat().st_size > 0
                   for name in ("hash.k2d", "opts.k2d", "taxo.k2d"))
    if key == "functional":
        return all(
            (target / name).is_file() and (target / name).stat().st_size > 0
            for name in (
                "eggnog.db",
                "eggnog.taxa.db",
                "eggnog_proteins.dmnd",
            )
        )
    if key == "amr":
        from metaquest.core.screening import ScreeningError, resolve_amr_database

        try:
            resolve_amr_database(target)
        except ScreeningError:
            return False
        return True
    return False


def _download(
    url: str,
    destination: Path,
    *,
    expected_bytes: int | None = None,
    max_requests: int = 64,
) -> None:
    """Download a file, continuing bounded HTTP range requests until complete."""
    target_size = expected_bytes
    last_error: Exception | None = None

    for _ in range(max_requests):
        existing = destination.stat().st_size if destination.exists() else 0
        if target_size is not None:
            if existing == target_size:
                return
            if existing > target_size:
                raise DatabaseSetupError(
                    f"Partial file is larger than expected for {destination.name}: "
                    f"expected {target_size} bytes, got {existing}"
                )

        range_end = target_size - 1 if target_size is not None else ""
        request = urllib.request.Request(
            url,
            headers={"Range": f"bytes={existing}-{range_end}"} if existing else {},
        )
        before_request = existing
        range_limit = None
        try:
            response = urllib.request.urlopen(request, timeout=30)
            status = getattr(response, "status", 200) or 200
            if status == 206:
                content_range = response.headers.get("Content-Range", "")
                match = re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)", content_range)
                if not match:
                    response.close()
                    raise DatabaseSetupError(
                        f"Upstream returned an invalid byte range for {destination.name}"
                    )
                start, end, total = (int(value) for value in match.groups())
                if start != existing or end < start or end >= total:
                    response.close()
                    raise DatabaseSetupError(
                        f"Upstream returned an invalid byte range for {destination.name}"
                    )
                if target_size is None:
                    target_size = total
                elif total != target_size:
                    response.close()
                    raise DatabaseSetupError(
                        f"Upstream size changed for {destination.name}: "
                        f"expected {target_size} bytes, got {total}"
                    )
                mode = "ab"
                range_limit = end - start + 1
            elif existing:
                response.close()
                raise DatabaseSetupError(
                    f"Upstream did not honor resume request for {destination.name}"
                )
            elif status != 200:
                response.close()
                raise DatabaseSetupError(
                    f"Upstream returned HTTP {status} for {destination.name}"
                )
            else:
                mode = "wb"
                if target_size is None:
                    content_length = response.headers.get("Content-Length")
                    if content_length and content_length.isdigit():
                        target_size = int(content_length)

            with response, destination.open(mode) as output:
                shutil.copyfileobj(response, output, length=1024 * 1024)
            last_error = None
        except urllib.error.HTTPError as exc:
            complete = re.fullmatch(r"bytes \*/(\d+)", exc.headers.get("Content-Range", ""))
            exc.close()
            if exc.code == 416 and complete and int(complete.group(1)) == existing:
                if target_size is None or target_size == existing:
                    return
            if 400 <= exc.code < 500 and exc.code not in (408, 429):
                raise DatabaseSetupError(f"HTTP {exc.code} downloading {destination.name}") from exc
            last_error = exc
        except (
            ConnectionError,
            TimeoutError,
            http.client.HTTPException,
            urllib.error.URLError,
        ) as exc:
            last_error = exc

        observed = destination.stat().st_size if destination.exists() else 0
        if range_limit is not None and observed - before_request > range_limit:
            with destination.open("r+b") as output:
                output.truncate(before_request)
            raise DatabaseSetupError(f"Response body exceeds advertised byte range for {destination.name}")
        if target_size is not None and observed > target_size:
            raise DatabaseSetupError(
                f"Downloaded more data than expected for {destination.name}: "
                f"expected {target_size} bytes, got {observed}"
            )
        if target_size is None:
            if last_error is not None:
                continue
            return
        if observed == target_size:
            return
        if observed <= before_request and last_error is None:
            last_error = DatabaseSetupError(
                f"Download made no progress for {destination.name}"
            )

    detail = f": {last_error}" if last_error else ""
    raise DatabaseSetupError(
        f"Download did not complete after {max_requests} requests for "
        f"{destination.name}{detail}"
    )


def _md5(path: Path) -> str:
    digest = hashlib.md5()  # nosec: upstream publishes MD5 for transfer integrity
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _expected_md5(checksum_file: Path, archive_name: str) -> str:
    """Read the checksum belonging to the archive, not an extracted DB file."""
    for line in checksum_file.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2 and Path(fields[-1]).name == archive_name:
            return fields[0].lower()
    raise DatabaseSetupError(
        f"Upstream checksum file has no entry for {archive_name}"
    )


def _safe_extract(archive: Path, destination: Path) -> None:
    destination_resolved = destination.resolve()
    with tarfile.open(archive, "r:gz") as bundle:
        for member in bundle.getmembers():
            member_path = (destination / member.name).resolve()
            if (
                not (member.isfile() or member.isdir())
                or member.issym()
                or member.islnk()
                or destination_resolved not in member_path.parents
            ):
                raise DatabaseSetupError(
                    f"Unsafe path in database archive: {member.name}"
                )
        bundle.extractall(destination)


def _promote_database(staged: Path, target: Path) -> None:
    """Retain the previous installation until promotion succeeds, including Ctrl+C."""
    previous = target.with_name(f".{target.name}.previous")
    if previous.exists():
        raise DatabaseSetupError(f"Unresolved database backup: {previous}")
    try:
        if target.exists():
            target.replace(previous)
        staged.replace(target)
    except BaseException:
        if previous.exists() and not target.exists():
            previous.replace(target)
        raise
    if previous.exists():
        shutil.rmtree(previous)


def _validate_sqlite(path: Path) -> None:
    """Perform a cheap read-only validation without scanning a large database."""
    with path.open("rb") as handle:
        if handle.read(16) != b"SQLite format 3\x00":
            raise DatabaseSetupError(f"Invalid SQLite database header: {path.name}")
    try:
        connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
        try:
            connection.execute("SELECT name FROM sqlite_master LIMIT 1").fetchone()
        finally:
            connection.close()
    except sqlite3.DatabaseError as exc:
        raise DatabaseSetupError(f"Invalid SQLite database: {path.name}") from exc


def _validate_functional_tools(target: Path, release: str) -> None:
    diamond = shutil.which("diamond")
    emapper = shutil.which("emapper.py")
    if not diamond or not emapper:
        raise DatabaseSetupError(
            "DIAMOND and eggNOG-mapper are required to validate functional data"
        )

    commands = (
        ([diamond, "dbinfo", "--db", str(target / "eggnog_proteins.dmnd")], "DIAMOND"),
        ([emapper, "--data_dir", str(target), "--version"], "eggNOG-mapper"),
    )
    for command, name in commands:
        result = run_command(command, capture_output=True, text=True)
        output = f"{result.stdout}\n{result.stderr}"
        if result.returncode or (name == "eggNOG-mapper" and release not in output):
            raise DatabaseSetupError(f"{name} could not validate the functional database")


def _extract_functional_artifact(
    archive: Path,
    destination: Path,
    artifact: dict[str, object],
) -> None:
    if artifact["format"] == "tar":
        _safe_extract(archive, destination)
        return
    output = destination / str(artifact["output"])
    try:
        with gzip.open(archive, "rb") as source, output.open("wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
    except (EOFError, gzip.BadGzipFile, zlib.error) as exc:
        raise DatabaseSetupError(f"Invalid gzip archive: {archive.name}") from exc


def _install_functional(
    spec: DatabaseSpec,
    db_root: Path,
    target: Path,
    *,
    force: bool,
    notify: Callable[[str], None] | None,
) -> Path:
    downloads = db_root / ".downloads" / "functional"
    downloads.mkdir(parents=True, exist_ok=True)
    provenance = []

    for artifact in FUNCTIONAL_ARTIFACTS:
        archive_name = str(artifact["archive"])
        partial = downloads / f"{archive_name}.partial"
        url = f"{spec.url}/{archive_name}"
        if notify:
            if partial.exists() and partial.stat().st_size:
                notify(
                    f"Downloading {archive_name} from existing partial file "
                    f"({partial.stat().st_size / 1024**3:.1f} GB)"
                )
            else:
                notify(f"Downloading {archive_name}")
        expected_bytes = int(artifact["compressed_bytes"])
        if not partial.exists() or partial.stat().st_size != expected_bytes:
            _download(url, partial, expected_bytes=expected_bytes)
        if partial.stat().st_size != expected_bytes:
            raise DatabaseSetupError(
                f"Unexpected size for {archive_name}: expected {expected_bytes} bytes, "
                f"got {partial.stat().st_size}"
            )
        if notify:
            notify(f"Calculating SHA-256 for {archive_name}")
        provenance.append(
            {
                "url": url,
                "compressed_bytes": expected_bytes,
                "sha256": _sha256(partial),
            }
        )

    with tempfile.TemporaryDirectory(prefix=".metaquest-functional-", dir=db_root) as tmp:
        extracted = Path(tmp) / "functional"
        extracted.mkdir()
        for artifact in FUNCTIONAL_ARTIFACTS:
            archive_name = str(artifact["archive"])
            if notify:
                notify(f"Extracting {archive_name}")
            partial = downloads / f"{archive_name}.partial"
            try:
                _extract_functional_artifact(partial, extracted, artifact)
            except (DatabaseSetupError, tarfile.TarError, EOFError):
                partial.unlink(missing_ok=True)
                raise

        required = tuple(str(item["output"]) for item in FUNCTIONAL_ARTIFACTS)
        for name in required:
            path = extracted / name
            if not path.is_file() or path.stat().st_size == 0:
                raise DatabaseSetupError(f"Functional database is missing required file: {name}")
        _validate_sqlite(extracted / "eggnog.db")
        _validate_sqlite(extracted / "eggnog.taxa.db")
        _validate_functional_tools(extracted, spec.release)

        manifest = {
            **asdict(spec),
            "installed_at": datetime.now(timezone.utc).isoformat(),
            "artifacts": provenance,
            "source": spec.provenance_url,
        }
        (extracted / "metaquest-db.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

        _promote_database(extracted, target)

    for artifact in FUNCTIONAL_ARTIFACTS:
        (downloads / f"{artifact['archive']}.partial").unlink()
    try:
        downloads.rmdir()
        downloads.parent.rmdir()
    except OSError:
        pass
    if notify:
        notify("Functional database manifest written and installation validated")
    return target


def _install_amr(db_root: Path, target: Path, notify) -> Path:
    """Use NCBI's updater, retaining its versioned layout and native thresholds."""
    from metaquest.core.screening import AMR_FILES, resolve_amr_database

    updater = shutil.which("amrfinder_update")
    tool = shutil.which("amrfinder")
    if not updater or not tool:
        raise DatabaseSetupError(
            "AMRFinderPlus is required; install environment/screening.yml separately and add its bin directory to PATH (see docs/installation.md)"
        )
    if shutil.disk_usage(db_root).free < 2 * 1024**3:
        raise DatabaseSetupError("AMRFinderPlus setup requires at least 2 GiB of free working space")
    with tempfile.TemporaryDirectory(prefix=".metaquest-amr-", dir=db_root) as temporary:
        staged = Path(temporary) / "amr"
        staged.mkdir()
        if notify:
            notify("Downloading and indexing the tool-compatible NCBI AMRFinderPlus database")
        with (staged / "setup.log").open("w", encoding="utf-8") as log:
            result = run_command([updater, "--database", str(staged)], stdout=log, stderr=log)
        if result.returncode:
            tail = (staged / "setup.log").read_text(errors="replace")[-4000:]
            raise DatabaseSetupError(f"AMRFinderPlus database setup failed: {tail}")
        database = resolve_amr_database(staged)
        check = run_command([tool, "--database", str(database), "--database_version"],
                               capture_output=True, text=True)
        if check.returncode:
            raise DatabaseSetupError(f"AMRFinderPlus database validation failed: {check.stderr}")
        manifest = {
            **asdict(DATABASES["amr"]), "installed_at": datetime.now(timezone.utc).isoformat(),
            "release": (database / "version.txt").read_text().strip(),
            "validation": (check.stdout + check.stderr).strip(),
            "sha256": {name: _sha256(database / name) for name in AMR_FILES},
        }
        (staged / "metaquest-db.json").write_text(json.dumps(manifest, indent=2) + "\n")
        _promote_database(staged, target)
    return target


def install_database(
    key: str,
    db_root: Path,
    *,
    force: bool = False,
    notify: Callable[[str], None] | None = None,
) -> Path:
    """Download, verify, extract, validate, and atomically install a database."""
    if key not in DATABASES:
        raise DatabaseSetupError(f"Unknown database: {key}")
    spec = DATABASES[key]
    if spec.status != "available" or not spec.url:
        raise DatabaseSetupError(
            f"Database '{key}' is planned but its installation is not finalized"
        )

    db_root = Path(db_root).expanduser().resolve()
    target = db_root / key
    previous = target.with_name(f".{target.name}.previous")
    if previous.exists():
        if not target.exists():
            previous.replace(target)
        elif is_installed(key, db_root):
            shutil.rmtree(previous)
        else:
            raise DatabaseSetupError(f"Incomplete database replacement; previous installation retained at {previous}")
    if is_installed(key, db_root) and not force:
        if notify:
            notify(f"{spec.name} is already installed and valid")
        return target
    if target.exists() and not force:
        raise DatabaseSetupError(
            f"Database directory exists but is not valid: {target}",
            suggestions=["Move it aside or rerun with --force to replace it."],
        )

    db_root.mkdir(parents=True, exist_ok=True)
    if key == "amr":
        return _install_amr(db_root, target, notify)
    required_bytes = int((spec.download_gb + spec.installed_gb + 1) * 1024**3)
    if notify:
        notify(
            f"Disk preflight: {spec.download_gb:g} GB download, "
            f"{spec.installed_gb:g} GB installed"
        )
    if shutil.disk_usage(db_root).free < required_bytes:
        raise DatabaseSetupError(
            f"Insufficient free space for {spec.name}; "
            f"at least {required_bytes / 1024**3:.1f} GB is required"
        )

    if key == "functional":
        return _install_functional(
            spec,
            db_root,
            target,
            force=force,
            notify=notify,
        )
    if not spec.checksum_url:
        raise DatabaseSetupError(f"No checksum source is configured for {spec.name}")

    downloads = db_root / ".downloads"
    downloads.mkdir(exist_ok=True)
    archive_name = Path(spec.url).name
    archive = downloads / f"{archive_name}.partial"
    checksum_file = downloads / f"{archive_name}.md5"

    if notify and archive.exists():
        notify(
            f"Resuming database archive at "
            f"{archive.stat().st_size / 1024**3:.1f} GB"
        )

    if notify:
        notify("Downloading upstream checksum")
    checksum_file.unlink(missing_ok=True)
    _download(spec.checksum_url, checksum_file)
    if notify:
        notify("Downloading database archive")
    _download(spec.url, archive)
    if notify:
        notify("Verifying archive checksum")
    expected = _expected_md5(checksum_file, archive_name)
    observed = _md5(archive)
    if observed != expected:
        archive.unlink(missing_ok=True)
        raise DatabaseSetupError(
            f"Checksum mismatch for {spec.name}: expected {expected}, got {observed}"
        )

    with tempfile.TemporaryDirectory(prefix=".metaquest-db-", dir=db_root) as tmp:
        work = Path(tmp)
        extracted = work / "extracted"
        extracted.mkdir()

        if notify:
            notify("Extracting database archive")
        _safe_extract(archive, extracted)
        required = ("hash.k2d", "opts.k2d", "taxo.k2d")
        if not all((extracted / name).is_file() and (extracted / name).stat().st_size for name in required):
            raise DatabaseSetupError(
                f"{spec.name} archive is missing required Kraken2 index files"
            )

        manifest = {
            **asdict(spec),
            "installed_at": datetime.now(timezone.utc).isoformat(),
            "archive_md5": observed,
            "source": spec.provenance_url,
        }
        (extracted / "metaquest-db.json").write_text(
            json.dumps(manifest, indent=2) + "\n",
            encoding="utf-8",
        )

        _promote_database(extracted, target)
    archive.unlink(missing_ok=True)
    checksum_file.unlink(missing_ok=True)
    try:
        downloads.rmdir()
    except OSError:
        pass
    if notify:
        notify("Database manifest written and installation validated")
    return target
