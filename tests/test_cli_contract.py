"""Observable CLI behaviour, exercised outside the working directory."""
import io
import os
from pathlib import Path
import subprocess
import sys

import metaquest
from metaquest.config import __version__
from metaquest.io.output_formatter import OutputFormatter, Spinner
import pytest
import yaml


def invoke(tmp_path, *arguments, env_updates=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(metaquest.__file__).resolve().parents[1])
    env.update(env_updates or {})
    return subprocess.run(
        [sys.executable, "-c", "from metaquest.cli import main; main()", *arguments],
        cwd=tmp_path, env=env, text=True, capture_output=True, timeout=30,
    )


@pytest.mark.parametrize("command", [
    [], ["run"], ["validate"], ["check"], ["compare"], ["databases"], ["init-config"],
])
@pytest.mark.parametrize("option", ["-h", "--help"])
def test_help_is_stdout(tmp_path, command, option):
    result = invoke(tmp_path, *command, option)
    assert result.returncode == 0
    assert "usage:" in result.stdout.lower()
    assert result.stderr == ""


@pytest.mark.parametrize("option", ["-v", "--version"])
def test_version_is_stdout(tmp_path, option):
    result = invoke(tmp_path, option)
    assert result.returncode == 0
    assert result.stdout.strip() == f"MetaQuest version {__version__}"
    assert result.stderr == ""


@pytest.mark.parametrize("arguments, code, message", [
    ([], 2, "required"),
    (["--unknown-option"], 2, "error:"),
    (["run"], 2, "required"),
    (["compare"], 2, "--sample-sheet"),
    (["validate", "--single", "missing.fastq"], 1, "Input file not found"),
    (["--quiet", "validate", "--single", "missing.fastq"], 1, "Input file not found"),
    (["databases"], 1, "No database selected"),
    (["check", "--config", "missing.yaml"], 1, "Config file not found"),
])
def test_failures_use_stderr(tmp_path, arguments, code, message):
    result = invoke(tmp_path, *arguments)
    assert result.returncode == code
    assert result.stdout == ""
    assert message in result.stderr
    assert "Traceback" not in result.stderr
    assert "\x1b" not in result.stderr


@pytest.mark.parametrize("command, option, value", [
    ("validate", "--min-quality", "-1"),
    ("validate", "--min-sequences", "0"),
    ("validate", "--overrep-threshold", "-0.1"),
    ("validate", "--overrep-threshold", "1.1"),
    ("validate", "--overrep-threshold", "nan"),
    ("validate", "--overrep-threshold", "inf"),
    ("run", "--annotation-threads", "0"),
])
def test_invalid_thresholds_fail_before_output_creation(tmp_path, command, option, value):
    result = invoke(tmp_path, command, "--single", "missing.fastq", option, value)
    assert result.returncode == 2
    assert option in result.stderr
    assert result.stdout == ""
    assert not any(tmp_path.iterdir())


@pytest.mark.parametrize("quiet", [False, True])
@pytest.mark.parametrize("content, message", [
    ("", "File is empty"),
    ("not a FASTQ record\n", "Invalid FASTQ format"),
    ("@one\nACGT\n+\nIIII\n", "Insufficient sequences"),
])
def test_invalid_fastq_errors_survive_quiet(tmp_path, quiet, content, message):
    reads = tmp_path / "invalid reads.fastq"
    reads.write_text(content)
    options = ["--quiet"] if quiet else []
    result = invoke(tmp_path, *options, "validate", "--single", str(reads))
    assert result.returncode == 1
    assert result.stdout == ""
    assert message in result.stderr
    assert "All files validated successfully" not in result.stderr


@pytest.mark.parametrize("quiet", [False, True])
def test_successful_validation_with_space_in_path(tmp_path, quiet):
    reads = tmp_path / "valid reads.fastq"
    reads.write_text("@one\nACGT\n+\nIIII\n")
    options = ["--quiet"] if quiet else ["--no-color"]
    result = invoke(tmp_path, *options, "validate", "--single", str(reads), "--min-sequences", "1")
    assert result.returncode == 0
    assert result.stdout == ""
    assert "\x1b" not in result.stderr
    assert "\r" not in result.stderr
    if quiet:
        assert result.stderr == ""
    else:
        assert "All files validated successfully" in result.stderr


def test_missing_tools_fail_in_preflight(tmp_path):
    reads = tmp_path / "reads.fastq"
    reads.write_text("@one\nACGT\n+\nIIII\n")
    output = tmp_path / "new results"
    result = invoke(
        tmp_path, "run", "--single", str(reads), "--taxonomy-only",
        "--output", str(output), "--db-dir", str(tmp_path / "empty-db"),
        env_updates={"PATH": str(tmp_path / "empty-bin")},
    )
    assert result.returncode == 1
    assert result.stdout == ""
    assert "Tool not found: 'fastp'" in result.stderr
    assert "Tools and databases are ready" not in result.stderr
    assert {path.name for path in output.iterdir()} == {"metaquest.log"}


def test_packaged_config_can_be_written_outside_checkout(tmp_path):
    config = tmp_path / "config with spaces.yaml"
    result = invoke(tmp_path, "init-config", "--output", str(config))
    assert result.returncode == 0
    assert result.stdout == ""
    assert "Configuration file created" in result.stderr
    assert {"databases", "assembly", "annotation"} <= yaml.safe_load(config.read_text()).keys()


def test_debug_tool_output_stays_on_stderr_and_in_log(tmp_path, capsys):
    log = tmp_path / "metaquest.log"
    formatter = OutputFormatter("debug", log_file=log)
    command = [sys.executable, "-c", "import sys; print('tool data'); print('tool warning', file=sys.stderr)"]
    code, stdout, stderr = formatter.run_subprocess(command, "fixture")
    formatter.log_handle.close()
    captured = capsys.readouterr()
    assert (code, stdout, stderr) == (0, "tool data\n", "tool warning\n")
    assert captured.out == ""
    assert "tool data" in captured.err and "tool warning" in captured.err
    assert "tool data" in log.read_text() and "tool warning" in log.read_text()


def test_redirected_stderr_disables_animation_even_when_stdout_is_terminal(monkeypatch):
    stdout = io.StringIO()
    stderr = io.StringIO()
    monkeypatch.setattr(stdout, "isatty", lambda: True)
    monkeypatch.setattr(sys, "stdout", stdout)
    monkeypatch.setattr(sys, "stderr", stderr)
    formatter = OutputFormatter()
    with formatter.spinner("fixture"):
        pass
    assert stdout.getvalue() == ""
    assert stderr.getvalue() == ""


def test_spinner_uses_diagnostic_stream(capsys):
    spinner = Spinner("fixture")
    spinner.start()
    spinner.stop("finished")
    captured = capsys.readouterr()
    assert captured.out == ""
    assert "finished" in captured.err
