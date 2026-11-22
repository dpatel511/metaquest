"""An unrelated caller repository must not become an eggNOG build identity."""
import os
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

from metaquest.core import functional_analysis


@pytest.mark.skipif(shutil.which("git") is None, reason="Git not installed")
@pytest.mark.parametrize("scratch_in_repository", [False, True])
@pytest.mark.parametrize("inherited_git_environment", [False, True])
def test_emapper_version_ignores_callers_git_repository(
    tmp_path, monkeypatch, scratch_in_repository, inherited_git_environment,
):
    caller = tmp_path / "unrelated-project"
    caller.mkdir()
    environment = {key: value for key, value in os.environ.items() if not key.startswith("GIT_")}
    subprocess.run(["git", "init", "--quiet", str(caller)], check=True, env=environment)
    subprocess.run(["git", "-C", str(caller), "-c", "user.name=Fixture",
                    "-c", "user.email=fixture@example.invalid", "-c", "commit.gpgsign=false",
                    "-c", f"core.hooksPath={tmp_path / 'no-hooks'}", "commit", "--quiet",
                    "--allow-empty", "-m", "fixture"], check=True, env=environment)
    probe = tmp_path / "emapper_probe.py"
    probe.write_text(
        "import subprocess\n"
        "result = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], "
        "capture_output=True, text=True)\n"
        "suffix = '-' + result.stdout.strip() if result.returncode == 0 else ''\n"
        "print('emapper-2.1.15' + suffix)\n",
        encoding="utf-8",
    )
    real_run = subprocess.run

    def run_probe(command, **kwargs):
        assert command[0] == "emapper.py"
        return real_run([sys.executable, str(probe)], **kwargs)

    monkeypatch.setattr(functional_analysis.subprocess, "run", run_probe)
    monkeypatch.chdir(caller)
    if inherited_git_environment:
        monkeypatch.setenv("GIT_DIR", str(caller / ".git"))
        monkeypatch.setenv("GIT_WORK_TREE", str(caller))
        monkeypatch.setenv("GIT_COMMON_DIR", str(caller / ".git"))
    if scratch_in_repository:
        monkeypatch.setattr(functional_analysis.tempfile, "tempdir", str(caller))

    assert functional_analysis._emapper_version(tmp_path / "database") == "emapper-2.1.15"
    assert Path.cwd() == caller
