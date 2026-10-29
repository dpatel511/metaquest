"""Exercise cancellation using owned fixture processes, never real analyses."""
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import metaquest
import pytest


def alive(pid):
    status = Path(f"/proc/{pid}/status")
    if not status.exists():
        return False
    try:
        return not any(line.startswith("State:") and "Z" in line for line in status.read_text().splitlines())
    except FileNotFoundError:
        return False


@pytest.mark.skipif(sys.platform != "linux", reason="POSIX process-group integration")
@pytest.mark.parametrize("verbosity", ["standard", "debug"])
@pytest.mark.parametrize("interrupt", [signal.SIGINT, signal.SIGTERM])
def test_interrupt_reaps_tool_and_its_child(tmp_path, verbosity, interrupt):
    pid_file = tmp_path / "pids"
    child = tmp_path / "child.py"
    child.write_text(
        "import os, subprocess, sys, time\n"
        "from pathlib import Path\n"
        "grandchild = subprocess.Popen([sys.executable, '-c', 'import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(60)'])\n"
        "Path(sys.argv[1]).write_text(f'{os.getpid()} {grandchild.pid}')\n"
        "time.sleep(60)\n"
    )
    driver = tmp_path / "driver.py"
    driver.write_text(
        "import sys\nfrom metaquest.io import utils\nfrom metaquest.cli import main\n"
        "verbosity, child, pids = sys.argv[1:]\n"
        "def check(formatter, **kwargs):\n"
        "    formatter.run_subprocess([sys.executable, child, pids], 'fixture')\n"
        "utils.run_system_check = check\n"
        "sys.argv = ['metaquest', *(['--debug'] if verbosity == 'debug' else []), 'check']\n"
        "main()\n"
    )
    env = os.environ.copy()
    env["PYTHONPATH"] = str(Path(metaquest.__file__).resolve().parents[1])
    process = subprocess.Popen([sys.executable, str(driver), verbosity, str(child), str(pid_file)],
                               env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                               start_new_session=True)
    pids = []
    try:
        deadline = time.monotonic() + 10
        while not pid_file.exists() and time.monotonic() < deadline:
            time.sleep(0.05)
        assert pid_file.exists(), "fixture did not start"
        pids = [int(value) for value in pid_file.read_text().split()]
        os.kill(process.pid, interrupt)
        assert process.wait(timeout=10) == 128 + interrupt
        deadline = time.monotonic() + 3
        while any(alive(pid) for pid in pids) and time.monotonic() < deadline:
            time.sleep(0.05)
        assert not any(alive(pid) for pid in pids), "tool or descendant survived cancellation"
    finally:
        for pid in pids:
            if alive(pid):
                os.kill(pid, signal.SIGKILL)
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGKILL)
        process.wait(timeout=5)
