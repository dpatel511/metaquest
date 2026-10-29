"""Own native process groups for interruption cleanup."""
from contextlib import contextmanager
import os
import signal
import subprocess


class TerminationRequested(KeyboardInterrupt):
    """Allow SIGTERM to unwind the same cleanup path as Ctrl+C."""


def handle_sigterm(_signum, _frame):
    raise TerminationRequested("SIGTERM received")


@contextmanager
def managed_process(command, **kwargs):
    process = subprocess.Popen(command, start_new_session=os.name == "posix", **kwargs)
    try:
        yield process
    finally:
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGTERM)
            except ProcessLookupError:
                pass
        elif process.poll() is None:
            process.terminate()
        try:
            process.wait(timeout=3)
        except subprocess.TimeoutExpired:
            process.kill()
        if os.name == "posix":
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        process.wait()
        for stream in (process.stdout, process.stderr):
            if stream is not None:
                stream.close()


def run_command(command, *, capture_output=False, **kwargs):
    if capture_output:
        kwargs.update(stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    with managed_process(command, **kwargs) as process:
        stdout, stderr = process.communicate()
        return subprocess.CompletedProcess(command, process.returncode, stdout, stderr)
