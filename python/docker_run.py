"""Run the Docker CLI as a child process with teardown on Ctrl+C (SIGINT)."""

from __future__ import annotations

import subprocess
import sys
from typing import Sequence

_KILL_TIMEOUT_S = 30.0


def run_docker(cmd: Sequence[str], *, kill_timeout_s: float = _KILL_TIMEOUT_S) -> int:
    """
    Run ``docker ...`` with inherited stdin/stdout/stderr.

    On KeyboardInterrupt, sends SIGTERM to the Docker client (TerminateProcess on
    Windows) so the daemon can stop the container; escalates to kill if needed.
    """
    proc = subprocess.Popen(list(cmd))
    try:
        return proc.wait()
    except KeyboardInterrupt:
        print("\nInterrupt received; stopping Docker...", file=sys.stderr)
        proc.terminate()
        try:
            proc.wait(timeout=kill_timeout_s)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        return 130
