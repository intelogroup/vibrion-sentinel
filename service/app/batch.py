"""Batch execution backends for the worker.

The worker runs pipeline commands through a BatchBackend. LocalBackend is
the current behavior (subprocess on the worker host). A future AWS
Batch / AWS HealthOmics backend implements the same three methods —
that is Jim's infra call, not this phase.

Timeouts are enforced by the backend: a runaway command is killed and
reported as a timeout, never left hanging.
"""

from __future__ import annotations

import abc
import subprocess
from typing import Optional


class InfraError(RuntimeError):
    """The execution infrastructure failed (not the pipeline).

    Download errors, R2 5xx, worker crashes, timeouts. These are the ONLY
    errors the worker retries (bounded: max 3, then failed).
    """


class JobTimeoutError(InfraError):
    """The command exceeded its wall-clock budget and was killed."""


class PipelineError(RuntimeError):
    """The pipeline itself failed (snakemake exit != 0).

    Never retried automatically — a failing pipeline will fail the same way
    on re-run. Human override: POST /jobs/{id}/retry.
    """


class BatchBackend(abc.ABC):
    @abc.abstractmethod
    def run(
        self, workdir: str, cmd: list[str], timeout_hours: float
    ) -> tuple[int, str]:
        """Run cmd in workdir. Returns (returncode, combined log).

        Raises JobTimeoutError if the wall clock is exceeded (the command
        is killed first). Raises InfraError for infrastructure failures.
        A nonzero returncode is NOT an exception — it is a PipelineError
        signal the caller converts (see worker).
        """


class LocalBackend(BatchBackend):
    """Run the command as a local subprocess with a wall-clock timeout."""

    def run(
        self, workdir: str, cmd: list[str], timeout_hours: float
    ) -> tuple[int, str]:
        header = f"$ {' '.join(cmd)}\n"
        try:
            proc = subprocess.run(
                cmd,
                cwd=workdir,
                capture_output=True,
                text=True,
                timeout=max(timeout_hours * 3600, 1),
            )
        except subprocess.TimeoutExpired as e:
            # The child is killed by subprocess on timeout.
            tail = ""
            if e.stdout:
                tail += e.stdout.decode() if isinstance(e.stdout, bytes) else e.stdout
            if e.stderr:
                tail += e.stderr.decode() if isinstance(e.stderr, bytes) else e.stderr
            raise JobTimeoutError(
                f"command exceeded {timeout_hours}h wall clock and was killed"
            ) from e
        except OSError as e:
            raise InfraError(f"failed to start command: {e}") from e
        log = f"{header}{proc.stdout}\n{proc.stderr}"
        return proc.returncode, log


def run_to_outcome(
    backend: BatchBackend, workdir: str, cmd: list[str], timeout_hours: float
) -> tuple[str, int, str]:
    """Run and classify: returns (outcome, returncode, log) where outcome is
    'ok' | 'pipeline_error' | 'timeout' | 'infra_error'."""
    try:
        rc, log = backend.run(workdir, cmd, timeout_hours)
    except JobTimeoutError as e:
        return "timeout", -1, str(e)
    except InfraError as e:
        return "infra_error", -1, str(e)
    if rc != 0:
        return "pipeline_error", rc, log
    return "ok", rc, log
