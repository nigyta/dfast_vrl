"""Subprocess helpers shared between test modules."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Iterable


def run_cli(cmd: list[str], cwd: Path | None = None, timeout: int = 1800) -> subprocess.CompletedProcess:
    """Run a CLI command and return the CompletedProcess.

    Captures stdout/stderr as text. Does not raise on non-zero exit; tests
    decide what to assert. ``timeout`` is in seconds — VADR end-to-end runs
    can take several minutes for large genomes (e.g. mpox ~200kb).
    """
    return subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        capture_output=True,
        text=True,
        timeout=timeout,
    )


def assert_files_exist(out_dir: Path, names: Iterable[str]) -> None:
    missing = [n for n in names if not (out_dir / n).exists()]
    assert not missing, f"expected files missing under {out_dir}: {missing}"


def assert_non_empty(path: Path) -> None:
    assert path.exists(), f"missing: {path}"
    assert path.stat().st_size > 0, f"empty file: {path}"


def load_json(path: Path) -> dict:
    assert_non_empty(path)
    with path.open() as f:
        return json.load(f)


def isolate_to_prefix(isolate: str) -> str:
    """Mirror ``dfv.common.get_isolate``'s prefix derivation."""
    return isolate.replace("/", "_")
