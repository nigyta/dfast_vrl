"""Shared fixtures and helpers for DFAST_VRL / vadr2mss.py integration tests.

Tests are end-to-end: they spawn the CLI scripts as subprocesses and inspect
the resulting output directory. Tests that require a specific VADR model are
auto-skipped if the corresponding ``*.minfo`` file is missing under
``$VADRMODELDIR`` (default ``/vadr_models``). This keeps the suite portable
between the Docker image, the NIG super-computer environment, and developer
laptops.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
DFAST_VRL_BIN = REPO_ROOT / "dfast_vrl"
VADR2MSS_BIN = REPO_ROOT / "vadr2mss.py"
EXAMPLES_DIR = REPO_ROOT / "examples"
FLU_DIR = REPO_ROOT / "flu"  # symlink to /data/flu inside the container


def _candidate_pythons() -> list[str]:
    """Order of preference for the Python that DFAST_VRL scripts will run under.

    The scripts' shebang is ``#!/usr/bin/env python`` and they ``import Bio``,
    so we must hand them an interpreter where ``biopython`` is importable.
    pytest itself often runs from an isolated venv that lacks Biopython, so
    relying on ``PATH`` inherited from pytest is unreliable.
    """
    out: list[str] = []
    env_override = os.environ.get("DFV_PYTHON")
    if env_override:
        out.append(env_override)
    out.extend(
        [
            "/opt/conda/bin/python",
            "/opt/conda/bin/python3",
            "/miniconda3/bin/python",
        ]
    )
    which_python = shutil.which("python")
    if which_python:
        out.append(which_python)
    out.append(sys.executable)
    # de-dup while preserving order
    seen: set[str] = set()
    return [p for p in out if not (p in seen or seen.add(p))]


def _python_can_import_bio(python: str) -> bool:
    try:
        r = subprocess.run(
            [python, "-c", "import Bio"],
            capture_output=True,
            text=True,
            timeout=15,
        )
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False
    return r.returncode == 0


def find_dfv_python() -> str | None:
    for cand in _candidate_pythons():
        if Path(cand).exists() and _python_can_import_bio(cand):
            return cand
    return None


def vadr_model_dir() -> Path:
    return Path(os.environ.get("VADRMODELDIR", "/vadr_models"))


# Map of model-name (as used by vadr2mss_config.models) -> minfo path template.
# Versions are read from environment variables to mirror the Dockerfile
# defaults; if a variable is unset we fall back to a sensible default that
# matches the image.
_MODEL_MINFO = {
    "mpox": ("vadr-models-mpxv-{ver}/mpxv.minfo", "VADR_MPXV_MODELS_VERSION", "1.4.2-1"),
    "sarscov2": ("vadr-models-sarscov2-{ver}/sarscov2.minfo", "VADR_SCOV2_MODELS_VERSION", "1.3-2"),
    "RSV": ("vadr-models-rsv-{ver}/rsv.minfo", "VADR_RSV_MODELS_VERSION", "1.5-2"),
    "Flu": ("vadr-models-flu-1.6.3-2/flu.minfo", None, None),
}


def model_minfo_path(model: str) -> Path:
    template, env_var, default = _MODEL_MINFO[model]
    if env_var is None:
        rel = template
    else:
        version = os.environ.get(env_var, default)
        rel = template.format(ver=version)
    return vadr_model_dir() / rel


def require_model(model: str) -> None:
    """Skip the current test if the VADR model is not installed."""
    minfo = model_minfo_path(model)
    if not minfo.exists():
        pytest.skip(f"VADR model for {model!r} not found at {minfo}")


def require_path(path: Path, what: str) -> None:
    if not path.exists():
        pytest.skip(f"{what} not available at {path}")


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO_ROOT


@pytest.fixture(scope="session")
def dfv_python() -> str:
    """Path to a Python interpreter that has Biopython available."""
    py = find_dfv_python()
    if py is None:
        pytest.skip(
            "No Python interpreter with Biopython found. "
            "Set $DFV_PYTHON to override."
        )
    return py


@pytest.fixture(scope="session")
def dfast_vrl_bin() -> Path:
    assert DFAST_VRL_BIN.exists(), f"missing executable: {DFAST_VRL_BIN}"
    return DFAST_VRL_BIN


@pytest.fixture(scope="session")
def vadr2mss_bin() -> Path:
    assert VADR2MSS_BIN.exists(), f"missing executable: {VADR2MSS_BIN}"
    return VADR2MSS_BIN
