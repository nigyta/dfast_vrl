"""Operational tests for the ``dfast_vrl`` CLI (SARS-CoV-2 pipeline)."""

from __future__ import annotations

from pathlib import Path

import pytest

from _helpers import (
    assert_files_exist,
    assert_non_empty,
    isolate_to_prefix,
    load_json,
    run_cli,
)
from conftest import EXAMPLES_DIR, require_model


SCOV2_FASTA = EXAMPLES_DIR / "nig_vrl_sample.fasta"
SCOV2_METADATA = EXAMPLES_DIR / "metadata.txt"
SCOV2_ISOLATE = "hCov-19/Japan/SZ-NIG-12345/2021"


def test_help(dfv_python: str, dfast_vrl_bin: Path) -> None:
    """``dfast_vrl --help`` prints usage and exits 0."""
    result = run_cli([dfv_python, str(dfast_vrl_bin), "--help"], timeout=30)
    assert result.returncode == 0, result.stderr
    assert "DFAST_VRL" in result.stdout
    assert "--input_fasta" in result.stdout


def test_version(dfv_python: str, dfast_vrl_bin: Path) -> None:
    """``dfast_vrl --version`` prints a version string."""
    result = run_cli([dfv_python, str(dfast_vrl_bin), "--version"], timeout=30)
    # argparse writes the version to stdout for some Python versions and
    # stderr for others; just look in both.
    assert result.returncode == 0, result.stderr
    assert "DFAST_VRL" in (result.stdout + result.stderr)


def test_missing_input_aborts(dfv_python: str, dfast_vrl_bin: Path, tmp_path: Path) -> None:
    """An input fasta path that does not exist should make the script fail."""
    out_dir = tmp_path / "out"
    result = run_cli(
        [
            dfv_python,
            str(dfast_vrl_bin),
            "-i",
            str(tmp_path / "does_not_exist.fasta"),
            "-o",
            str(out_dir),
        ],
        timeout=60,
    )
    assert result.returncode != 0


def test_full_pipeline_scov2(dfv_python: str, dfast_vrl_bin: Path, tmp_path: Path) -> None:
    """End-to-end run on the SARS-CoV-2 sample produces every expected artifact."""
    require_model("sarscov2")
    if not SCOV2_FASTA.exists():
        pytest.skip(f"SARS-CoV-2 sample fasta missing: {SCOV2_FASTA}")
    if not SCOV2_METADATA.exists():
        pytest.skip(f"SARS-CoV-2 metadata missing: {SCOV2_METADATA}")

    out_dir = tmp_path / "OUT_scov2"
    result = run_cli(
        [
            dfv_python,
            str(dfast_vrl_bin),
            "-i",
            str(SCOV2_FASTA),
            "-m",
            str(SCOV2_METADATA),
            "-o",
            str(out_dir),
        ],
    )
    assert result.returncode == 0, (
        f"dfast_vrl failed (rc={result.returncode})\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )

    # Top-level outputs.
    prefix = isolate_to_prefix(SCOV2_ISOLATE)
    assert_files_exist(
        out_dir,
        [
            "annotation.gbk",
            "annotation.gff",
            "genome.fna",
            "metadata.txt",
            "dfv_report.json",
            "variants.json",
            f"{prefix}.annt.tsv",
            f"{prefix}.seq.fa",
        ],
    )
    # Sub-working directories.
    assert (out_dir / "preprocessing").is_dir()
    assert (out_dir / "vadr").is_dir()
    assert (out_dir / "variants").is_dir()

    # Annotation file should contain at least one CDS feature.
    gbk_text = (out_dir / "annotation.gbk").read_text()
    assert "LOCUS" in gbk_text
    assert "CDS" in gbk_text

    # Report JSON should describe both the preprocessing and VADR stages.
    report = load_json(out_dir / "dfv_report.json")
    assert "preprocessing" in report
    assert "vadr" in report
    assert "warnings" in report
    vadr_block = report["vadr"]
    assert vadr_block.get("seq_status") in {"complete", "nearly complete", "draft"}
    assert vadr_block.get("number_of_sequence", 0) >= 1

    # Variant JSON should be valid JSON dict (may be empty for clean samples).
    variants = load_json(out_dir / "variants.json")
    assert isinstance(variants, dict)

    # MSS files (Mass Submission System) — must not be empty.
    assert_non_empty(out_dir / f"{prefix}.annt.tsv")
    assert_non_empty(out_dir / f"{prefix}.seq.fa")


def test_force_overwrite(dfv_python: str, dfast_vrl_bin: Path, tmp_path: Path) -> None:
    """Running into a pre-existing output directory without ``--force`` is refused."""
    if not SCOV2_FASTA.exists() or not SCOV2_METADATA.exists():
        pytest.skip("SARS-CoV-2 sample data missing")

    out_dir = tmp_path / "OUT_force"
    out_dir.mkdir()
    (out_dir / "stale.txt").write_text("pre-existing")

    result = run_cli(
        [
            dfv_python,
            str(dfast_vrl_bin),
            "-i",
            str(SCOV2_FASTA),
            "-m",
            str(SCOV2_METADATA),
            "-o",
            str(out_dir),
        ],
        timeout=60,
    )
    assert result.returncode != 0
    assert "already exists" in (result.stderr + result.stdout)
