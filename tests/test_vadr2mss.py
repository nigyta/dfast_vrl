"""Operational tests for ``vadr2mss.py`` across mpox / sarscov2 / RSV / Flu."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from _helpers import (
    assert_files_exist,
    assert_non_empty,
    isolate_to_prefix,
    load_json,
    run_cli,
)
from conftest import EXAMPLES_DIR, FLU_DIR, require_model, require_path


@dataclass(frozen=True)
class ModelCase:
    """Describes one virus / model combination to run end-to-end."""

    id: str             # pytest parametrize id
    model: str          # value passed via -M
    fasta: Path
    metadata: Path
    isolate: str        # used to derive the MSS file prefix


CASES: list[ModelCase] = [
    ModelCase(
        id="mpox",
        model="mpox",
        fasta=EXAMPLES_DIR / "mpox" / "LC756923.fasta",
        metadata=EXAMPLES_DIR / "mpox" / "metadata_mpox.txt",
        isolate="MPXV/Japan/NIG-xxxxx/2023",
    ),
    ModelCase(
        id="sarscov2",
        model="sarscov2",
        fasta=EXAMPLES_DIR / "MW850590.fasta",
        metadata=EXAMPLES_DIR / "metadata.txt",
        isolate="hCov-19/Japan/SZ-NIG-12345/2021",
    ),
    ModelCase(
        id="rsv",
        model="RSV",
        fasta=EXAMPLES_DIR / "rsv" / "PX986933.fasta",
        metadata=EXAMPLES_DIR / "rsv" / "metadata_rsv.tsv",
        isolate="HRSV/Human/B/Japan/XX-XXX-001/2026",
    ),
    ModelCase(
        id="flu_A",
        model="Flu",
        fasta=FLU_DIR / "flu_test_seq.fa",
        metadata=FLU_DIR / "metadata_flu.txt",
        isolate="A/Japan/NIG-xxxxx/2023",
    ),
    ModelCase(
        id="flu_B",
        model="Flu",
        fasta=FLU_DIR / "fluB_PV958151-8.fasta",
        metadata=FLU_DIR / "metadata_fluB.txt",
        isolate="flu/Japan/NIG-xxxxx/2023",
    ),
]


# ---------- CLI surface tests (do not need VADR) ----------

def test_help(dfv_python: str, vadr2mss_bin: Path) -> None:
    """``vadr2mss.py --help`` prints usage."""
    result = run_cli([dfv_python, str(vadr2mss_bin), "--help"], timeout=30)
    assert result.returncode == 0, result.stderr
    assert "VADR2MSS" in result.stdout or "vadr2mss" in result.stdout
    assert "--model" in result.stdout


def test_invalid_model_rejected(dfv_python: str, vadr2mss_bin: Path, tmp_path: Path) -> None:
    """An unknown ``-M`` value is rejected by argparse before any work happens."""
    out_dir = tmp_path / "out"
    result = run_cli(
        [
            dfv_python,
            str(vadr2mss_bin),
            "-i",
            "/tmp/whatever.fa",
            "-o",
            str(out_dir),
            "-M",
            "NotAVirus",
        ],
        timeout=30,
    )
    assert result.returncode != 0
    assert "invalid choice" in result.stderr or "invalid choice" in result.stdout


def test_missing_input_aborts(dfv_python: str, vadr2mss_bin: Path, tmp_path: Path) -> None:
    """A non-existent input fasta should be reported and exit non-zero."""
    out_dir = tmp_path / "out"
    result = run_cli(
        [
            dfv_python,
            str(vadr2mss_bin),
            "-i",
            str(tmp_path / "does_not_exist.fasta"),
            "-o",
            str(out_dir),
            "-M",
            "mpox",
        ],
        timeout=30,
    )
    assert result.returncode != 0


# ---------- End-to-end pipeline tests ----------

@pytest.mark.parametrize("case", CASES, ids=[c.id for c in CASES])
def test_pipeline(dfv_python: str, vadr2mss_bin: Path, tmp_path: Path, case: ModelCase) -> None:
    """Full VADR → MSS pipeline for each virus model.

    Skipped when the corresponding VADR model directory or sample data is
    missing — this lets the same suite run in environments where only a
    subset of models has been downloaded.
    """
    require_model(case.model)
    require_path(case.fasta, f"{case.id} input fasta")
    require_path(case.metadata, f"{case.id} metadata")

    out_dir = tmp_path / f"OUT_{case.id}"
    result = run_cli(
        [
            dfv_python,
            str(vadr2mss_bin),
            "-i",
            str(case.fasta),
            "-m",
            str(case.metadata),
            "-o",
            str(out_dir),
            "-M",
            case.model,
        ],
    )
    assert result.returncode == 0, (
        f"vadr2mss.py {case.id} failed (rc={result.returncode})\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )

    # Always-produced top-level outputs.
    assert_files_exist(
        out_dir,
        [
            "annotation.gbk",
            "metadata.txt",
            "dfv_report.json",
        ],
    )
    assert (out_dir / "vadr").is_dir()

    # MSS files are named after the isolate (slashes → underscores). We can't
    # assume the input metadata still carries the documented isolate, so we
    # accept either the documented prefix or whichever single ``*.annt.tsv``
    # was produced.
    expected_prefix = isolate_to_prefix(case.isolate)
    annt_candidates = sorted(out_dir.glob("*.annt.tsv"))
    seq_candidates = sorted(out_dir.glob("*.seq.fa"))
    assert annt_candidates, "no MSS .annt.tsv file produced"
    assert seq_candidates, "no MSS .seq.fa file produced"

    annt_path = out_dir / f"{expected_prefix}.annt.tsv"
    if not annt_path.exists():
        # Fall back to whatever was produced — useful when the metadata
        # isolate string drifts (e.g. flu B example).
        annt_path = annt_candidates[0]
    seq_path = annt_path.with_suffix("").with_suffix(".seq.fa")
    if not seq_path.exists():
        seq_path = seq_candidates[0]
    assert_non_empty(annt_path)
    assert_non_empty(seq_path)

    # GenBank file sanity.
    gbk_text = (out_dir / "annotation.gbk").read_text()
    assert "LOCUS" in gbk_text

    # Report JSON shape (see vadr2mss.py L156-L158).
    report = load_json(out_dir / "dfv_report.json")
    assert "annotation" in report
    assert "warnings" in report
    annotation = report["annotation"]
    # Single-genome models report a textual status; Flu (segmented) uses
    # ``n.a.`` because per-segment evaluation differs from one-genome viruses.
    valid_statuses = {"complete", "nearly complete", "draft", "incomplete", "n.a."}
    assert annotation.get("status") in valid_statuses, annotation
    assert annotation.get("number_of_sequence", 0) >= 1
