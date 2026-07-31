"""Operational tests for ``cox1_to_ddbj.py`` (COX1 barcode-gene pipeline).

Unlike ``vadr2mss.py``, this script fixes the VADR model to COX1 and hardcodes
the MSS file prefix to ``cox1``. Its metadata format is also different: a single
file holding a ``##COMMON`` JSON block plus a ``##SPECIFIC`` per-entry TSV block
(see ``examples/cox1/metadata_example.tsv``), parsed by ``dfv/cox1_helper.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from _helpers import assert_files_exist, assert_non_empty, load_json, run_cli
from conftest import EXAMPLES_DIR, require_model, require_path

COX1_DIR = EXAMPLES_DIR / "cox1"
COX1_FASTA = COX1_DIR / "multi.fa"            # 2 records: cox1_test, cox1_test2
COX1_METADATA = COX1_DIR / "metadata_example.tsv"

# cox1_to_ddbj.py L145: the MSS prefix is not derived from the isolate.
MSS_PREFIX = "cox1"


# ---------- MSS helpers ----------

def parse_annt(path: Path) -> list[tuple[str, str, str, str, str]]:
    """Flatten an MSS ``.annt.tsv`` into 5-tuples with entry/feature carried down.

    MSS leaves the entry, feature and location columns blank on continuation
    rows; filling them in makes the rows independently assertable.
    """
    rows: list[tuple[str, str, str, str, str]] = []
    entry = feature = location = ""
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        cols = (line.split("\t") + ["", "", "", "", ""])[:5]
        if cols[0]:
            entry = cols[0]
        if cols[1]:
            feature, location = cols[1], cols[2]
        rows.append((entry, feature, location, cols[3], cols[4]))
    return rows


def qualifiers(rows, entry: str, feature: str) -> dict[str, str]:
    return {k: v for e, f, _loc, k, v in rows if e == entry and f == feature and k}


# ---------- CLI surface tests (do not need VADR) ----------

def test_help(dfv_python: str, cox1_bin: Path) -> None:
    """``cox1_to_ddbj.py --help`` prints usage and exits cleanly."""
    result = run_cli([dfv_python, str(cox1_bin), "--help"], timeout=30)
    assert result.returncode == 0, result.stderr
    assert "cox1_to_ddbj" in result.stdout
    # The model is fixed to COX1, so -M must not be offered.
    assert "--model" not in result.stdout


def test_missing_input_aborts(dfv_python: str, cox1_bin: Path, tmp_path: Path) -> None:
    """A non-existent input fasta is reported and exits non-zero."""
    result = run_cli(
        [
            dfv_python,
            str(cox1_bin),
            "-i",
            str(tmp_path / "does_not_exist.fa"),
            "-m",
            str(COX1_METADATA),
            "-o",
            str(tmp_path / "out"),
        ],
        timeout=60,
    )
    assert result.returncode != 0
    assert "does not exist" in result.stderr


def test_existing_out_dir_requires_force(dfv_python: str, cox1_bin: Path, tmp_path: Path) -> None:
    """Refuse to write into an existing output directory unless ``--force``."""
    require_path(COX1_FASTA, "cox1 input fasta")
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    result = run_cli(
        [
            dfv_python,
            str(cox1_bin),
            "-i",
            str(COX1_FASTA),
            "-m",
            str(COX1_METADATA),
            "-o",
            str(out_dir),
        ],
        timeout=60,
    )
    assert result.returncode != 0
    # Not argparse's usage line, which also mentions --force.
    assert "already exists" in result.stderr


# ---------- End-to-end pipeline ----------

@pytest.fixture(scope="module")
def cox1_out(dfv_python: str, tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Run the full COX1 pipeline once and share the output dir module-wide."""
    from conftest import COX1_BIN

    require_model("COX1")
    require_path(COX1_FASTA, "cox1 input fasta")
    require_path(COX1_METADATA, "cox1 metadata")

    out_dir = tmp_path_factory.mktemp("cox1") / "OUT_cox1"
    result = run_cli(
        [
            dfv_python,
            str(COX1_BIN),
            "-i",
            str(COX1_FASTA),
            "-m",
            str(COX1_METADATA),
            "-o",
            str(out_dir),
        ],
    )
    assert result.returncode == 0, (
        f"cox1_to_ddbj.py failed (rc={result.returncode})\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )
    return out_dir


def test_pipeline_outputs(cox1_out: Path) -> None:
    """The documented set of output files is produced."""
    assert_files_exist(
        cox1_out,
        [
            "annotation.gbk",
            "metadata.txt",
            "dfv_report.json",
            f"{MSS_PREFIX}.annt.tsv",
            f"{MSS_PREFIX}.seq.fa",
        ],
    )
    assert (cox1_out / "vadr").is_dir()
    assert_non_empty(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    assert "LOCUS" in (cox1_out / "annotation.gbk").read_text()


def test_seq_fa_holds_both_entries(cox1_out: Path) -> None:
    """Both records of the multi-FASTA survive into the MSS sequence file."""
    headers = [
        line[1:].strip()
        for line in (cox1_out / f"{MSS_PREFIX}.seq.fa").read_text().splitlines()
        if line.startswith(">")
    ]
    assert headers == ["cox1_test", "cox1_test2"]


def test_common_block_from_metadata(cox1_out: Path) -> None:
    """The ``##COMMON`` JSON block lands in the MSS COMMON entry."""
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    submitter = qualifiers(rows, "COMMON", "SUBMITTER")
    assert submitter["contact"] == "Masanori Arita"
    assert submitter["institute"] == "National Institute of Genetics"
    assert qualifiers(rows, "COMMON", "REFERENCE")["title"] == "COX1 sequences for XXXXX"


def test_source_features_are_cox1_specific(cox1_out: Path) -> None:
    """``fix_cox1_mss`` adds the COX1/mitochondrion source qualifiers."""
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    for entry in ("cox1_test", "cox1_test2"):
        source = qualifiers(rows, entry, "source")
        assert source["mol_type"] == "genomic DNA"
        assert source["organism"] == "Genus species"
        assert source["organelle"] == "mitochondrion"

        cds = qualifiers(rows, entry, "CDS")
        assert cds["product"] == "cytochrome c oxidase subunit I"
        assert cds["gene"] == "COX1"
        assert cds["transl_table"] == "5"  # invertebrate mitochondrial code


def test_dblink_only_for_entries_with_accessions(cox1_out: Path) -> None:
    """DBLINK is emitted for cox1_test (has BioProject/BioSample/SRA) only."""
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    dblink = qualifiers(rows, "cox1_test", "DBLINK")
    assert dblink["project"] == "PRJDB9999"
    assert dblink["biosample"] == "SAMD999999"
    assert dblink["sequence read archive"] == "SRR99999999"
    # cox1_test2 leaves all three blank in the SPECIFIC block.
    assert not qualifiers(rows, "cox1_test2", "DBLINK")


def test_report_json(cox1_out: Path) -> None:
    report = load_json(cox1_out / "dfv_report.json")
    assert "annotation" in report and "warnings" in report
    assert report["annotation"]["number_of_sequence"] == 2


def test_no_empty_qualifier_values(cox1_out: Path) -> None:
    """No qualifier may be written with an empty value.

    Columns a submitter left blank for a given entry (``strain`` for cox1_test,
    ``isolate`` for cox1_test2) and unfilled ``##COMMON`` keys (REFERENCE
    author/year) must be dropped, not emitted with an empty value.
    """
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    empty = [(e, f, k) for e, f, _loc, k, v in rows if k and not v]
    assert not empty, f"qualifiers with empty values: {empty}"


# ---------- Mandatory per-entry metadata ----------
#
# organism is required by DDBJ MSS and, COX1 being a barcode gene, can only come
# from the ##SPECIFIC block. These cases used to produce a clean exit(0) and an
# MSS file whose source feature held nothing but mol_type and organelle.

def _write_metadata(path: Path, specific_rows: str) -> Path:
    """Build a metadata file with a valid COMMON block and the given SPECIFIC rows."""
    path.write_text(
        "##COMMON - DO NOT CHANGE THIS LINE\n"
        '"COMMON": {\n'
        '\t"SUBMITTER": {"submitter": "Tanizawa,Y.", "contact": "Y Tanizawa",\n'
        '\t\t"email": "mss@ddbj.nig.ac.jp", "institute": "NIG",\n'
        '\t\t"department": "DDBJ", "country": "Japan"},\n'
        '\t"REFERENCE": {"reference": "COX1 test", "author": "", "year": ""}\n'
        "}\n"
        "##SPECIFIC - DO NOT CHANGE THIS LINE\n" + specific_rows
    )
    return path


def _run_cox1(dfv_python: str, tmp_path: Path, metadata: Path) -> "object":
    from conftest import COX1_BIN

    require_model("COX1")
    require_path(COX1_FASTA, "cox1 input fasta")
    return run_cli(
        [
            dfv_python,
            str(COX1_BIN),
            "-i",
            str(COX1_FASTA),
            "-m",
            str(metadata),
            "-o",
            str(tmp_path / "out"),
        ],
    )


def test_entry_not_in_specific_block_aborts(dfv_python: str, tmp_path: Path) -> None:
    """Entry names that match no sequence (and vice versa) abort the run."""
    metadata = _write_metadata(
        tmp_path / "mismatch.tsv",
        "# entry\torganism\n"
        "KJ948760\tGenus species\n"
        "KJ948761\tGenus species\n",
    )
    result = _run_cox1(dfv_python, tmp_path, metadata)
    assert result.returncode != 0, "a full entry-name mismatch must not exit 0"
    combined = result.stdout + result.stderr
    # Both directions are reported: unmatched sequences and unmatched rows.
    assert "cox1_test" in combined
    assert "KJ948760" in combined


def test_missing_organism_aborts(dfv_python: str, tmp_path: Path) -> None:
    """An entry that is listed but carries no organism aborts the run."""
    metadata = _write_metadata(
        tmp_path / "no_organism.tsv",
        "# entry\torganism\tnote\n"
        "cox1_test\tGenus species\tfine\n"
        "cox1_test2\t\tno organism here\n",
    )
    result = _run_cox1(dfv_python, tmp_path, metadata)
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "organism" in combined
    assert "cox1_test2" in combined


def test_metadata_file_is_required(dfv_python: str, cox1_bin: Path, tmp_path: Path) -> None:
    """``-m`` is mandatory: organism has no other source.

    It used to be advertised as optional and then crash with a TypeError, after
    VADR had already run.
    """
    result = run_cli(
        [dfv_python, str(cox1_bin), "-i", str(COX1_FASTA), "-o", str(tmp_path / "out")],
        timeout=60,
    )
    assert result.returncode != 0
    assert "metadata_file" in result.stderr
    assert "Traceback" not in result.stderr
