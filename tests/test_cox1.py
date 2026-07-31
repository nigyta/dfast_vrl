"""Operational tests for ``cox1_to_ddbj.py`` (COX1 barcode-gene pipeline).

Unlike ``vadr2mss.py``, this script fixes the VADR model to COX1 and hardcodes
the MSS file prefix to ``cox1``. Its metadata format is also different: a common
metadata JSON plus a per-sequence TSV with two header rows, matching
ddbj_mss_tools' batch_wgs_builder (see ``examples/cox1/``), parsed by
``dfv/cox1_helper.py``.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from _helpers import assert_files_exist, assert_non_empty, load_json, run_cli
from conftest import EXAMPLES_DIR, require_model, require_path

COX1_DIR = EXAMPLES_DIR / "cox1"
COX1_FASTA = COX1_DIR / "multi.fa"            # 2 records: cox1_test, cox1_test2
COX1_COMMON = COX1_DIR / "common_example.json"
COX1_SPECIFIC = COX1_DIR / "specific_example_with_dblink.tsv"
COX1_SPECIFIC_NO_DBLINK = COX1_DIR / "specific_example_without_dblink.tsv"


def read_common_example() -> dict:
    import json

    return json.loads(COX1_COMMON.read_text())


def read_specific_example(path: Path = COX1_SPECIFIC) -> dict[str, dict[str, dict[str, str]]]:
    """Parse an example specific TSV into {entry: {feature: {qualifier: value}}}.

    Expectations are derived from the example files rather than hardcoded, so
    editing the examples does not break these tests.
    """
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    features, qualifiers = lines[0].split("\t"), lines[1].split("\t")
    out: dict[str, dict[str, dict[str, str]]] = {}
    for line in lines[2:]:
        values = line.split("\t")
        row: dict[str, dict[str, str]] = {}
        entry = ""
        for feature, qualifier, value in zip(features, qualifiers, values):
            if not value.strip():
                continue
            if qualifier == "_entry":
                entry = value.strip()
            else:
                row.setdefault(feature, {})[qualifier] = value.strip()
        out[entry] = row
    return out

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
            str(COX1_COMMON),
            "-s",
            str(COX1_SPECIFIC),
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
            str(COX1_COMMON),
            "-s",
            str(COX1_SPECIFIC),
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
    require_path(COX1_COMMON, "cox1 common metadata")
    require_path(COX1_SPECIFIC, "cox1 specific metadata")

    out_dir = tmp_path_factory.mktemp("cox1") / "OUT_cox1"
    result = run_cli(
        [
            dfv_python,
            str(COX1_BIN),
            "-i",
            str(COX1_FASTA),
            "-m",
            str(COX1_COMMON),
            "-s",
            str(COX1_SPECIFIC),
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
    """The common metadata JSON lands in the MSS COMMON entry."""
    common = read_common_example()
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")

    submitter = qualifiers(rows, "COMMON", "SUBMITTER")
    assert submitter["contact"] == common["SUBMITTER"]["contact"]
    assert submitter["institute"] == common["SUBMITTER"]["institute"]
    # JSON "zip" -> metadata key "ZIP" -> MSS qualifier "zip".
    assert submitter["zip"] == common["SUBMITTER"]["zip"]
    submitter_names = [v for _e, f, _l, k, v in rows if f == "SUBMITTER" and k == "ab_name"]
    assert submitter_names == common["SUBMITTER"]["ab_name"]

    reference = common["REFERENCE"][0]
    assert qualifiers(rows, "COMMON", "REFERENCE")["title"] == reference["title"]
    # ab_name is a JSON list; every author must survive the join/split round-trip.
    authors = [v for _e, f, _l, k, v in rows if f == "REFERENCE" and k == "ab_name"]
    assert authors == reference["ab_name"]

    expected_comment = common["COMMENT"][0]["line"]
    assert qualifiers(rows, "COMMON", "COMMENT")["line"] == "; ".join(expected_comment)


def test_source_features_are_cox1_specific(cox1_out: Path) -> None:
    """``fix_cox1_mss`` adds the COX1/mitochondrion source qualifiers."""
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    expected = read_specific_example()
    for entry in expected:
        source = qualifiers(rows, entry, "source")
        assert source["mol_type"] == "genomic DNA"
        assert source["organism"] == expected[entry]["source"]["organism"]
        assert source["organelle"] == "mitochondrion"

        # DDBJ builds the DEFINITION line from this; @@[...]@@ is expanded from
        # the qualifiers of the same source feature.
        assert source["ff_definition"] == (
            "@@[organism]@@ @@[isolate]@@ mitochondrial COX1 gene for "
            "cytochrome c oxidase subunit 1, partial cds"
        )
        # Mandatory for a COX1 barcode submission.
        for mandatory in ("organism", "isolate", "collection_date", "geo_loc_name"):
            assert source.get(mandatory), f"{entry} lacks {mandatory}"

        cds = qualifiers(rows, entry, "CDS")
        assert cds["product"] == "cytochrome c oxidase subunit I"
        assert cds["gene"] == "COX1"
        assert cds["transl_table"] == "5"  # invertebrate mitochondrial code


def test_dblink_only_for_entries_with_accessions(cox1_out: Path) -> None:
    """DBLINK appears exactly for the entries whose TSV row carries accessions."""
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    expected = read_specific_example()
    assert any("DBLINK" in row for row in expected.values()), (
        f"{COX1_SPECIFIC.name} no longer exercises DBLINK"
    )
    for entry, row in expected.items():
        assert qualifiers(rows, entry, "DBLINK") == row.get("DBLINK", {})


def test_report_json(cox1_out: Path) -> None:
    report = load_json(cox1_out / "dfv_report.json")
    assert "annotation" in report and "warnings" in report
    assert report["annotation"]["number_of_sequence"] == 2


def test_no_empty_qualifier_values_in_per_entry_features(cox1_out: Path) -> None:
    """A blank TSV column must be dropped, not emitted as an empty qualifier.

    Scoped to the per-entry features. COMMON is exempt on purpose: metadataUtil
    renders every mss_required field even when unset, as a template for the
    submitter to fill in, and dr_tools rejects the record if those rows are
    missing.
    """
    rows = parse_annt(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    empty = [
        (e, f, k) for e, f, _loc, k, v in rows if k and not v and e != "COMMON"
    ]
    assert not empty, f"qualifiers with empty values: {empty}"


def test_dfast_record_json_written_when_organism_is_shared(
    dfv_python: str, tmp_path: Path
) -> None:
    """mss2json must succeed when every entry shares one organism.

    Also the regression guard for the REFERENCE ``year`` placeholder: dropping it
    from COMMON made dr_tools reject the record with
    ``COMMON.REFERENCE.0.year: Field required`` and no JSON was written.
    """
    specific = write_specific_tsv(
        tmp_path / "shared_organism.tsv",
        [
            {"_entry": "cox1_test", **COMPLETE_ROW},
            {"_entry": "cox1_test2", **COMPLETE_ROW, "isolate": "isolate-2"},
        ],
    )
    result = _run_cox1(dfv_python, tmp_path, specific)
    assert result.returncode == 0, result.stdout + result.stderr
    log = (tmp_path / "out" / "application.log").read_text()
    assert "Failed to convert MSS to JSON" not in log, log[-2000:]
    assert (tmp_path / "out" / "dfast_record.json").exists()


def test_differing_organisms_still_produce_valid_mss(cox1_out: Path) -> None:
    """Per-species organisms are the normal COX1 case and must not fail the run.

    dr_tools hoists the source qualifiers shared by every entry into
    COMMON_SOURCE and its DdbjRecord model requires ``organism`` there, so a
    submission of several species -- the whole point of barcoding -- cannot
    produce dfast_record.json. The MSS files are unaffected (ddbj-validator
    passes them), so the pipeline warns and carries on rather than failing.
    """
    organisms = {
        row["source"]["organism"] for row in read_specific_example().values()
    }
    log = (cox1_out / "application.log").read_text()
    if len(organisms) > 1:
        assert "Failed to convert MSS to JSON" in log
        assert not (cox1_out / "dfast_record.json").exists()
    # Either way the MSS deliverables exist and the run succeeded.
    assert_non_empty(cox1_out / f"{MSS_PREFIX}.annt.tsv")
    assert_non_empty(cox1_out / f"{MSS_PREFIX}.seq.fa")


def test_without_dblink_example(dfv_python: str, tmp_path: Path) -> None:
    """The DBLINK-less example runs and emits no DBLINK feature at all."""
    require_path(COX1_SPECIFIC_NO_DBLINK, "cox1 specific metadata (no DBLINK)")
    result = _run_cox1(dfv_python, tmp_path, COX1_SPECIFIC_NO_DBLINK)
    assert result.returncode == 0, result.stdout + result.stderr
    annt = tmp_path / "out" / f"{MSS_PREFIX}.annt.tsv"
    rows = parse_annt(annt)
    assert not [r for r in rows if r[1] == "DBLINK"]
    for entry in read_specific_example(COX1_SPECIFIC_NO_DBLINK):
        assert qualifiers(rows, entry, "source")["organelle"] == "mitochondrion"


# ---------- Mandatory per-entry metadata ----------
#
# organism / isolate / collection_date / geo_loc_name are required by DDBJ for a
# COX1 barcode submission and, COX1 being a barcode gene, can only come from the
# specific TSV. These cases used to produce a clean exit(0) and an MSS file
# whose source feature held nothing but mol_type and organelle.

MANDATORY_COLUMNS = ("organism", "isolate", "collection_date", "geo_loc_name")

COMPLETE_ROW = {
    "organism": "Genus species",
    "isolate": "isolate-1",
    "collection_date": "2023-01-01",
    "geo_loc_name": "Japan:Tokyo",
}


def write_specific_tsv(path: Path, rows: list[dict], row_key: str = "_entry") -> Path:
    """Render a two-header-row specific TSV, blanking columns a row omits."""
    columns = [("_", row_key)] + [("source", c) for c in MANDATORY_COLUMNS]
    lines = [
        "\t".join(feature for feature, _ in columns),
        "\t".join(qualifier for _, qualifier in columns),
    ]
    lines += ["\t".join(row.get(q, "") for _, q in columns) for row in rows]
    return _written(path, "\n".join(lines) + "\n")


def _written(path: Path, text: str) -> Path:
    path.write_text(text)
    return path


def _run_cox1(dfv_python: str, tmp_path: Path, specific: Path, fasta: Path | None = COX1_FASTA):
    from conftest import COX1_BIN

    require_model("COX1")
    require_path(COX1_FASTA, "cox1 input fasta")
    cmd = [dfv_python, str(COX1_BIN), "-m", str(COX1_COMMON), "-s", str(specific),
           "-o", str(tmp_path / "out")]
    if fasta is not None:
        cmd[2:2] = ["-i", str(fasta)]
    return run_cli(cmd)


def test_entry_not_in_specific_block_aborts(dfv_python: str, tmp_path: Path) -> None:
    """Entry names that match no sequence (and vice versa) abort the run."""
    specific = write_specific_tsv(
        tmp_path / "mismatch.tsv",
        [
            {"_entry": "KJ948760", **COMPLETE_ROW},
            {"_entry": "KJ948761", **COMPLETE_ROW},
        ],
    )
    result = _run_cox1(dfv_python, tmp_path, specific)
    assert result.returncode != 0, "a full entry-name mismatch must not exit 0"
    combined = result.stdout + result.stderr
    # Both directions are reported: unmatched sequences and unmatched rows.
    assert "cox1_test" in combined
    assert "KJ948760" in combined


@pytest.mark.parametrize("column", MANDATORY_COLUMNS)
def test_missing_mandatory_column_aborts(dfv_python: str, tmp_path: Path, column: str) -> None:
    """Blanking any one mandatory column aborts and names both column and entry."""
    specific = write_specific_tsv(
        tmp_path / f"missing_{column}.tsv",
        [
            {"_entry": "cox1_test", **COMPLETE_ROW},
            {"_entry": "cox1_test2", **COMPLETE_ROW, column: ""},
        ],
    )
    result = _run_cox1(dfv_python, tmp_path, specific)
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert column in combined
    assert "cox1_test2" in combined


def test_metadata_files_are_required(dfv_python: str, cox1_bin: Path, tmp_path: Path) -> None:
    """``-m`` and ``-s`` are mandatory: organism has no other source.

    The metadata file used to be advertised as optional and then crash with a
    TypeError, after VADR had already run.
    """
    result = run_cli(
        [dfv_python, str(cox1_bin), "-i", str(COX1_FASTA), "-o", str(tmp_path / "out")],
        timeout=60,
    )
    assert result.returncode != 0
    assert "--common" in result.stderr and "--specific" in result.stderr
    assert "Traceback" not in result.stderr


# ---------- batch_wgs_builder-compatible TSV handling ----------

def test_file_path_row_key(dfv_python: str, tmp_path: Path) -> None:
    """``_file_path`` keys rows by FASTA, as batch_wgs_builder does, with no -i."""
    require_path(COX1_FASTA, "cox1 input fasta")
    from Bio import SeqIO

    paths = []
    for record in SeqIO.parse(str(COX1_FASTA), "fasta"):
        single = tmp_path / f"{record.id}.fa"
        SeqIO.write([record], str(single), "fasta")
        paths.append(single)

    specific = write_specific_tsv(
        tmp_path / "by_file.tsv",
        [{"_file_path": str(p), **COMPLETE_ROW, "isolate": f"iso-{i}"}
         for i, p in enumerate(paths, 1)],
        row_key="_file_path",
    )
    result = _run_cox1(dfv_python, tmp_path, specific, fasta=None)
    assert result.returncode == 0, result.stdout + result.stderr

    rows = parse_annt(tmp_path / "out" / f"{MSS_PREFIX}.annt.tsv")
    assert qualifiers(rows, "cox1_test", "source")["isolate"] == "iso-1"
    assert qualifiers(rows, "cox1_test2", "source")["isolate"] == "iso-2"


def test_file_path_rejects_input_fasta(dfv_python: str, tmp_path: Path) -> None:
    """``-i`` together with ``_file_path`` is ambiguous and is refused."""
    specific = write_specific_tsv(
        tmp_path / "by_file.tsv",
        [{"_file_path": str(COX1_FASTA), **COMPLETE_ROW}],
        row_key="_file_path",
    )
    result = _run_cox1(dfv_python, tmp_path, specific)
    assert result.returncode != 0
    assert "_file_path" in result.stdout + result.stderr


def test_common_only_feature_in_tsv_is_refused(dfv_python: str, tmp_path: Path) -> None:
    """COMMENT is shared by every entry, so a per-row column cannot be honoured."""
    specific = _written(
        tmp_path / "with_comment.tsv",
        "_\tsource\tCOMMENT\n_entry\torganism\tline\ncox1_test\tGenus species\thi\n",
    )
    result = _run_cox1(dfv_python, tmp_path, specific)
    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "COMMENT" in combined and "--common" in combined


def test_malformed_tsv_fails_before_vadr(dfv_python: str, tmp_path: Path) -> None:
    """A TSV with neither row key aborts in seconds, not after annotating."""
    specific = _written(
        tmp_path / "no_key.tsv", "_\tsource\n_nope\torganism\nx\tGenus species\n"
    )
    result = _run_cox1(dfv_python, tmp_path, specific)
    assert result.returncode != 0
    assert "_entry" in result.stdout + result.stderr
    # VADR never ran.
    assert not (tmp_path / "out" / "vadr").exists()
