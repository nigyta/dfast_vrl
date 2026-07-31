#!/usr/bin/env python

"""Metadata input for cox1_to_ddbj, in the ddbj_mss_tools batch_wgs_builder format.

Two files, matching what ``batch_wgs_builder`` accepts so that the same
metadata can be produced by the mss_tools web UI:

``--common`` JSON
    SUBMITTER / REFERENCE / COMMENT / DATE shared by every entry. Trailing
    commas before ``}`` / ``]`` are tolerated (JSON5-style), as in mss_tools.

``--specific`` TSV
    Two header rows -- row 1 names the feature, row 2 the qualifier -- and one
    row per sequence.

The JSON is flattened onto dfast_vrl's ``metadata.txt`` key->value form, which
``dfv.metadataUtil.Metadata`` consumes. Those keys are dfast_vrl's own names,
not the MSS qualifier names (``submitter`` becomes ``ab_name``, ``ZIP`` becomes
``zip``, ...); see ``metadataDefinition.tsv`` for the full table.
"""

import csv
import json
import logging
import os
import re

logger = logging.getLogger(__name__)

# common JSON SUBMITTER key -> metadata.txt key. Same name unless noted.
_SUBMITTER_KEYS = {
    "consrtm": "consrtm",
    "contact": "contact",
    "email": "email",
    "url": "url",
    "institute": "institute",
    "department": "department",
    "country": "country",
    "state": "state",
    "city": "city",
    "street": "street",
    "zip": "ZIP",  # dfast_vrl spells the metadata key upper-case
}

# common JSON REFERENCE key -> metadata.txt key.
_REFERENCE_KEYS = {
    "title": "reference",
    "consrtm": "refconsrtm",
    "status": "status",
    "year": "year",
    "journal": "journal",
    "volume": "volume",
    "start_page": "start_page",
    "end_page": "end_page",
}

# TSV features we can place. Everything in the MSS COMMON block is shared by
# every entry here (cox1_to_ddbj emits one annotation file for the whole run),
# so per-row COMMENT / ST_COMMENT could not be represented and are refused
# rather than silently dropped.
_TSV_FEATURES = ("_", "source", "DBLINK")
_COMMON_ONLY_FEATURES = ("COMMENT", "ST_COMMENT", "SUBMITTER", "REFERENCE", "DATE")

_ROW_KEYS = ("_entry", "_file_path")


def load_common_json(path):
    """Read the common metadata JSON, tolerating JSON5-style trailing commas."""
    with open(path, encoding="utf-8") as f:
        raw = f.read()
    try:
        return json.loads(re.sub(r",\s*([}\]])", r"\1", raw))
    except json.JSONDecodeError as exc:
        raise ValueError(f"{path} is not valid JSON: {exc}") from exc


def common_json_to_metadata(common):
    """Flatten the common JSON into dfast_vrl's metadata.txt key->value mapping."""
    metadata = {"projectType": "cox1"}

    submitter = common.get("SUBMITTER", {})
    ab_name = submitter.get("ab_name", [])
    if isinstance(ab_name, str):  # tolerate a pre-joined string
        ab_name = [ab_name]
    if ab_name:
        metadata["submitter"] = "; ".join(ab_name)
    for json_key, metadata_key in _SUBMITTER_KEYS.items():
        value = submitter.get(json_key)
        if value:
            metadata[metadata_key] = str(value)
    if submitter.get("phone"):
        # DDBJ dropped phone/fax from the MSS submitter block.
        logger.warning("SUBMITTER.phone is not used in MSS files and is ignored.")

    references = common.get("REFERENCE", [])
    if isinstance(references, dict):  # a lone object instead of a list
        references = [references]
    for index, reference in enumerate(references, start=1):
        # Metadata.addReferences() looks up the 2nd reference onwards under
        # "<key>:<n>"; the first uses the bare key.
        suffix = "" if index == 1 else f":{index}"
        authors = reference.get("ab_name", [])
        if isinstance(authors, str):
            authors = [authors]
        if authors:
            metadata[f"author{suffix}"] = "; ".join(authors)
        for json_key, metadata_key in _REFERENCE_KEYS.items():
            value = reference.get(json_key)
            if value:
                metadata[f"{metadata_key}{suffix}"] = str(value)

    for index, comment in enumerate(_comment_blocks(common), start=1):
        suffix = "" if index == 1 else f":{index}"
        metadata[f"comment{suffix}"] = "; ".join(comment)

    hold_date = common.get("DATE", {}).get("hold_date")
    if hold_date:
        metadata["holdDate"] = str(hold_date)

    return metadata


def _comment_blocks(common):
    """Yield each COMMENT block as a list of lines, accepting mss_tools' shapes."""
    comments = common.get("COMMENT", [])
    if isinstance(comments, dict):
        comments = [comments]
    elif isinstance(comments, str):
        comments = [{"line": [comments]}]
    for comment in comments:
        lines = comment.get("line", []) if isinstance(comment, dict) else comment
        if isinstance(lines, str):
            lines = [lines]
        lines = [line for line in lines if str(line).strip()]
        if lines:
            yield [str(line) for line in lines]


def write_metadata_file(work_dir, common):
    """Write metadata.txt for MSS2 and return its path."""
    output_metadata_file = os.path.join(work_dir, "metadata.txt")
    with open(output_metadata_file, "w", encoding="utf-8") as f:
        for key, value in common_json_to_metadata(common).items():
            f.write(f"{key}\t{value}\n")
    return output_metadata_file


def parse_specific_tsv(path):
    """Parse the two-header-row TSV into one {feature: {qualifier: value}} per row.

    Returns (rows, row_key) where row_key is whichever of ``_entry`` /
    ``_file_path`` the file uses.
    """
    with open(path, encoding="utf-8", newline="") as f:
        table = [row for row in csv.reader(f, delimiter="\t") if any(cell.strip() for cell in row)]
    if len(table) < 2:
        raise ValueError(f"{path}: expected two header rows (features, then qualifiers).")

    features, qualifiers = table[0], table[1]
    if len(features) != len(qualifiers):
        raise ValueError(
            f"{path}: the two header rows have different widths "
            f"({len(features)} vs {len(qualifiers)})."
        )
    # Row 1 repeats the feature name on every column it covers, but tolerate it
    # being written only above the first column of a group.
    filled_features, current = [], ""
    for feature in features:
        current = feature.strip() or current
        filled_features.append(current)

    columns = list(zip(filled_features, [q.strip() for q in qualifiers]))
    _reject_unusable_features(path, columns)
    row_key = _resolve_row_key(path, columns)

    rows = []
    for line_no, values in enumerate(table[2:], start=3):
        row = {feature: {} for feature in _TSV_FEATURES}
        for (feature, qualifier), value in zip(columns, values):
            value = value.strip()
            if value and qualifier:
                row[feature][qualifier] = value
        if not row["_"].get(row_key):
            raise ValueError(f"{path} line {line_no}: '{row_key}' is empty.")
        rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no data rows.")
    return rows, row_key


def _reject_unusable_features(path, columns):
    for feature, qualifier in columns:
        if feature in _COMMON_ONLY_FEATURES:
            raise ValueError(
                f"{path}: '{feature}' cannot be set per row -- cox1_to_ddbj writes one "
                f"COMMON block for the whole run. Move it to the --common JSON."
            )
        if feature not in _TSV_FEATURES:
            raise ValueError(
                f"{path}: unsupported feature '{feature}' (expected one of "
                f"{', '.join(_TSV_FEATURES)})."
            )
        if feature == "_" and qualifier not in _ROW_KEYS:
            # _submission_category selects genome-registration DATATYPE/DIVISION
            # /KEYWORD, none of which apply to a COX1 gene submission.
            logger.warning("Ignoring '_/%s' column: not applicable to COX1.", qualifier)


def resolve_entries(rows, row_key, input_fasta, work_dir):
    """Map TSV rows onto sequences; return (fasta_to_annotate, specific_metadata).

    With ``_entry`` the row key is a sequence id inside the single ``-i`` FASTA.
    With ``_file_path`` each row names its own FASTA, mirroring
    batch_wgs_builder; those are concatenated into one file so that VADR still
    runs once, and every sequence in a row's file inherits that row's metadata.
    """
    if row_key == "_entry":
        if not input_fasta:
            raise ValueError("The TSV uses '_entry', so an input FASTA (-i) is required.")
        specific_metadata = {}
        for row in rows:
            entry = row["_"]["_entry"]
            if entry in specific_metadata:
                raise ValueError(f"Duplicate '_entry' in the TSV: {entry}")
            specific_metadata[entry] = {"source": row["source"], "DBLINK": row["DBLINK"]}
        return input_fasta, specific_metadata

    if input_fasta:
        raise ValueError("The TSV uses '_file_path', so an input FASTA (-i) must not be given.")

    from Bio import SeqIO

    merged_fasta = os.path.join(work_dir, "input.fasta")
    specific_metadata = {}
    records = []
    for row in rows:
        file_path = row["_"]["_file_path"]
        if not os.path.exists(file_path):
            raise ValueError(f"FASTA file not found: {file_path}")
        found = list(SeqIO.parse(file_path, "fasta"))
        if not found:
            raise ValueError(f"No sequence found in {file_path}")
        for record in found:
            if record.id in specific_metadata:
                raise ValueError(f"Duplicate sequence id across the TSV's FASTA files: {record.id}")
            specific_metadata[record.id] = {"source": row["source"], "DBLINK": row["DBLINK"]}
        records.extend(found)
    SeqIO.write(records, merged_fasta, "fasta")
    logger.info("Merged %d sequence(s) from %d file(s) into %s",
                len(records), len(rows), merged_fasta)
    return merged_fasta, specific_metadata


def _resolve_row_key(path, columns):
    present = [key for key in _ROW_KEYS if ("_", key) in columns]
    if len(present) != 1:
        raise ValueError(
            f"{path}: exactly one of the '_' columns {' / '.join(_ROW_KEYS)} is required, "
            f"found {present or 'neither'}."
        )
    return present[0]
