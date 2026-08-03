import logging
import os
import sys

logger = logging.getLogger(__name__)

# DBLINK qualifiers, in output order. The TSV names them exactly as MSS does.
_DBLINK_QUALIFIERS = ("project", "biosample", "sequence read archive")

# Required by DDBJ for a COX1 barcode submission. organism and isolate also
# feed the ff_definition placeholders below, so a blank one would surface as a
# malformed DEFINITION line rather than an obvious error.
_MANDATORY_COLUMNS = ("organism", "isolate", "collection_date", "geo_loc_name")

# Qualifiers that are the same for every COX1 record. @@[...]@@ is expanded by
# the DDBJ toolchain from the qualifiers of the same source feature.
_FIXED_SOURCE_QUALIFIERS = (
    ("mol_type", "genomic DNA"),
    ("organelle", "mitochondrion"),
    (
        "ff_definition",
        "@@[organism]@@ @@[isolate]@@ mitochondrial COX1 gene for "
        "cytochrome c oxidase subunit 1, partial cds",
    ),
)

_RESERVED_COLUMNS = {key for key, _ in _FIXED_SOURCE_QUALIFIERS}


def fix_cox1_mss(work_dir, mss_file_prefix, specific_metadata, out_mss_file=None):
    """Rewrite the generic MSS annotation file with COX1-specific source features.

    ``specific_metadata`` maps a sequence id to ``{"source": {...}, "DBLINK":
    {...}}``, as built by ``cox1_helper.resolve_entries`` from the specific TSV.
    Each ``source`` feature is rebuilt from the matching row, because the generic
    conversion has no per-entry metadata and emits ``organism``/``mol_type`` with
    empty values.

    ``organism`` is mandatory in DDBJ MSS, and COX1 -- unlike the fixed-organism
    virus models (cf. ``Mpox.organism`` in ``vadr2mss_config``) -- can only get
    it from that per-entry row, because a barcode gene is sequenced from
    arbitrary species. A sequence with no matching row is therefore reported as
    an error instead of silently yielding a record without an organism.
    """

    mss_file = os.path.join(work_dir, f"{mss_file_prefix}.annt.tsv")
    if not os.path.exists(mss_file):
        logger.error(f"MSS file not found: {mss_file}")
        return

    errors = []
    matched_entries = set()
    out_buffer = []

    for feature in read_mss_file(mss_file):
        if feature[0][1] != "source":
            out_buffer.append(drop_empty_qualifiers(feature))
            continue

        seq_id, location = feature[0][0], feature[0][2]
        row = specific_metadata.get(seq_id)
        if row is None:
            errors.append(f"Sequence '{seq_id}' has no matching row in the specific TSV.")
            row = {}
        else:
            matched_entries.add(seq_id)
            # Only worth reporting when a row was found; a missing row is
            # already covered by the error above.
            missing = [
                c for c in _MANDATORY_COLUMNS if not row.get("source", {}).get(c, "").strip()
            ]
            if missing:
                errors.append(
                    f"Sequence '{seq_id}' is missing mandatory source "
                    f"column(s): {', '.join(missing)}."
                )
        metadata_dict = row.get("source", {})
        dblink = row.get("DBLINK", {})

        dblink_feature = [
            ["", "", "", qualifier, dblink[qualifier]]
            for qualifier in _DBLINK_QUALIFIERS
            if dblink.get(qualifier, "").strip()
        ]
        head = [seq_id, "source", location]
        if dblink_feature:
            dblink_feature[0][0], dblink_feature[0][1] = seq_id, "DBLINK"
            out_buffer.append(dblink_feature)
            # The entry name is already carried by the DBLINK feature above.
            head[0] = ""

        first_qualifier, *other_qualifiers = _FIXED_SOURCE_QUALIFIERS
        source_feature = [head + list(first_qualifier)]
        for key, value in metadata_dict.items():
            # Column names become qualifier names verbatim, and are deliberately
            # not checked against a list of legal qualifiers -- that is a
            # separate validator's job. The mandatory four are the exception:
            # they are looked up by name above, so a misspelt one is reported as
            # missing.
            #
            # A blank column means "not provided for this entry". MSS cannot
            # express a qualifier with an empty value, so such columns are
            # dropped rather than emitted.
            if key not in _RESERVED_COLUMNS and value.strip():
                source_feature.append(["", "", "", key, value])
        for key, value in other_qualifiers:
            source_feature.append(["", "", "", key, value])
        out_buffer.append(source_feature)

    unmatched = sorted(set(specific_metadata) - matched_entries)
    if unmatched:
        errors.append(
            "Specific TSV rows matched no sequence -- check the entry names against "
            f"the FASTA headers: {', '.join(unmatched)}"
        )

    if errors:
        for message in errors:
            logger.error(message)
        logger.error("Aborted: the MSS file would be missing mandatory qualifiers.")
        sys.exit(1)

    if out_mss_file is None:
        out_mss_file = mss_file
        logger.info(f"Overwriting original MSS file: {mss_file}")
    with open(out_mss_file, "w") as f:
        out_str = "\n".join(["\t".join(row) for feature in out_buffer for row in feature])
        f.write(out_str + "\n")
    logger.info(f"Fixed MSS file for COX1: {out_mss_file}")


def drop_empty_qualifiers(feature):
    """Drop rows whose qualifier value is blank, preserving the feature header.

    metadataUtil.renderCommonEntry emits every mss_required field even when
    unset (MetadataField.render writes three blank ab_name rows on purpose),
    which suits a template for a submitter to fill in but not a finished file:
    the DDBJ validator raises ANN2645 "Missing value for the qualifier" on any
    qualifier row without a value, unless the qualifier is declared value-less.

    No value-less qualifier can reach here. The seven in metadataDefinition.tsv
    (environmental_sample, focus, germline, macronuclear, proviral, rearranged,
    transgenic) are all source qualifiers, and the source feature is rebuilt
    from the specific TSV rather than passed through this function.
    """
    kept = [row for row in feature if row[4].strip()]
    if not kept:
        return []
    # Entry / feature / location only ever appear on a feature's first row, so
    # they have to move across if that row was the one dropped.
    kept[0] = feature[0][:3] + kept[0][3:]
    return kept


def read_mss_file(mss_file):
    """
    Read MSS file and return list of lines
    """
    def iterEntries(fileName):
        Buffer = []
        for line in open(fileName):
            rows = line.strip("\n").split("\t")
            if len(rows) < 5:
                rows += [""] * (5 - len(rows))
            if rows[0] != "" and len(Buffer) > 0:
                yield Buffer
                Buffer = []
            Buffer.append(rows)
        yield Buffer

    def iterFeatures(entry):
        Buffer = []
        for rows in entry:
            if rows[1] != "" and len(Buffer) > 0:
                yield Buffer
                Buffer = []
            Buffer.append(rows)
        yield Buffer
    for entry in iterEntries(mss_file):
        for feature in iterFeatures(entry):
            yield feature
