"""Per-entry content checks for a COX1 barcode submission.

Replaces check_annotation_stats for this pipeline. That module reports one set
of figures for the whole run against a single reference model, which is
meaningless here: a COX1 submission is many independent barcode sequences, each
matched to its own model, so `query_coverage` came out as the sum of all query
lengths over one model length (106.55% for two sequences) and `cds_completeness`
counted intact/partial across every sequence against one model's CDS count.

What is worth reporting instead is per entry, and centres on the translation --
the one thing that silently goes wrong when the genetic code is picked wrongly.
"""

import logging

from Bio.Data import CodonTable
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

_AMBIGUOUS_WARN_FRACTION = 0.01  # 1% is the usual barcode-database ceiling


def check_cox1_annotation(records, models_by_entry=None):
    """Return a per-entry report for the annotated records.

    ``records`` are the SeqRecords written to annotation.gbk. Nothing here fails
    the run: the findings are advisory and go into dfv_report.json, while the
    hard requirements (mandatory qualifiers, entry/metadata agreement) are
    enforced in fix_cox1_mss.
    """
    models_by_entry = models_by_entry or {}
    entries = []
    for record in records:
        entry = {
            "entry": record.id,
            "model": models_by_entry.get(record.id),
            "length": len(record.seq),
            "ambiguous_bases": _count_ambiguous(record.seq),
        }
        entry.update(_check_cds(record))
        entries.append(entry)
        _log_findings(entry)
    return entries


def _count_ambiguous(seq):
    return sum(1 for base in str(seq).upper() if base not in "ACGT")


def _check_cds(record):
    """Translate the CDS under its own transl_table and count internal stops."""
    for feature in record.features:
        if feature.type != "CDS":
            continue
        transl_table = _first(feature.qualifiers.get("transl_table"), default=1)
        try:
            protein = feature.translate(record.seq, cds=False)
        except (CodonTable.TranslationError, ValueError) as exc:
            logger.warning("Could not translate the CDS of '%s': %s", record.id, exc)
            return {"transl_table": transl_table, "internal_stop_codons": None}
        # A stop in the final position is the terminator, not a defect. A partial
        # CDS (the usual case for a barcode) simply has none.
        return {
            "transl_table": transl_table,
            "internal_stop_codons": str(protein).rstrip("*").count("*"),
        }
    return {"transl_table": None, "internal_stop_codons": None}


def _first(values, default=None):
    if not values:
        return default
    value = values[0]
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def _log_findings(entry):
    name = entry["entry"]
    stops = entry.get("internal_stop_codons")
    if stops:
        # Points at a frameshift, the wrong codon_start, or a model whose genetic
        # code does not suit the sequence.
        logger.warning(
            "%s: %d internal stop codon(s) translating with transl_table %s (model %s).",
            name, stops, entry.get("transl_table"), entry.get("model"),
        )
    ambiguous, length = entry["ambiguous_bases"], entry["length"]
    if length and ambiguous / length > _AMBIGUOUS_WARN_FRACTION:
        logger.warning(
            "%s: %d/%d bases (%.1f%%) are ambiguous.",
            name, ambiguous, length, 100 * ambiguous / length,
        )
