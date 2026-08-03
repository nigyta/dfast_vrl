"""Per-entry translation table for COX1, taken from the model VADR selected.

The COX1 model set spans six genetic codes -- 52 of its 86 models use table 5,
12 use table 2, 11 use 4, 9 use 9, and one each 13 and 14 -- and v-annotate.pl
picks a model per sequence. Annotating every record with one hardcoded table
therefore mistranslates whichever sequences land on a different one, and does so
silently: tables that differ only in ATA (Ile/Met) and AGA/AGG (Arg/Ser) yield no
stop codon to trip over, so neither VADR nor the DDBJ validator objects.

``vadr.vadr.ftr`` gives sequence -> model and the model info file gives model ->
transl_table, so no extra input is needed beyond what the run already produced.
"""

import logging
import os
import re

logger = logging.getLogger(__name__)

_FTR_SEQ_NAME = 1
_FTR_PASS_FAIL = 3
_FTR_MODEL = 4


def parse_ftr_models(vadr_dir):
    """Map sequence name -> VADR model name, from ``vadr.vadr.ftr``."""
    ftr_file = os.path.join(vadr_dir, "vadr.vadr.ftr")
    if not os.path.exists(ftr_file):
        logger.warning("FTR file not found: %s", ftr_file)
        return {}

    models = {}
    with open(ftr_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.split()
            if len(cols) <= _FTR_MODEL:
                continue
            seq_name, model = cols[_FTR_SEQ_NAME], cols[_FTR_MODEL]
            previous = models.setdefault(seq_name, model)
            if previous != model:
                # One sequence is annotated against one model; disagreement means
                # the file is not what we think it is, so do not guess.
                logger.warning(
                    "Sequence '%s' maps to more than one model (%s, %s).",
                    seq_name, previous, model,
                )
    return models


def parse_minfo_transl_tables(minfo_file):
    """Map VADR model name -> transl_table, from the model info file."""
    tables = {}
    with open(minfo_file) as f:
        for line in f:
            if not line.startswith("MODEL "):
                continue
            cols = line.split(None, 2)
            match = re.search(r'transl_table:"(\d+)"', cols[2]) if len(cols) > 2 else None
            if match:
                tables[cols[1]] = int(match.group(1))
    return tables


def transl_table_by_entry(vadr_dir, model):
    """Return {sequence name: transl_table} for the sequences VADR annotated.

    Falls back to an empty mapping if either input is unreadable; callers then
    keep the model-level default.
    """
    minfo_file = os.path.expandvars(model.minfo_file)
    if not os.path.exists(minfo_file):
        logger.warning("Model info file not found, keeping the default transl_table: %s",
                       minfo_file)
        return {}

    tables = parse_minfo_transl_tables(minfo_file)
    by_entry = {}
    for seq_name, model_name in parse_ftr_models(vadr_dir).items():
        transl_table = tables.get(model_name)
        if transl_table is None:
            logger.warning("No transl_table for model '%s'; keeping the default.", model_name)
            continue
        by_entry[seq_name] = transl_table

    default = getattr(model, "transl_table", None)
    differing = {e: t for e, t in by_entry.items() if t != default}
    if differing:
        logger.info(
            "transl_table taken from the VADR model instead of the default (%s): %s",
            default, ", ".join(f"{e}={t}" for e, t in sorted(differing.items())),
        )
    return by_entry
