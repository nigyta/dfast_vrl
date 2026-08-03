#!/usr/bin/env python

import os
import sys
import argparse
import logging
import json
from dfv.cox1_helper import (
    load_common_json,
    parse_specific_tsv,
    resolve_entries,
    write_metadata_file,
)
from dfv.fix_mss_for_cox1 import fix_cox1_mss

VERSION = "1.7.0-0.1"


parser = argparse.ArgumentParser(prog="cox1_to_ddbj.py",
                                 description=f'COX1 to DDBJ converter: Pipeline for COX1 gene submission to DDBJ. (ver. {VERSION})'
                                 )
parser.add_argument(
    "-i", "--input", metavar="FILE",
    help="Input single/multi-FASTA file. Required when the specific TSV keys rows "
         "by '_entry'; omit it when the TSV keys rows by '_file_path'.")
parser.add_argument("-o", "--out_dir", metavar="FILE", help="Output Directory", required=True)
parser.add_argument(
    '--isolate',
    type=str,
    help='Isolate name [Optional, but recommendable to specify a unique value to distinguish the sample]',
    metavar="STR")
parser.add_argument(
        "-m",
        "--common",
        type=str,
        help="Common metadata JSON (SUBMITTER, REFERENCE, COMMENT, DATE), "
             "in the ddbj_mss_tools batch_wgs_builder format",
        metavar="PATH",
        required=True
    )
parser.add_argument(
        "-s",
        "--specific",
        type=str,
        help="Per-sequence metadata TSV with two header rows "
             "(row 1 = feature, row 2 = qualifier)",
        metavar="PATH",
        required=True
    )
# model is fixed to COX1, so no need to specify model argument
# parser.add_argument(
#         "-M",
#         "--model",
#         type=str,
#         help="Reference model for VADR [mpox, sarscov2, corona, RSV, Noro, Calici, Dengue, Flavi, COX1, Flu]",
#         metavar="MODEL",
#         required=True,
#         choices=["mpox", "sarscov2", "corona", "RSV", "Noro", "Calici", "Dengue", "Flavi", "COX1", "Flu"]
#     )
parser.add_argument(
        '--force',
        action='store_true',
        help='Force overwriting result'
    )
parser.add_argument(
        '--debug',
        action='store_true',
        help='Debug mode'
    )


if len(sys.argv)==1:
    parser.print_help()
    exit()

args = parser.parse_args()

from Bio import SeqIO

from dfv.vadr_generic import main as run_vadr
from dfv.tbl2genbank import tbl2genbank
from dfv.genbank2mss import MSS2
from dfv.common import get_isolate  # , get_logger
from dfv.vadr2mss_config import models
from dfv.check_cox1_annotation import check_cox1_annotation
from dfv.cox1_transl_table import parse_ftr_models, transl_table_by_entry

# cheking input files
input_fasta = args.input
if input_fasta and not os.path.exists(input_fasta):
    sys.stderr.write(f"cox1_to_ddbj.py: error: Specified FASTA file does not exist. [{input_fasta}]\n")
    exit(1)
for label, path in (("common metadata", args.common), ("specific metadata", args.specific)):
    if not os.path.exists(path):
        sys.stderr.write(f"cox1_to_ddbj.py: error: Specified {label} file does not exist. [{path}]\n")
        exit(1)


# setting output directory
work_dir = args.out_dir
work_dir = work_dir.rstrip("/")
if os.path.exists(work_dir) and not args.force:
    if not args.force:
        sys.stderr.write("cox1_to_ddbj.py: error: Output directory already exists. Use '--force' to overwrite.\n")
        exit(1)
else:
    os.makedirs(work_dir, exist_ok=True)


def get_logger(name=None, debug=False):
    if debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logger = logging.getLogger(__name__)
    sh = logging.StreamHandler(stream=sys.stdout)
    log_file = os.path.join(work_dir, "application.log")
    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8", delay=True)

    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        level=log_level,
        handlers=[sh, fh]) 
    logger = logging.getLogger(__name__)
    return logger

logger = get_logger(name=__name__, debug=args.debug)
logger.info("cox1_to_ddbj started. version=%s", VERSION)


# setting model to used
model = models["COX1"]  # fixed to COX1 model
logger.info("Selected VADR Model: %s", "COX1")




# Read the metadata before running VADR: a malformed TSV should fail in
# seconds rather than after a full annotation run.
try:
    common_metadata = load_common_json(args.common)
    specific_rows, row_key = parse_specific_tsv(args.specific)
    input_fasta, specific_metadata = resolve_entries(
        specific_rows, row_key, input_fasta, work_dir)
except ValueError as e:
    logger.error("%s", e)
    exit(1)

output_metadata_file = write_metadata_file(work_dir, common_metadata)
logger.debug("common_metadata: %s", common_metadata)
logger.debug("specific_metadata: %s", specific_metadata)


# Run VADR
vadr_dir = os.path.join(work_dir, "vadr")
vadr_out_tbl_pass, vadr_out_fasta_pass, vadr_warnings = run_vadr(input_fasta, vadr_dir, model, cpu=1)

# Convert VADR result into .gbk. The COX1 model set spans six genetic codes, so
# the translation table comes from whichever model VADR matched each sequence to
# rather than from one constant.
models_by_entry = parse_ftr_models(vadr_dir)
table_by_entry = transl_table_by_entry(vadr_dir, model)
out_gbk_file = os.path.join(work_dir, "annotation.gbk")
tbl2genbank(vadr_out_tbl_pass, vadr_out_fasta_pass, out_gbk_file, model, table_by_entry)

# metadata fileからisolate名を取得、出力されるMSSファイルのprefixはisolate名に基づいて決定
# isolate, mss_file_prefix = get_isolate(metadata_file, args)

# check_annotation_stats is not used here: it reports one set of figures for the
# whole run against a single reference model, which does not describe a
# submission of independent barcode sequences (its query_coverage summed every
# query length over one model length). check_cox1_annotation reports per entry.
annotated_records = list(SeqIO.parse(out_gbk_file, "genbank"))
entry_reports = check_cox1_annotation(annotated_records, models_by_entry)
logger.debug("Entry reports: %s", entry_reports)

# No model-level organism for COX1: it is a barcode gene sequenced from
# arbitrary species, so the organism differs per entry and comes from the
# specific TSV. fix_cox1_mss applies it and rejects entries that lack one.
# update_metadata_file() is likewise unusable here -- it writes a single
# organism for the whole run and a "complete genome" ff_definition.


# Convert .gbk file and metadata file into MSS format.
mss = MSS2(out_gbk_file, output_metadata_file)
mss_file_prefix = "cox1"
mss.convert(work_dir, mss_file_prefix)

fix_cox1_mss(work_dir, mss_file_prefix, specific_metadata)

# fix MSS file for influenza virus (adding segment information, etc.)
# if args.model == "Flu":
#     from dfv.fix_mss_for_ful import fix_flu_mss
#     fix_flu_mss(work_dir, mss_file_prefix)

warnings = [warning.to_tuple() for warning in vadr_warnings]
if vadr_warnings:
    for warning in vadr_warnings:
        logger.warning("VADR Warning: %s", warning)

out_report_file = os.path.join(work_dir, "dfv_report.json")
logger.info(f"Writing report json to {out_report_file}")
with open(out_report_file, "w") as f:
    json.dump(
        {
            "annotation": {"number_of_sequence": len(entry_reports)},
            "entries": entry_reports,
            "warnings": warnings,
        },
        f, indent=4,
    )

# Convert MSS file into JSON
# ver3.10以上ならモジュールインポートしてjson出力
# python_version = sys.version_info
# if python_version.major >= 3 and python_version.minor >= 10:
#     try:
#         from dr_tools import drt_ann2json
#         ann_file = os.path.join(work_dir, f"{mss_file_prefix}.annt.tsv")
#         seq_file = os.path.join(work_dir, f"{mss_file_prefix}.seq.fa")
#         out_json_file = os.path.join(work_dir, "dfast_record.json")
#         drt_ann2json(ann_file, seq_file, out_json_file, division="VRL")
#         logger.info(f"Converted MSS file into JSON. {ann_file} --> {out_json_file}")

#     except ImportError:
#         logger.warning("Failed to import dr_tools. Skip converting MSS to JSON.")


# else:
#     logger.warning("Python version is less than 3.10. Skip converting MSS to JSON.")
from dfv.common import mss2json
mss2json(work_dir, mss_file_prefix)

logger.info("cox1_to_ddbj finished.")