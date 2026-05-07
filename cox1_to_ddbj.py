#!/usr/bin/env python

import os
import sys
import argparse
import logging
import json
from dfv.cox1_helper import create_metadata_file
from dfv.fix_mss_for_cox1 import fix_cox1_mss

VERSION = "1.6.4-0.11"


parser = argparse.ArgumentParser(prog="cox1_to_ddbj.py",
                                 description=f'COX1 to DDBJ converter: Pipeline for COX1 gene submission to DDBJ. (ver. {VERSION})'
                                 )
parser.add_argument("-i", "--input", metavar="FILE", help="Input single/multi-FASTA file", required=True)
parser.add_argument("-o", "--out_dir", metavar="FILE", help="Output Directory", required=True)
parser.add_argument(
    '--isolate',
    type=str,
    help='Isolate name [Optional, but recommendable to specify a unique value to distinguish the sample]',
    metavar="STR")
parser.add_argument(
        "-m",
        "--metadata_file",
        type=str,
        help="Metadata file (Tab-separated table) [Optional]",
        metavar="PATH"
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

from dfv.vadr_generic import main as run_vadr
from dfv.tbl2genbank import tbl2genbank
from dfv.genbank2mss import MSS2
from dfv.common import get_isolate  # , get_logger
from dfv.vadr2mss_config import models
from dfv.check_annotation_stats import check_annotation_stats

# cheking input fasta
input_fasta = args.input
if not os.path.exists(input_fasta):
    sys.stderr.write(f"cox1_to_ddbj.py: error: Specified FASTA file does not exist. [{input_fasta}]\n")
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




# Run VADR
vadr_dir = os.path.join(work_dir, "vadr")
vadr_out_tbl_pass, vadr_out_fasta_pass, vadr_warnings = run_vadr(input_fasta, vadr_dir, model, cpu=1)

# Convert VADR result into .gbk
out_gbk_file = os.path.join(work_dir, "annotation.gbk")
tbl2genbank(vadr_out_tbl_pass, vadr_out_fasta_pass, out_gbk_file, model)



# prepare metadata file
# metadata_file = copy_or_create_metadata_file(work_dir, args)
output_metadata_file, common_metadata, specific_metadata = create_metadata_file(work_dir, args)

print("common_metadata:", common_metadata)
print("specific_metadata:", specific_metadata)

# metadata fileからisolate名を取得、出力されるMSSファイルのprefixはisolate名に基づいて決定
# isolate, mss_file_prefix = get_isolate(metadata_file, args)

annotation_stats = check_annotation_stats(work_dir, vadr_dir, model)
# {'status': 'complete', 'total_length': 10735, 'model_length': 10735, 'query_coverage': '100.00%', 'qap_length': 0, 'cds_completeness': '1 / 0 / 1 [intact/partial/expected]'}
logger.debug("Annotation stats: %s", annotation_stats)
seq_status, number_of_sequence = annotation_stats["status"], annotation_stats["number_of_sequence"]

organism = model.organism if hasattr(model, "organism") else None
    
# update_metadata_file(metadata_file, seq_status=seq_status, number_of_sequence=number_of_sequence, mol_type=model.mol_type, organism=organism)


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
    json.dump({"annotation": annotation_stats, "warnings": warnings}, f, indent=4)

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