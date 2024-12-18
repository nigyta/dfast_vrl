#!/usr/bin/env python

import os
import sys
import argparse
import logging
import json
VERSION = "1.6.4-0.9"


parser = argparse.ArgumentParser(prog="vadr2mss.py",
                                 description=f'VADR2MSS: Pipeline for various kind of viruses. (ver. {VERSION})'
                                 )
parser.add_argument("-i", "--input", metavar="FILE", help="Input FASTA file", required=True)
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
parser.add_argument(
        "-M",
        "--model",
        type=str,
        help="Reference model for VADR",
        metavar="MODEL",
        required=True,
        choices=["mpox", "sarscov2", "corona", "RSV", "Noro", "Calici", "Dengue", "Flavi", "COX1"]
    )
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
from dfv.common import update_metadata_file, get_isolate, copy_or_create_metadata_file  # , get_logger
from dfv.vadr2mss_config import models
from dfv.check_annotation_stats import check_annotation_stats

# cheking input fasta
input_fasta = args.input
if not os.path.exists(input_fasta):
    sys.stderr.write(f"vadr2mss.py: error: Specified FASTA file does not exist. [{input_fasta}]\n")
    exit(1)


# setting output directory
work_dir = args.out_dir
work_dir = work_dir.rstrip("/")
if os.path.exists(work_dir) and not args.force:
    if not args.force:
        sys.stderr.write("vadr2mss.py: error: Output directory already exists. Use '--force' to overwrite.\n")
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
    log_file = os.path.join(work_dir, "dfast_vrl.log")
    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8", delay=True)

    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        level=log_level,
        handlers=[sh, fh]) 
    logger = logging.getLogger(__name__)
    return logger

logger = get_logger(name=__name__, debug=args.debug)
logger.info("vadr2mss started. version=%s", VERSION)


# setting model to used
model = models[args.model]
logger.info("Selected VADR Model: %s", args.model)




# Run VADR
vadr_dir = os.path.join(work_dir, "vadr")
vadr_out_tbl_pass, vadr_out_fasta_pass, vadr_warnings = run_vadr(input_fasta, vadr_dir, model, cpu=1)

# Convert VADR result into .gbk
out_gbk_file = os.path.join(work_dir, "annotation.gbk")
tbl2genbank(vadr_out_tbl_pass, vadr_out_fasta_pass, out_gbk_file, model)



# prepare metadata file
metadata_file = copy_or_create_metadata_file(work_dir, args)


isolate, mss_file_prefix = get_isolate(metadata_file, args)

annotation_stats = check_annotation_stats(work_dir, vadr_dir, model)
# {'status': 'complete', 'total_length': 10735, 'model_length': 10735, 'query_coverage': '100.00%', 'qap_length': 0, 'cds_completeness': '1 / 0 / 1 [intact/partial/expected]'}
seq_status, number_of_sequence = annotation_stats["status"], annotation_stats["number_of_sequence"]
logger.info("Annotation stats: %s", annotation_stats)

organism = model.organism if hasattr(model, "organism") else None
    
update_metadata_file(metadata_file, seq_status=seq_status, number_of_sequence=number_of_sequence, mol_type=model.mol_type, organism=organism)


# Convert .gbk file and metadata file into MSS format.
mss = MSS2(out_gbk_file, metadata_file)
mss.convert(work_dir, mss_file_prefix)

warnings = [warning.to_tuple() for warning in vadr_warnings]
if vadr_warnings:
    for warning in vadr_warnings:
        logger.warning("VADR Warning: %s", warning)

out_report_file = os.path.join(work_dir, "dfv_report.json")
logger.info(f"Writing report json to {out_report_file}")
with open(out_report_file, "w") as f:
    json.dump({"annotation": annotation_stats, "warnings": warnings}, f, indent=4)

