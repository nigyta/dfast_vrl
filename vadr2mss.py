#!/usr/bin/env python

import os
import sys
import argparse

parser = argparse.ArgumentParser(prog="vadr2mss.py")
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
        metavar="MODEL"
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
from dfv.common import update_metadata_file, get_isolate, copy_or_create_metadata_file, get_logger
from dfv.vadr2mss_config import models

# cheking input fasta
input_fasta = args.input
if not os.path.exists(input_fasta):
    sys.stderr.write(f"vadr2mss.py: error: Specified FASTA file does not exist. [{input_fasta}]\n")
    exit(1)


# setting output directory
work_dir = args.out_dir
if os.path.exists(work_dir) and not args.force:
    if not args.force:
        sys.stderr.write("vadr2mss.py: error: Output directory already exists. Use '--force' to overwrite.\n")
        exit(1)
else:
    os.makedirs(work_dir, exist_ok=True)

# setting model to used
model = models[args.model]

logger = get_logger(name=__name__, debug=args.debug, work_dir=work_dir)

# Run VADR
vadr_out_tbl_pass, vadr_out_fasta_pass = run_vadr(input_fasta, work_dir, model, cpu=1)

# Convert VADR result into .gbk
out_gbk_file = os.path.join(work_dir, "annotation.gbk")
tbl2genbank(vadr_out_tbl_pass, vadr_out_fasta_pass, out_gbk_file, model)



# prepare metadata file
metadata_file = copy_or_create_metadata_file(work_dir, args)


isolate, mss_file_prefix = get_isolate(metadata_file, args)
# TODO: implement function to check status and num of sequence 
# TDDO: refer to .minfo file for CDSs, total length
update_metadata_file(metadata_file, seq_status="nearly complete", number_of_sequence=1, mol_type=model.mol_type)

# Convert .gbk file and metadata file into MSS format.
mss = MSS2(out_gbk_file, metadata_file)
mss.convert(work_dir, mss_file_prefix)
