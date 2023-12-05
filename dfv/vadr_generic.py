import os
import sys
import subprocess

import logging

from .dfv_warnings import VADR_ANNOTATION_FAILED, VADR_ANNOTATION_FAILED_WITH_CRITICAL_WARNINGS, create_VADR_warning

logger = logging.getLogger(__name__)

# commands = {
#     "mpox": "v-annotate.pl --split --cpu {cpu} --glsearch --minimap2 -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --s_overhang 150 -f --mdir $VADRINSTALLDIR/vadr-models-mpxv-1.4.2-1 {fasta} {outdir}",
#     "sarscov2": "v-annotate.pl --split --cpu {cpu} --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn --mdir $VADRINSTALLDIR/vadr-models-sarscov2-1.3-2 {fasta} {outdir}"
#     }


def main(input_file, output_dir, model, cpu=1):

    logging.info("Running VADR...")
    cmd = model.command
    cmd = cmd.format(fasta=input_file, outdir=output_dir, cpu=cpu)

    warnings = []

    logger.debug(f"VADR command: {cmd}")
    p = subprocess.run(cmd, shell=True, encoding="UTF-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logger.info(f"{'='*30}  VADR stdout  {'='*30}\n{p.stdout}\n{'='*80}\n")
    if p.stderr:
        logger.info(f"{'='*30} VADR stderr {'='*30}\n{p.stderr}\n{'='*80}\n")
    logger.debug(f"VADR return status={p.returncode}")
    vadr_out_tbl_pass, vadr_out_fasta_pass, vadr_warnings = get_vadr_tbl_and_fasta(output_dir)
    warnings += vadr_warnings
    logger.info(f"VADR finished. Output files were generated in {output_dir}")
    return vadr_out_tbl_pass, vadr_out_fasta_pass, warnings

def get_vadr_tbl_and_fasta(output_dir):
    vadr_warnings = []
    vadr_prefix = os.path.basename(output_dir)  # same as the directory name
    vadr_out_tbl_pass = os.path.join(output_dir, vadr_prefix + ".vadr.pass.tbl")
    vadr_out_fasta_pass = os.path.join(output_dir, vadr_prefix + ".vadr.pass.fa")
    vadr_out_fail_list = os.path.join(output_dir, vadr_prefix + ".vadr.fail.list")
    if os.path.getsize(vadr_out_fasta_pass) == 0 and os.path.getsize(vadr_out_fail_list) > 0:
        vadr_out_tbl_pass = vadr_out_tbl_pass.replace(".pass.", ".fail.")
        vadr_out_fasta_pass = vadr_out_fasta_pass.replace(".pass.", ".fail.")
        logger.warning("VADR failed. Will use [%s] for downstream processes.", vadr_out_tbl_pass)
        vadr_warnings.append(VADR_ANNOTATION_FAILED_WITH_CRITICAL_WARNINGS)
        # logger.error("See fail file. [%s]", vadr_out_fasta_pass.replace(".pass.", ".fail."))
    alt_file = os.path.join(output_dir, vadr_prefix + ".vadr.alt")
    vadr_warnings += check_alert(alt_file)
    return vadr_out_tbl_pass, vadr_out_fasta_pass, vadr_warnings


def check_alert(alt_file):
    alerts = []
    for line in open(alt_file):
        cols = line.strip().split(maxsplit=9)
        if not cols[0].startswith("#"):
            fail = cols[7]
            alert_desc = cols[8]
            alert_detail = cols[9]
            target = f"{cols[1]}:{cols[3]}:{cols[4]}"
            VADR_WARNING = create_VADR_warning(fail, alert_desc, alert_detail)
            VADR_WARNING.add(target)
            alerts.append(VADR_WARNING)
    return alerts

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    # logger.setLevel(logging.DEBUG)
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    vadr_out_tbl_pass, vadr_out_fasta_pass = main(input_fasta, output_dir)
    print(vadr_out_tbl_pass, vadr_out_fasta_pass)