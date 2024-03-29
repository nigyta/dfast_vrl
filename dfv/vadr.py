import os
import sys
import subprocess

import logging

logger = logging.getLogger(__name__)


def run(input_file, output_dir, cpu=1):
    # v-annotate.pl --split --cpu 4 --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5term 2 --lowsim3term 2 
    # --alt_fail lowscore,fstukcnf,insertnn,deletinn [--mdir vadr1.2/vadr-models-sarscov2-1.2-2] input.fasta out_dir
    logging.info("Running VADR...")
    # for vadr1.2
    # cmd = f"v-annotate.pl --split --cpu {cpu} --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5term 2 --lowsim3term 2 " + \
    #    f"--alt_fail lowscore,fstukcnf,insertnn,deletinn {input_file} {output_dir} -f" # --out_stk --out_afa --out_rpafa  --out_rpstk --out_allfasta "

    # For DFAST_VRL, model files are under $VADRINSTALLDIR. For VADR2MSS, under 
    # cmd = f"v-annotate.pl --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 " + \
    #    f"--alt_fail lowscore,insertnn,deletinn -f --mdir $VADRMODELDIR/vadr-models-sarscov2-$VADR_SCOV2_MODELS_VERSION {input_file} {output_dir}" # --out_stk --out_afa --out_rpafa  --out_rpstk --out_allfasta "

    # cmd = f"v-annotate.pl --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 " + \
    #    f"--alt_fail lowscore,insertnn,deletinn -f --mdir $VADRINSTALLDIR/vadr-models-sarscov2-$VADR_SCOV2_MODELS_VERSION {input_file} {output_dir}" # --out_stk --out_afa --out_rpafa  --out_rpstk --out_allfasta "

    # 20231028modified
    cmd = f"v-annotate.pl --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 " + \
       f"--alt_fail lowscore,insertnn,deletinn -f --mdir $VADRMODELDIR/vadr-models-sarscov2-$VADR_SCOV2_MODELS_VERSION {input_file} {output_dir}" # --out_stk --out_afa --out_rpafa  --out_rpstk --out_allfasta "

    logger.debug(f"VADR command: {cmd}")
    p = subprocess.run(cmd, shell=True, encoding="UTF-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    logger.info(f"{'='*30}  VADR stdout  {'='*30}\n{p.stdout}\n{'='*80}\n")
    if p.stderr:
        logger.info(f"{'='*30} VADR stderr {'='*30}\n{p.stderr}\n{'='*80}\n")
    logger.debug(f"VADR return status={p.returncode}")
    vadr_result_fasta, vadr_result_ftr = get_vadr_result(output_dir)
    logger.info(f"VADR finished. Output files were generated in {output_dir}")
    return vadr_result_fasta, vadr_result_ftr

def get_vadr_result(output_dir):
    vadr_prefix = os.path.basename(output_dir)  # same as the directory name
    vadr_result_ftr = os.path.join(output_dir, vadr_prefix + ".vadr.ftr")
    vadr_out_fasta_pass = os.path.join(output_dir, vadr_prefix + ".vadr.pass.fa")
    vadr_out_fasta_fail = os.path.join(output_dir, vadr_prefix + ".vadr.fail.fa")
    vadr_result_fasta = os.path.join(output_dir, vadr_prefix + ".result.fa")
    with open(vadr_result_fasta, "w") as f:
        f.write(open(vadr_out_fasta_pass).read())
        f.write(open(vadr_out_fasta_fail).read())
    return vadr_result_fasta, vadr_result_ftr

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    # logger.setLevel(logging.DEBUG)
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    vadr_result_fasta, vadr_result_ftr = run(input_fasta, output_dir)