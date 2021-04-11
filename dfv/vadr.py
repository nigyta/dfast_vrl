import os
import subprocess

def run(input_file, output_dir, vadr_model="NC_045512"):
    cmd = f"v-annotate.pl --mxsize 16000 -s -r --nomisc --mkey {vadr_model} --lowsim5term 2 " + \
       f"--lowsim3term 2 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf,insertnn,deletinn " + \
       f"{input_file} {output_dir} -f"
    p = subprocess.run(cmd, shell=True, encoding="UTF-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(p.stdout)
    print(p.stderr)
    print(p.returncode)
    vadr_result_fasta, vadr_result_ftr = get_vadr_result(output_dir)
    return vadr_result_fasta, vadr_result_ftr

def get_vadr_result(output_dir):
    vadr_prefix = os.path.basename(output_dir)
    vadr_result_ftr = os.path.join(output_dir, vadr_prefix + ".vadr.ftr")
    vadr_out_fasta_pass = os.path.join(output_dir, vadr_prefix + ".vadr.pass.fa")
    vadr_out_fasta_fail = os.path.join(output_dir, vadr_prefix + ".vadr.fail.fa")
    vadr_result_fasta = os.path.join(output_dir, vadr_prefix + ".result.fa")
    with open(vadr_result_fasta, "w") as f:
        f.write(open(vadr_out_fasta_pass).read())
        f.write(open(vadr_out_fasta_fail).read())
    return vadr_result_fasta, vadr_result_ftr
