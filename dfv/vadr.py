import subprocess

def run(input_file, output_dir, vadr_model="NC_045512"):
    cmd = f"v-annotate.pl --mxsize 16000 -s -r --nomisc --mkey {vadr_model} --lowsim5term 2 " + \
       f"--lowsim3term 2 --fstlowthr 0.0 --alt_fail lowscore,fsthicnf,fstlocnf,insertnn,deletinn " + \
       f"{input_file} {output_dir} -f"
    p = subprocess.run(cmd, shell=True, encoding="UTF-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    print(p.stdout)
    print(p.stderr)
    print(p.returncode)

