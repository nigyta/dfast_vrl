
class Mpox():
    command = "v-annotate.pl --split --cpu {cpu} --glsearch --minimap2 -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --s_overhang 150 -f --mdir $VADRINSTALLDIR/vadr-models-mpxv-$VADR_MPXV_MODELS_VERSION {fasta} {outdir}"
    mol_type = "DNA"
    length = 200000

class Sarscov2():
    command = "v-annotate.pl --split --cpu {cpu} --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn -f --mdir $VADRINSTALLDIR/vadr-models-sarscov2-$VADR_SCOV2_MODELS_VERSION {fasta} {outdir}"
    mol_type = "RNA"
    length = 30000

models = {
    "mpox": Mpox,
    "sarscov2": Sarscov2
}

if __name__ == "__main__":
    print(models["mpox"].command)