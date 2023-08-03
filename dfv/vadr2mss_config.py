
class Mpox():
    command = "v-annotate.pl --split --cpu {cpu} --glsearch --minimap2 -s -r --nomisc --mkey mpxv --r_lowsimok --r_lowsimxd 100 --r_lowsimxl 2000 --alt_pass discontn,dupregin --s_overhang 150 -f --mdir $VADRINSTALLDIR/vadr-models-mpxv-$VADR_MPXV_MODELS_VERSION {fasta} {outdir}"
    mol_type = "DNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-mpxv-$VADR_MPXV_MODELS_VERSION/mpxv.minfo"
    organism = "Monkeypox virus"

class Sarscov2():
    command = "v-annotate.pl --split --cpu {cpu} --glsearch -s -r --nomisc --mkey sarscov2 --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn -f --mdir $VADRINSTALLDIR/vadr-models-sarscov2-$VADR_SCOV2_MODELS_VERSION {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-sarscov2-$VADR_SCOV2_MODELS_VERSION/sarscov2.minfo"

class Corona():
    command = "v-annotate.pl --split --cpu {cpu} --glsearch -s -r --nomisc --mkey corona --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn -f --mdir $VADRINSTALLDIR/vadr-models-corona-$VADR_CORONA_MODELS_VERSION {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-corona-$VADR_CORONA_MODELS_VERSION/corona.minfo"

class RSV():
    command = "v-annotate.pl --split --cpu {cpu} -r --mkey rsv --xnocomp -f --mdir $VADRINSTALLDIR/vadr-models-rsv-$VADR_RSV_MODELS_VERSION {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-rsv-$VADR_RSV_MODELS_VERSION/rsv.minfo"

class Noro():
    command = "v-annotate.pl --split --cpu {cpu} --group Norovirus --mkey calici --nomisc --noprotid -f -mdir $VADRINSTALLDIR/vadr-models-calici {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-calici/noro.minfo"

class Calisi():
    command = "v-annotate.pl --split --cpu {cpu} --mkey calici --nomisc --noprotid -f -mdir $VADRINSTALLDIR/vadr-models-calici {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-calisi/calisi.minfo"

class Dengue():
    command = "v-annotate.pl --split --cpu {cpu} --group Dengue --nomisc --noprotid --mkey flavi -f -mdir $VADRINSTALLDIR/vadr-models-flavi {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-flavi/dengue.minfo"

class Flavi():
    command = "v-annotate.pl --split --cpu {cpu} --nomisc --noprotid --mkey flavi -f -mdir $VADRINSTALLDIR/vadr-models-flavi {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-flavi/flavi.minfo"

class COX1():
    command = "v-annotate.pl --split --cpu {cpu} --xmaxdel 3 --xmaxins 3 --xlongest --mkey cox1 --alt_pass lowcovrg --alt_fail fstlocfi,fstlocf5,fstlocf3 --fstminnti 5 --fstminnt5 5 --fstminnt3 5 --nomisc --noprotid --xsub $VADRINSTALLDIR/vadr-models-cox1-$VADR_COX1_MODELS_VERSION/cox1.phy.xsub --mdir $VADRINSTALLDIR/vadr-models-cox1-$VADR_COX1_MODELS_VERSION {fasta} {outdir}"
    mol_type = "RNA"
    minfo_file = "$VADRINSTALLDIR/vadr-models-cox1-$VADR_COX1_MODELS_VERSION/cox1.minfo"







models = {
    "mpox": Mpox,
    "sarscov2": Sarscov2,
    "corona": Corona,
    "RSV": RSV,
    "Noro": Noro,
    "Calici": Calisi,
    "Dengue": Dengue,
    "Flavi": Flavi,
    "COX1": COX1
}

if __name__ == "__main__":
    print(models["mpox"].command)
