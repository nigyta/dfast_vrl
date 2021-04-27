from Bio import SeqIO

def gbk2fasta(gbk_file, fasta_file):
    R = list(SeqIO.parse(gbk_file, "genbank"))
    with open(fasta_file, "w") as f:
        SeqIO.write(R, f, "fasta")
