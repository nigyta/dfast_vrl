"""annotation.gbk から annotation.gff を作る

GFF3 への変換は dr_tools (seqrecords_to_gff) に任せる。DFAST の Web サービスが DDBJ Record から作る GFF3 と、
同じ変換で同じ形になる。以前はここに旧 DFAST の genbank2gff.py の写しがあった。dr_tools の変換は、それと同じ列を書き、
join などの複数の区間からなる位置を区間ごとの行にし、##sequence-region を足す
"""

from Bio import SeqIO
from dr_tools import seqrecords_to_gff


def gbk2gff(inputGenbankFileName, outputGffName):
    with open(inputGenbankFileName) as f:
        gff = seqrecords_to_gff(SeqIO.parse(f, "genbank"))
    with open(outputGffName, "w") as f:
        f.write(gff)
