from .fix_blank_in_products import fix_blank_in_products

from Bio import SeqIO
# import urllib

def getPhase(feature):
    if feature.type == "CDS":
        return str(int(feature.qualifiers.get("codon_start", [1])[0]) - 1)
    else:
        return "."

def getSource(feature):
    q = feature.qualifiers.get("inference", [])
    for val in q:
        if not val.startswith("similar to") or not val.startswith("nucleotide motif:Rfam"):
            return val.split(":")[-2] + ":" + val.split(":")[-1]
    return "dfast"

# def removeKeys(D, keys):
#     for key in keys:
#         if key in D:
#             del D[key]

def quote_attribute(x):
    return str(x).replace("%", "%25").replace(",", "%2C").replace("=", "%3D").replace(";", "%3B").replace("&", "%26").replace("\t", "%09")

def createAttribute(feature):
    q = feature.qualifiers
    # removeKeys(q, ["codon_start", "transl_table", "#inference"])
    ignored_key = ["codon_start", "transl_table"]
    # attributes = "; ".join([key + "=" + quote_attribute("; ".join(value)) for key, value in q.items() if key not in ignored_key])
    attributes = ";".join([key + "=" + ",".join(map(quote_attribute, value))
                               for key, value in q.items() if key not in ignored_key])

    # attributes = ";".join([key + "=" + quote_attribute(",".join(value)) for key, value in q.items()])

    locusTag = q.get("locus_tag", [None])[0]
    if locusTag:
        return "ID=" + locusTag + ";" + attributes
    else:
        return attributes

def gbk2gff(inputGenbankFileName, outputGffName):
    stringBuffer = "##gff-version 3\n"
    records = [record for record in SeqIO.parse(open(inputGenbankFileName), "genbank")]
    for record in records:
        for feature in record.features:
            if feature.type not in ["gene", "source", "assembly_gap"]:

                seqid = record.name  # col0
                inferenceSource = getSource(feature)  # col1
                featureType = feature.type  # col2
                start_pos = str(feature.location.start + 1)  # col3
                end_pos = str(int(feature.location.end))  # col4
                score = "."  # col5
                strand = "-" if feature.location.strand == -1 else "+"  # col6
                phase = getPhase(feature)
                attributes = createAttribute(feature)
                line = "\t".join([seqid, inferenceSource, featureType, start_pos, end_pos, score, strand, phase, attributes]) + "\n"
                stringBuffer += line
    # stringBuffer += "###\n"
    stringBuffer += "##FASTA\n"
    for record in records:
        seqID = record.name
        stringBuffer += ">" + seqID + "\n" + str(record.seq) + "\n"
    stringBuffer = fix_blank_in_products(stringBuffer)
    with open(outputGffName, "w") as f:
        f.write(stringBuffer)

if __name__ == '__main__':
    fileName = "annotation.gbk"
    outFile = "test.gff"
    gbk2gff(fileName, outFile)
