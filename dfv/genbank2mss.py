#! /usr/bin/env python
# coding:utf8


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.InsdcIO import _insdc_location_string
try:
    from metadataUtil import Metadata
except ModuleNotFoundError:
    from .metadataUtil import Metadata

# from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import re
import logging

logger = logging.getLogger(__name__)

class MSS2():

    def __init__(self, genbankFileName, metadataFileName, dfc_annotation=False):

        def checkGap():
            S = set()
            records = SeqIO.parse(open(self.genbankFileName), "genbank")
            for r in records:
                for f in r.features:
                    S.add(f.type)
            return 'assembly_gap' in S

        self.genbankFileName = genbankFileName
        self.metadata = Metadata.readFromTSV(metadataFileName)
        # print self.metadata.fields["seqNames"].value
        self.locusTagPrefix = self.metadata.getValue("locusTagPrefix", "LOCUS")


        if checkGap():
            self.metadata.setValue("seqRank", "scaffold")
        else:
            self.metadata.setValue("seqRank", "contig")
        self.dfc_annotation = dfc_annotation


    def convert(self, outDir=".", prefix=None):
        if os.path.exists(outDir):
            if os.path.isfile(outDir):
                logger.error("Error. Could not create a directory.")
                return None
        else:
            os.makedirs(outDir)
        if prefix is None:
            prefix = ".".join(os.path.basename(self.genbankFileName).split(".")[:-1])

        logger.info("Converting GenBank file (%s) to DDBJ data submission files." % self.genbankFileName)
        fastaFileName = os.path.join(outDir, prefix + ".seq.fa")
        annotationFileName = os.path.join(outDir, prefix + ".annt.tsv")

        self.createSeqFile(fastaFileName)
        self.createTSV(annotationFileName)
        logger.info("Created DDBJ data submission files")


    def createSeqFile(self, outputFileName):
        def _toFasta(record):
            tmpRecord = SeqRecord(record.seq.lower(), id=record.name, description="")
            return tmpRecord.format("fasta")

        with open(outputFileName, "w") as f:
            records = SeqIO.parse(open(self.genbankFileName), "genbank")
            for record in records:
                f.write(_toFasta(record))
                f.write("//\n")
        logger.info("Created DDBJ seq file. :%s", outputFileName)

    def createTSV(self, outputFileName):

        def render_feature_and_qualifiers(feature, rec_length):
            ret_feature = []
            for qualifier, value in feature.qualifiers.items():
                if (feature.type == "assembly_gap" or feature.type == "gap") and qualifier == "estimated_length":
                    value = ["known"]
                if qualifier == "locus_tag" and self.locusTagPrefix:
                    value[0] = self.locusTagPrefix + "_" + value[0].split("_")[-1]
                elif qualifier == "transl_except":
                    transl_except_val = value[0]
                    if "complement" in transl_except_val:
                        value[0] = transl_except_val.replace("complement(", "").replace("),", ",")
                if qualifier not in ['translation', 'note', 'inference', 'EC_number', "ID"]:
                    ret_feature.append(["", "", "", qualifier, value[0]])
            for value in feature.qualifiers.get("inference", []):
                if "MBGD:" in value:
                    continue
                ret_feature.append(["", "", "", "inference", value])
            pseudo_flags = [None, None]  # internal stop codon, frameshift
            for value in feature.qualifiers.get("note", []):
                if value.startswith("MBGD"):
                    continue
                if value.endswith("]"):
                    if ", Eval:" in value:
                        continue
                if "_" in value and " " not in value:
                    prefix, tag = value.split("_", 1)
                    if tag.isdigit():
                        continue
                if "internal stop codon" in value:
                    pseudo_flags[0] = "internal stop codon"
                    continue
                if "frameshifted;" in value:
                    pseudo_flags[1] = "frameshift"
                    continue
                value = value.replace("  ", " ")
                ret_feature.append(["", "", "", "note", value])
            pseudo_flags = "/".join([x for x in pseudo_flags if x])
            if pseudo_flags:
                ret_feature.append(["", "", "", "note", "possible pseudo due to " + pseudo_flags])

            if len(ret_feature) == 0:  # Add row with empty values in col4 and col5 for features without any qualifier
                ret_feature.append(["", "", "", "", ""])
            ret_feature[0][1] = feature.type
            ret_feature[0][2] = _insdc_location_string(feature.location, rec_length)
            return ret_feature

        def renderOtherFeatures(record):
            ret = []
            rec_length = len(record)
            for feature in record.features:
                if feature.type in ['source', 'gene']:
                    continue
                else:
                    ret += render_feature_and_qualifiers(feature, rec_length)
            return ret

        ret = []
        ret += self.metadata.render4DFAST(dfc_annotation=self.dfc_annotation)
        records = SeqIO.parse(open(self.genbankFileName), "genbank")
        for i, record in enumerate(records):

            ret += self.metadata.renderSourceFeature(seqNum=i, seqRecord=record) 
            ret += renderOtherFeatures(record)


        out = "\n".join(["\t".join(cols) for cols in ret]) + "\n"
        with open(outputFileName, "w") as f:
            f.write(out)
        logger.info("Created DDBJ ann file. :%s", outputFileName)


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1:
        print(" Usage: python genbank2mss.py <genbankFileName> <metadataFileName> <outDir> <filePrefix>")
        print(" [OLD-MSS] Usage: python genbank2mss.py <GenbankFileName> <TemplateFileName> <metadataFileName> <OutputDirectory> <FilePrefix>")
        exit()
    else:
        genbankFileName = sys.argv[1]
        # templateFileName = sys.argv[2]
        # metadataFileName = sys.argv[3]
        metadataFileName = sys.argv[2]
        outputDirectory = sys.argv[3]
        filePrefix = sys.argv[4]
        # mss = MSS(genbankFileName, templateFileName, metadataFileName)
        # mss.convert(outputDirectory, filePrefix)

        # workDir = sys.argv[1]
        mss = MSS2(genbankFileName, metadataFileName)
        mss.convert(outputDirectory, filePrefix)
