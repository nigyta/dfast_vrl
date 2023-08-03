import sys
import re
from dataclasses import dataclass, field
from typing import List

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition
from Bio.SeqIO.InsdcIO import _insdc_location_string
from Bio.Data.CodonTable import TranslationError
from .fix_product_for_mpox import fix_product_for_mpox
from .vadr2mss_config import Mpox
import logging

logger = logging.getLogger(__name__)

@dataclass
class FeatureData():
    """
    Class to store feature data in a .tbl file
    """

    type: str = None
    locations : List[List[str]] = field(default_factory=list)
    qualifiers : List[List[str]] = field(default_factory=list)
        
    
    def get_location(self):
        """
        Create FeatureLocation from self.locations
        If consisting of multiple locations, compound location will be created
        """

        def _get_location_for_each_seqment(start, end):
            '''
            Create FeatureLocation object for a given feature segment.

            for partial feature
             The “<” symbol always appears in column 1 and “>” always appears in column 2, regardless of the strandedness of the feature. 
            '''
            five_prime_partial = "<" in start            
            three_prime_partial = ">" in end

            int_start = int(start.replace("<", ""))
            int_end = int(end.replace(">", ""))
            if int_start < int_end:
                left, right, strand = int_start, int_end, 1
                left_pos = BeforePosition(left - 1) if five_prime_partial else ExactPosition(left - 1)
                right_pos = AfterPosition(right) if three_prime_partial else ExactPosition(right)

            elif int_start > int_end:
                left, right, strand = int_end, int_start, -1
                left_pos = BeforePosition(left - 1) if three_prime_partial else ExactPosition(left - 1)
                right_pos = AfterPosition(right) if five_prime_partial else ExactPosition(right)
            else:
                    logger.error(f"start position == end position, start={start}, end={end}\n")
                    logger.error("Sorry, currently not implemented. Aborted.")
                    exit(1)
            return FeatureLocation(left_pos, right_pos, strand)

        list_feature_locations = [_get_location_for_each_seqment(start, end) for start, end in self.locations]
        return sum(list_feature_locations)
    
    def to_seq_feature(self):
        """
        Convert into BioPython/SeqFeature object
        """
        location = self.get_location()
        qualifiers = {k: [v] for k, v in self.qualifiers}
        return SeqFeature(location=location, type=self.type, qualifiers=qualifiers)


def parse_tbl(tbl_file, fasta_file):
    """
    Generator function to yield SeqRecord with SeqFeatures converted from .tbl file
    """

    def _create_seq_features(rows):
        """
        convert .tbl feature rows into a list of SeqFeature
        """
        L = []
        feature_data = None
        for cols in rows:
            if len(cols) == 3:
                if feature_data:
                    L.append(feature_data)
                start, end, feature_type = cols
                feature_data = FeatureData(type=feature_type, locations=[[start, end]], qualifiers=[])
            elif len(cols) == 2:
                feature_data.locations.append(cols)  # cols should be [start, end]
            else:
                *_, q_key, q_value = cols
                if q_key == "protein_id":
                    continue  # skip unwanted protein_id
                feature_data.qualifiers.append([q_key, q_value])
        L.append(feature_data)
        
        # convert to SeqFeature
        return [feature_data.to_seq_feature() for feature_data in L]


    # Read FASTA files into dictionary {seq_id: SeqRecord}
    dict_seq = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    current_record = None
    rows = []
    for line in open(tbl_file):
        if line.startswith(">") and "\t" not in line:
            if current_record:
                features = _create_seq_features(rows)  # create SeqFeature list from rows
                current_record.features = features
                yield current_record
            seq_id = line.strip().split()[-1]
            current_record = dict_seq[seq_id]
            rows = []
        else:
            cols = line.strip("\n").split("\t")
            rows.append(cols)
    features = _create_seq_features(rows)  # create SeqFeature list from rows
    current_record.features = features
    yield current_record


def remove_gene_feature_and_add_gene_qualifier(r):
    """
    Remove gene features and add gene qualifiers to CDSs, mat_peptides, and RNAs.
    """
    target = ["CDS", "mat_peptide", "misc_feature", "stem_loop"]  # and *RNAs
    def _is_target(f_type):
        return f_type in target or f_type.endswith("RNA")
    
    gene_features = [f for f in r.features if f.type == "gene"]
    other_features = [f for f in r.features if f.type != "gene"]
    for feature in other_features:
        if not _is_target(feature.type):
            continue
        parent_gene_features = [f for f in gene_features if f.location.start <= feature.location.start and feature.location.end <= f.location.end and feature.location.strand == f.location.strand]
        if len(parent_gene_features) > 1:
            logger.warning(f"multiple gene features found for {str(feature.location)}")
        elif len(parent_gene_features) == 1:
            gene_qualifier = parent_gene_features[0].qualifiers.get("gene", [""])[0]
            if gene_qualifier:
                feature.qualifiers["gene"] = [gene_qualifier]
        else:
            logger.warning(f"No gene features found for {str(feature.location)}")
    r.features = other_features

def add_qualifiers_to_cds_feature(r, model):
    """
    Add translation qualifier to CDS features, as well as transl_table and codon_start.
    Product names will be modified to curated names (Mpox)
    """

    for feature in r.features:
        if feature.type != "CDS":
            continue
        
        # add tranl_table
        if "codon_start" not in feature.qualifiers:
            feature.qualifiers["codon_start"] = [1]
        transl_table = model.transl_table if hasattr(model, "transl_table") else 1
        if "transl_table" not in feature.qualifiers:
            feature.qualifiers["transl_table"] = [transl_table]

        try:
            translation = feature.translate(r.seq)
        except TranslationError as e:
            logger.warning(str(feature))
            logger.warning(f"length: {len(feature)}")
            logger.warning("translation error. Will try to translate again with 'cds=false'")
            print(e)
            translation = feature.translate(r.seq, cds=False).rstrip("*")
        feature.qualifiers["translation"] = [translation]

        # fix product
        if model.__name__ == "Mpox":  # isinstance(model, Mpox):
            prodict_original = feature.qualifiers.get("product", [""])[0]
            product_fixed = fix_product_for_mpox(prodict_original)
            if prodict_original != product_fixed:
                logger.debug(f"Fixing product name for Mpox annotation. {prodict_original} --> {product_fixed}")
            feature.qualifiers["product"] = [product_fixed]


def add_gaps(r, min_length=10):
    """
    add gap feature to SeqRecord.


    qualifiers = {"estimated_length": LENGTH}
    """
    pat_gap = re.compile(r"[Nn]+")
    seq = str(r.seq)
    matches = pat_gap.finditer(seq)

    for m in matches:
        start, end = m.span()
        length = end - start
        if length >= min_length:
            qualifiers = {"estimated_length": [length]}
            r.features.append(SeqFeature(FeatureLocation(start, end, 1), type="gap", qualifiers=qualifiers))


def sort_features(r):
    """Sort features by start position"""
    r.features = sorted(r.features, key=lambda x: x.location.start)


def tbl2genbank(tbl_file, fasta_file, out_gbk_file, model):
    logger.info(f"Converting {tbl_file} and {fasta_file} into {out_gbk_file}")
    R = [r for r in parse_tbl(tbl_file, fasta_file)]
    for r in R:
        r.annotations["molecule_type"] = model.mol_type
        remove_gene_feature_and_add_gene_qualifier(r)
        add_qualifiers_to_cds_feature(r, model)  # model is provided to specify transl_table number.
        add_gaps(r)
        sort_features(r)
    logger.info("Writing .gbk file [%s]", out_gbk_file)
    SeqIO.write(R, open(out_gbk_file, "w"), "genbank")



if __name__ == "__main__":
    import argparse
    # parser = argparse.ArgumentParser(prog="tbl2genbank.py")
    # parser.add_argument("-t", "--tbl", metavar="FILE", help="Input .tbl file", required=True)
    # parser.add_argument("-f", "--fasta", metavar="FILE", help="Input FASTA file", required=True)
    # parser.add_argument("-o", "--output", metavar="FILE", help="Output GenBank file", required=True)
    # args = parser.parse_args()
    # tbl2genbank(args.tbl, args.fasta, args.output, model)
