import os
import sys
import re
from logging import getLogger
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition
from Bio.Data.CodonTable import TranslationError
from Bio.SeqIO.InsdcIO import _insdc_location_string

from dataclasses import dataclass
from datetime import datetime
# Requires Biopython 1.78 and higher


logger = getLogger(__name__)

@dataclass
class VadrFeature:
    """
    Class for VADR feature table file.
    The result file (.vadr.ftr) is a table consisting of 25 columns.
    """
    idx: str
    seq_name: str
    seq_len: int  # col2
    pass_or_fail: str
    model: str
    ftr_type: str
    ftr_name: str  # replace _ with " "  
    ftr_len: int  # col7
    ftr_idx: int  # col8
    par_idx: str
    strand: str
    n_from: int  # col11
    n_to: int  # col12
    n_instp: int # col13, int or None
    trc: str
    five_prime_n: int  # col15 
    three_prime_n: int  # col16
    p_from: int  # col17, int or None
    p_to: int  # col18, int or None
    p_instp: int  # col19, int or None
    p_sc: int  # col20, int or None
    nsa: int  # col21
    nsn: int  # col22
    seq_coords: str
    model_coords: str
    ftr_alerts: str  # col25 str or None

    @staticmethod
    def read(file_name):
        for line in open(file_name):
            if line.startswith("#"):
                continue
            cols = line.strip("\n").split()
            assert len(cols) == 26
            for idx in [2, 7, 8, 11, 12, 15, 16, 21, 22]:
                cols[idx] = int(cols[idx])
            for idx in [13, 17, 18, 19, 20]:
                cols[idx] = None if cols[idx] == "-" else int(cols[idx])
            cols[25] = None if cols[25] == "-" else cols[25]
            cols[6] = cols[6].replace("_", " ")
            vadr_feature = VadrFeature(*cols)
            yield vadr_feature

    def __str__(self):
        return f"<{self.ftr_type} {self.seq_name}:{self.n_from}-{self.n_to}>"

    def get_feature_location(self):
        coords = self.seq_coords.split(",")
        locations = [MyLocation(cord) for cord in coords]
        if self.trc == "5'" or self.five_prime_n > 0:
            locations[0].set_left_partial(self.five_prime_n)
        if self.trc == "3'" or self.three_prime_n > 0:
            locations[-1].set_right_partial(self.three_prime_n)
        locations = [my_loc.get_feature_location() for my_loc in locations]
        return sum(locations)  # ret is CompoundLocation for joined locations

    def get_seq_feature(self):
        location = self.get_feature_location()
        if self.ftr_type == "stem_loop":
            qualifiers = {}
        else:
            qualifiers = {"product": [self.ftr_name]}
        if self.ftr_type == "CDS":
            # Add extra qualifiers for CDS
            if isinstance(location, CompoundLocation):
                qualifiers["ribosomal_slippage"] = [None]
            # set codon start (assuming strand is +), TODO: convert to misc_feature if both ends are ambiguous
            if isinstance(location.start, BeforePosition) and isinstance(location.end, ExactPosition):
                length = location.end - location.start
                codon_start = (length % 3) + 1
                qualifiers["codon_start"] = [codon_start]
            else:
                qualifiers["codon_start"] = [1]
        additional_qualifier = additional_qualifiers.get((self.ftr_type, self.ftr_name))
        if additional_qualifier:
            qualifiers.update(additional_qualifier)
        feature = SeqFeature(location=location, type=self.ftr_type, strand=location.strand, id=self.idx, qualifiers=qualifiers)
        return feature


class VADR2DDBJ:
    def __init__(self, fasta_file, ftr_file, metadata=None, 
                 organism="Severe acute respiratory syndrome coronavirus 2", isolate="(isolate)", complete=True,
                 mol_type="genomic RNA", transl_table=1, topology="linear", linkage_evidence="align genus"):
        
        self.seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.metadata = metadata
        self.organism = organism
        self.isolate = isolate
        self.mol_type = mol_type
        self.transl_table = transl_table
        self.complete = complete
        self.linkage_evidence = linkage_evidence
        self.set_source_feature()
        self.set_gap_features()
        self.set_features(ftr_file)

        today = datetime.now().strftime('%d-%b-%Y').upper()
        annotations = {
            "organism": self.organism, "isolate": self.isolate,
            "complete": self.complete, "date": today, "topology": topology,
            "data_file_division": "VRL", "molecule_type": self.mol_type,
            # "sequence_version": 1, "data_file_division": "BCT",
            # "taxonomy":['Eukaryota', 'Viridiplantae', 'Streptophyta', 'Embryophyta', 'Paphiopedilum'],
        }
        for seq_record in self.seq_dict.values():
            seq_record.annotations = annotations

    def set_source_feature(self):
        for seq_record in self.seq_dict.values():
            location = FeatureLocation(0, len(seq_record))
            qualifiers = {
                "organism": [self.organism],
                "isolate": [self.isolate],
                "mol_type": [self.mol_type],
            }
            source_feature = SeqFeature(location=location, type="source", strand=location.strand, qualifiers=qualifiers)
            seq_record.features.append(source_feature)

    def set_gap_features(self, len_cutoff=10):
        for record in self.seq_dict.values():
            startPosition = 0
            seq = str(record.seq).upper()
            pat = "(" + "N" * len_cutoff + "+)"
            for fragment in re.split(pat, seq):
                endPosition = startPosition + len(fragment)
                if fragment.startswith("N"):
                    qualifiers = {"estimated_length": [len(fragment)], "gap_type": ["within scaffold"],
                                  "linkage_evidence": [self.linkage_evidence]}
                    location = FeatureLocation(startPosition, endPosition, strand=1)
                    feature = SeqFeature(location, type="assembly_gap", qualifiers=qualifiers)

                    assert str(feature.extract(record).seq).upper() == fragment
                    record.features.append(feature)
                startPosition = endPosition

    def set_features(self, ftr_file):
        """
        Read Feature file from VADR, and set annotated features to SeqRecord.
        """
        def _get_gene(feature, gene_features):
            # get /gene qualifier
            gene_features_filtered = [gf for gf in gene_features if gf.n_from <= feature.n_from <= feature.n_to <= gf.n_to and gf.seq_name == feature.seq_name]
            if len(gene_features_filtered) != 1:
                logger.warning(f"Cannot find gene feature for {feature.idx}:{feature.ftr_name}")
                return None
            else:
                return gene_features_filtered[0].ftr_name

        def _set_translation(seq_feature, seq_record):
            try:
                translation = seq_feature.translate(seq_record.seq)
            except TranslationError as e:
                print("[Warning!] Translation Error.")
                print(seq_feature)
                translation = seq_feature.translate(seq_record.seq, cds=False, to_stop=True)
            seq_feature.qualifiers["translation"] = [str(translation)]

        def _set_five_prime_utr(features):
            assert features[0].type == "source"

            cds_features = [f for f in features if f.type == "CDS"]
            first_feature = cds_features[0]
            if first_feature.qualifiers.get("gene", [""])[0] == "ORF1ab" and first_feature.location.start > 1:
                utr_feature = SeqFeature(FeatureLocation(0, first_feature.location.start, 1), type="5'UTR")
                features.insert(1, utr_feature)  # insert UTR next to the source feature


        def _set_three_prime_utr(features, seq_length):
            cds_features = [f for f in features if f.type == "CDS"]
            last_feature = cds_features[-1]
            if last_feature.qualifiers.get("gene", [""])[0] == "ORF10" and last_feature.location.end < seq_length:
                utr_feature = SeqFeature(FeatureLocation(last_feature.location.end, seq_length, 1), type="3'UTR")
                features.append(utr_feature)


        features = list(VadrFeature.read(ftr_file))
        appended_features = set()
        gene_features = [f for f in features if f.ftr_type == "gene"]
        other_features = [f for f in features if f.ftr_type != "gene"]
        for feature in other_features:
            if str(feature) in appended_features:
                logger.info(f"Skipping duplicated features {feature}")
                continue
            seq_feature = feature.get_seq_feature()  # Create Biopython SeqFeature object from MyFeature
            gene = _get_gene(feature, gene_features)
            if gene:  # assing gene if possible
                seq_feature.qualifiers["gene"] = [gene]
            if seq_feature.type == "CDS":
                seq_feature.qualifiers["transl_table"] = [self.transl_table]
                _set_translation(seq_feature, self.seq_dict[feature.seq_name])

            # add newly-created feature to seq_record
            self.seq_dict[feature.seq_name].features.append(seq_feature)
            appended_features.add(str(feature))

        for seq_record in self.seq_dict.values():  # iterate seq records
            seq_record.features.sort(key=lambda f: f.location.start)
            _set_five_prime_utr(seq_record.features)
            _set_three_prime_utr(seq_record.features, len(seq_record))

    def to_gbk(self, output_file):
        SeqIO.write(self.seq_dict.values(), output_file, "genbank")


@dataclass
class MyLocation:
    """
    Class to handle coords in the VADR result.
    Coordinate is 1-based.
    e.g. 26467..27135:+
    """
    left: int
    right: int
    strand: int
    left_partial: bool
    right_partial: bool

    def __init__(self, coords):
        left_and_right, strand = coords.split(":")
        strand = 1 if strand == "+" else -1
        left, right = map(int, left_and_right.split(".."))
        self.left, self.right, self.strand = left, right, strand
        self.left_partial, self.right_partial = False, False

    def set_left_partial(self, five_prime_n):
        # Strand is always +, assuming singe-stranded genome of SARS-Cov2. 
        self.left = self.left + five_prime_n
        self.left_partial = True

    def set_right_partial(self, three_prime_n):
        # Strand is always +, assuming singe-stranded genome of SARS-Cov2. 
        self.right = self.right - three_prime_n
        self.right_partial = True

    def get_feature_location(self):
        # Coordinate is 0-based in Biopython object
        start = BeforePosition(self.left - 1) if self.left_partial else ExactPosition(self.left - 1)
        end = AfterPosition(self.right) if self.right_partial else ExactPosition(self.right)
        return FeatureLocation(start, end, strand=self.strand)

additional_qualifiers = {
    ("stem_loop","stem loop.1"): {"note": ["Coronavirus 3' stem-loop II-like motif (s2m)"]},
    ("stem_loop","stem loop.2"): {"note": ["Coronavirus 3' UTR pseudoknot stem-loop 1"]},
    ("stem_loop","stem loop.3"): {"note": ["Coronavirus 3' UTR pseudoknot stem-loop 2"]},
    ("stem_loop","stem loop.4"): {"note": ["Coronavirus frameshifting stimulation element stem-loop 1"]},
    ("stem_loop","stem loop.5"): {"note": ["Coronavirus frameshifting stimulation element stem-loop 2"]},
}


if __name__ == '__main__':
    
    # ftr_file = "VADR_MW850590/VADR_MW850590.vadr.ftr"
    # input_fasta = "VADR_MW850590/VADR_MW850590.vadr.pass.fa"
    ftr_file = "vadr_with_truncated/1_vadr_finished.vadr.ftr"
    input_fasta = "vadr_with_truncated/1_vadr_finished.vadr.pass.fa"
    output_gbk = "test.vadr.gbk"
    # for feature in VadrFeature.read(input_file):
    #     print(feature)
    #     print(feature.get_feature_location())
    v2d = VADR2DDBJ(input_fasta, ftr_file)
    v2d.to_gbk(output_gbk)


"""
python /Users/tanizawa/dfast_web/app/scripts/genbank2mss.py
python genbank2mss.py <genbankFileName> <metadataFileName> <outDir> <filePrefix>
"""