# TO BE IMPLEMENTED

import os
import sys
import re
import json
from logging import getLogger
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation, ExactPosition, BeforePosition, AfterPosition
from Bio.Data.CodonTable import TranslationError
from Bio.SeqIO.InsdcIO import _insdc_location_string
from uuid import uuid4
from dataclasses import dataclass
from datetime import datetime
# if __name__ == '__main__':
#     from reference_models import get_reference_model
#     from dfv_warnings import INCOMPLETE_CDS_WARNING, INCOMPLETE_GENOME_WARNING, \
#         MORE_THAN_ONE_MODEL_USED, create_VADR_warning, VADR_ANNOTATION_FAILED, \
#         MISSING_FEATURES, PARTIAL_FEATURES, DUPLICATED_FEATURES, FRAGMENTED_FEATURES, INFO_NEARLY_COMPLETE_GENOME, INFO_DRAFT_GENOME
#     from vadr2ddbj import VadrFeature, MyLocation, get_reference_model
# else:
#     from .reference_models import get_reference_model
#     from .dfv_warnings import INCOMPLETE_CDS_WARNING, INCOMPLETE_GENOME_WARNING, \
#         MORE_THAN_ONE_MODEL_USED, create_VADR_warning, VADR_ANNOTATION_FAILED, \
#         MISSING_FEATURES, PARTIAL_FEATURES, DUPLICATED_FEATURES, FRAGMENTED_FEATURES, INFO_NEARLY_COMPLETE_GENOME, INFO_DRAFT_GENOME
#     from .vadr2ddbj import VadrFeature, MyLocation, get_reference_model

# Requires Biopython 1.78 and higher

logger = getLogger(__name__)


# This file seems to be not used anymore. Will be removed in the future. 2024.07.01


def gbk2faa(file_name):
    R = SeqIO.parse(file_name, "gbk")
    aa = []
    



class VADR2MSA:
    def __init__(self, fasta_file, vadr_dir, transl_table=1):
        
        self.seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.vadr_dir = vadr_dir
        self.transl_table = transl_table
        self.models = get_reference_model(name="sars_cov_2")
        self.report = {}
        self.warnings = []
        vadr_prefix = os.path.basename(self.vadr_dir)
        ftr_file = os.path.join(self.vadr_dir, vadr_prefix + ".vadr.ftr")
        self.set_features(ftr_file)
        today = datetime.now().strftime('%d-%b-%Y').upper()


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
                logger.warning("Translation Error. Will try again in a lenient mode.")
                logger.warning(str(seq_feature))
                translation = seq_feature.translate(seq_record.seq, cds=False, to_stop=True)
            seq_feature.qualifiers["translation"] = [str(translation)]


        features = list(VadrFeature.read(ftr_file))
        if len(features) == 0:
            self.warnings.append(VADR_ANNOTATION_FAILED)
            return

        appended_features = set()
        gene_features = [f for f in features if f.ftr_type == "gene"]
        other_features = [f for f in features if f.ftr_type != "gene"]
        for feature in other_features:
            model_feature = self.models[feature.model][feature.ftr_idx]
            model_feature.hits.append(feature.check_intactness())
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
            if "note" in model_feature.attrs:
                seq_feature.qualifiers.setdefault("note", []).append(model_feature.attrs["note"])
            if "exception" in model_feature.attrs and isinstance(seq_feature.location, CompoundLocation):
                assert model_feature.attrs["exception"] == "ribosomal slippage"
                seq_feature.qualifiers["ribosomal_slippage"] = [None]

            # add newly-created feature to seq_record
            self.seq_dict[feature.seq_name].features.append(seq_feature)
            appended_features.add(str(feature))




    # def to_gbk(self, output_file, with_feature_id=True):
    #     if with_feature_id:
    #         for r in self.seq_dict.values():
    #             for f in r.features:
    #                 f.qualifiers["ID"] = [f.id]
    #     SeqIO.write(self.seq_dict.values(), output_file, "genbank")

    def make_report(self, output_file, format="json"):
        pass





if __name__ == '__main__':
    
    ftr_file = "test/vadr/vadr.vadr.ftr"
    input_fasta = "../meta_vry_result_rc.fa"
    output_gbk = "test/test.vadr.gbk"

    # ftr_file = sys.argv[1]
    # vadr_dir = sys.argv[2]
    # output_gbk = sys.argv[3]



    # v2d = VADR2DDBJ(input_fasta, vadr_dir)
    # v2d.to_gbk(output_gbk)
    # v2d.make_report()
    # print(vsd.warnings)
"""
python /Users/tanizawa/dfast_web/app/scripts/genbank2mss.py
python genbank2mss.py <genbankFileName> <metadataFileName> <outDir> <filePrefix>
"""

