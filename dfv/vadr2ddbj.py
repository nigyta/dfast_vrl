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
if __name__ == '__main__':
    from reference_models import get_reference_model
    from dfv_warnings import INCOMPLETE_CDS_WARNING, INCOMPLETE_GENOME_WARNING, \
        MORE_THAN_ONE_MODEL_USED, create_VADR_warning, VADR_ANNOTATION_FAILED, \
        MISSING_FEATURES, PARTIAL_FEATURES, DUPLICATED_FEATURES, FRAGMENTED_FEATURES, INFO_NEARLY_COMPLETE_GENOME, INFO_DRAFT_GENOME
else:
    from .reference_models import get_reference_model
    from .dfv_warnings import INCOMPLETE_CDS_WARNING, INCOMPLETE_GENOME_WARNING, \
        MORE_THAN_ONE_MODEL_USED, create_VADR_warning, VADR_ANNOTATION_FAILED, \
        MISSING_FEATURES, PARTIAL_FEATURES, DUPLICATED_FEATURES, FRAGMENTED_FEATURES, INFO_NEARLY_COMPLETE_GENOME, INFO_DRAFT_GENOME

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
        # ribosomal_slippageを含むmat_peptideでCDSの前半部分がすべてNになっている場合の例外処理 -----
        if len(locations) > 1:
            first_cds, second_cds = locations[0], locations[1]
            len_of_1st_cds = first_cds.right - first_cds.left + 1
            offset_5_prime_n = self.five_prime_n - len_of_1st_cds
            if offset_5_prime_n >= 0:
                logger.warning(f"First part of ribosomal_slippage is removed as it is located within a gap [{self.ftr_type}:{self.ftr_name}]")
                self.five_prime_n = offset_5_prime_n
                locations = [second_cds]
            # print("offset_5_prime_n", offset_5_prime_n)
        # -----
        if "5'" in self.trc or self.five_prime_n > 0:   # self.trc can be either of 5' 3' 5'&3'
            locations[0].set_left_partial(self.five_prime_n)
        if "3'" in self.trc or self.three_prime_n > 0:
            locations[-1].set_right_partial(self.three_prime_n)
        locations = [my_loc.get_feature_location() for my_loc in locations]
        return sum(locations)  # ret is CompoundLocation for joined locations

    def get_seq_feature(self):
        if self.ftr_len == self.five_prime_n == self.three_prime_n:  # feature is within a gap
            return None
        location = self.get_feature_location()
        if self.ftr_type == "stem_loop":
            qualifiers = {}
        else:
            qualifiers = {"product": [self.ftr_name]}
        if self.ftr_type == "CDS":
            if self.n_instp:
                if self.n_instp < self.n_to:
                    self.ftr_type = "misc_feature"
                    del qualifiers["product"]
                    note = [f"similar to {self.ftr_name}", f"internal stop at {self.n_instp}"]
                else:
                    note = [f"ambiguous 3'-end"]
                qualifiers["note"] = note
            else:
                # set codon start (assuming strand is +), TODO: convert to misc_feature if both ends are ambiguous
                if isinstance(location.start, BeforePosition) and isinstance(location.end, ExactPosition):
                    # length2 = location.end - location.start
                    length = len(location)
                    codon_start = (length % 3) + 1
                    # print("LENGTH, LENGTH2, CODON_START", length, length2, codon_start)
                    # print(location)
                    qualifiers["codon_start"] = [codon_start]
                else:
                    qualifiers["codon_start"] = [1]
        feature = SeqFeature(location=location, type=self.ftr_type, strand=location.strand, id=self.idx, qualifiers=qualifiers)
        return feature

    def check_intactness(self):
        if self.trc != "no" or self.five_prime_n > 0 or self.three_prime_n > 0:
            return "partial"
        else:
            return "intact"

class VADR2DDBJ:
    def __init__(self, fasta_file, vadr_dir, metadata=None, 
                 organism="Severe acute respiratory syndrome coronavirus 2", isolate=None, complete=True,
                 mol_type="genomic RNA", transl_table=1, topology="linear", linkage_evidence="align genus"):
        
        self.seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        self.vadr_dir = vadr_dir
        self.metadata = metadata
        self.organism = organism
        self.isolate = isolate or "(isolate)"
        self.mol_type = mol_type
        self.transl_table = transl_table
        self.complete = complete
        self.linkage_evidence = linkage_evidence
        self.models = get_reference_model(name="sars_cov_2")
        self.report = {}
        self.warnings = []
        self.set_source_feature()
        self.set_gap_features()
        vadr_prefix = os.path.basename(self.vadr_dir)
        ftr_file = os.path.join(self.vadr_dir, vadr_prefix + ".vadr.ftr")
        self.set_features(ftr_file)
        self.modify_incomplete_cds()
        self.change_mat_peptide_to_misc()
        self.check_completeness()
        self.add_UTRs_to_complete_genome()
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
            source_feature = SeqFeature(id=uuid4(), location=location, type="source", strand=location.strand, qualifiers=qualifiers)
            seq_record.features.append(source_feature)

    def set_gap_features(self, len_cutoff=10):
        num_assembly_gap = 0
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
                    feature = SeqFeature(location, id=uuid4(), type="assembly_gap", qualifiers=qualifiers)

                    assert str(feature.extract(record).seq).upper() == fragment
                    record.features.append(feature)
                    num_assembly_gap += 1
                startPosition = endPosition
        self.report["num_assembly_gap"] = num_assembly_gap

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
            if seq_feature is None:
                logger.warning(f"{feature.ftr_type}:{feature.n_from}-{feature.n_to}({feature.ftr_name}) is located within an assembly_gap.")
                continue
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



    def add_UTRs_to_complete_genome(self):

        def _set_five_prime_utr(features):
            assert features[0].type == "source"
            cds_features = [f for f in features if f.type == "CDS"]
            first_feature = cds_features[0]
            if first_feature.qualifiers.get("gene", [""])[0] == "ORF1ab" and first_feature.location.start > 1:
                utr_feature = SeqFeature(FeatureLocation(0, first_feature.location.start, 1), id=uuid4(), type="5'UTR")
                features.insert(1, utr_feature)  # insert UTR next to the source feature


        def _set_three_prime_utr(features, seq_length):
            cds_features = [f for f in features if f.type == "CDS"]
            last_feature = cds_features[-1]
            if last_feature.qualifiers.get("gene", [""])[0] == "ORF10" and last_feature.location.end < seq_length:
                utr_feature = SeqFeature(FeatureLocation(last_feature.location.end, seq_length, 1), id=uuid4(), type="3'UTR")
                features.append(utr_feature)

        seq_status = self.report.get("seq_status", "")
        for seq_record in self.seq_dict.values():  # iterate seq records
            seq_record.features.sort(key=lambda f: f.location.start)
            if seq_status == "complete":
                logger.info("Query genome is qualified as complete. Will try to annotate UTR features.")
                _set_five_prime_utr(seq_record.features)
                _set_three_prime_utr(seq_record.features, len(seq_record))


    def check_completeness(self):
        models = get_vadr_models(self.vadr_dir)
        num_seqs = len(self.seq_dict.values())
        query_length = 0
        query_non_N_length = 0
        has_gap = False
        for seq_record in self.seq_dict.values():
            query_length += len(seq_record)
            query_N_length = str(seq_record.seq).upper().count("N")
            query_non_N_length += len(seq_record) - query_N_length
            for feature in seq_record.features:
                if feature.type == "assembly_gap":
                    has_gap = True

        self.report["query_length"] = query_length
        self.report["query_N_length"] = query_N_length
        self.report["query_non_N_length"] = query_non_N_length
        logger.debug(f"Checking completeness: num_seqs={num_seqs}, has_gap={has_gap}")
        if len(models) == 1:
            num_total_cds = 0
            num_cds_intact = 0
            num_cds_partial = 0
            num_cds_multi = 0
            num_cds_missing = 0

            model_features = self.models[models[0]]
            reference_length = model_features[0].attrs.get("length")
            if reference_length:
                reference_length = int(reference_length)
                self.report["reference_length"] = reference_length

            for model_feature in model_features:
                if model_feature.type != "CDS":
                    continue
                num_total_cds += 1
                if len(model_feature.hits) == 0:
                    num_cds_missing += 1
                elif len(model_feature.hits) > 1:
                    num_cds_multi += 1
                else:
                    if model_feature.hits[0] == "intact":
                        num_cds_intact += 1
                    elif model_feature.hits[0] == "partial":
                        num_cds_partial += 1
            cds_completeness = f"{num_cds_intact} / {num_cds_partial} / {num_cds_multi} / {num_cds_missing} [intact/partial/multi/missing]"
            cds_completeness_percentage = 100 * num_cds_intact / num_total_cds
            gap_ratio = (reference_length - query_non_N_length) / reference_length * 100
            if num_seqs == 1 and has_gap == False and int(cds_completeness_percentage) == 100 and gap_ratio < 10:
                seq_status = "complete"
            elif num_seqs == 1 and cds_completeness_percentage > 80 and gap_ratio < 20:
                seq_status = "nearly complete"
                INFO_NEARLY_COMPLETE_GENOME.add(f"CDS completeness: {cds_completeness_percentage:.2f}%")
                INFO_NEARLY_COMPLETE_GENOME.add(f"Gap ratio: {gap_ratio:.2f}%")
                self.warnings.append(INFO_NEARLY_COMPLETE_GENOME)
            else:
                seq_status = "draft"
                INFO_DRAFT_GENOME.add(f"CDS completeness: {cds_completeness_percentage:.2f}%")
                INFO_DRAFT_GENOME.add(f"Gap ratio: {gap_ratio:.2f}%")
                self.warnings.append(INFO_DRAFT_GENOME)
        else:
            cds_completeness = "na"
            cds_completeness_percentage = "na"
            seq_status = "draft"
            self.warnings.append(INFO_DRAFT_GENOME)

        self.report.update({"number_of_sequence": num_seqs,
                "has_assembly_gap": has_gap,
                "cds_completeness": cds_completeness, 
                "cds_completeness_percentage": cds_completeness_percentage,
                "seq_status": seq_status
            })

    def modify_incomplete_cds(self):
        def get_n_ratio(seq, feature):
            sequence = str(feature.extract(seq)).upper()
            n_count = sequence.count("N")
            total_len = len(sequence)
            return n_count / total_len 

        max_n_ratio = 0.5
        for seq_record in self.seq_dict.values():
            for feature in seq_record.features:
                if feature.type != "CDS":
                    continue
                incomplete = False
                note = []
                product = feature.qualifiers.get("product", [""])[0]
                translation = feature.qualifiers.get("translation", [""])[0]
                if isinstance(feature.location.start, BeforePosition):
                    INCOMPLETE_CDS_WARNING.add("no start codon:" + product)
                    note.append("no start codon")
                    incomplete = True
                if isinstance(feature.location.end, AfterPosition):
                    INCOMPLETE_CDS_WARNING.add("no stop codon:" + product)
                    note.append("no stop codon")
                    incomplete = True
                # if "XXXX" in translation:
                    # INCOMPLETE_CDS_WARNING.add("too many ambiguous AA:" + product)
                    # incomplete = True
                    # note.append("too many ambiguous AA")
                n_ratio = get_n_ratio(seq_record.seq, feature)
                if n_ratio > max_n_ratio:
                    feature.type = "misc_feature"
                    del feature.qualifiers["codon_start"]
                    del feature.qualifiers["transl_table"]
                    del feature.qualifiers["translation"]
                    del feature.qualifiers["product"]
                    if "ribosomal_slippage" in feature.qualifiers:
                        del feature.qualifiers["ribosomal_slippage"]
                    INCOMPLETE_CDS_WARNING.add("too many ambiguous AA:" + product)
                    incomplete = True
                if incomplete:
                    # The part below is currently disabled, but might be revived in the future
                    # feature.type = "misc_feature"
                    # del feature.qualifiers["codon_start"]
                    # del feature.qualifiers["transl_table"]
                    # del feature.qualifiers["translation"]
                    # del feature.qualifiers["product"]
                    # if "ribosomal_slippage" in feature.qualifiers:
                    #     del feature.qualifiers["ribosomal_slippage"]
                    if product:
                        feature.qualifiers.setdefault("note", []).append(f"Incomplete CDS similar to {product}")
                    if note:
                        # note = f"incomplete CDS: " + "; ".join(note)
                        # note = f"incomplete CDS: product={product} [" + "; ".join(note) + "]"
                        feature.qualifiers.setdefault("note", []).append("; ".join(note))
        if len(INCOMPLETE_CDS_WARNING.targets) > 0:
            self.warnings.append(INCOMPLETE_CDS_WARNING)

    def change_mat_peptide_to_misc(self):
        def _feature_is_in_CDS(feature, CDS_locations):
            for cds_loc in CDS_locations:
                if cds_loc.start <= feature.location.start and feature.location.end <= cds_loc.end:
                    return True
            return False

        CDS_locations = []
        for seq_record in self.seq_dict.values():
            for feature in seq_record.features:
                if feature.type == "CDS":
                    CDS_locations.append(feature.location)
        for seq_record in self.seq_dict.values():
            features = []
            for feature in seq_record.features:
                if feature.type == "mat_peptide":
                    if not _feature_is_in_CDS(feature, CDS_locations):
                        product = feature.qualifiers.get("product", ["mat_peptide"])[0] 
                        logger.warning(f"mat_peptide is changed to misc_feature because its parent CDS is not properly annotated. [{product}]")
                        feature.type = "misc_feature"
                        feature.qualifiers.setdefault("note", []).append(f"mat_peptide:{product}")
                        del feature.qualifiers["product"]
                   

    def to_gbk(self, output_file, with_feature_id=True):
        if with_feature_id:
            for r in self.seq_dict.values():
                for f in r.features:
                    f.qualifiers["ID"] = [f.id]
        SeqIO.write(self.seq_dict.values(), output_file, "genbank")

    def make_report(self, output_file, format="json"):

        def _check_alert(vadr_dir):
            vadr_prefix = os.path.basename(self.vadr_dir)
            alt_file = os.path.join(vadr_dir, vadr_prefix + ".vadr.alt")
            alerts = []
            for line in open(alt_file):
                cols = line.strip().split(maxsplit=9)
                if not cols[0].startswith("#"):
                    fail = cols[7]
                    alert_desc = cols[8]
                    alert_detail = cols[9]
                    target = f"{cols[1]}:{cols[3]}:{cols[4]}"
                    VADR_WARNING = create_VADR_warning(fail, alert_desc, alert_detail)
                    VADR_WARNING.add(target)
                    alerts.append(VADR_WARNING)
            return alerts

        # warnings = []

        # check number of models
        models = get_vadr_models(self.vadr_dir)
        if len(models) > 1:
            logger.warning("More than 1 reference model was used in VADR. Some of the validation steps of DFAST_VRL may not work properly. Please check the results and logs carefully.")
            MORE_THAN_ONE_MODEL_USED.targets = models
            self.warnings.append(MORE_THAN_ONE_MODEL_USED)

        # check VADR alert
        alerts = _check_alert(self.vadr_dir)
        self.warnings += alerts

        intact, duplicated, fragmented, partial, missing = [], [], [], [], []

        for model_name in models:
            for model_feature in self.models[model_name]:
                if model_feature.type is None or model_feature.type in ["stem_loop", "mat_peptide", "gene"]:
                    continue  # Don't check gene feature
                num_intact = model_feature.hits.count("intact")
                num_partial = model_feature.hits.count("partial")
                num_total = len(model_feature.hits)
                if num_total == 1:
                    if num_intact == 1:
                        intact.append(repr(model_feature))
                    if num_partial == 1:
                        partial.append(repr(model_feature))
                elif num_total > 1:
                    if num_intact >= 1:
                        duplicated.append(repr(model_feature))
                    if num_partial == num_total:
                        fragmented.append(repr(model_feature))
                elif num_total == 0:
                    if len(models) == 1:
                        missing.append(repr(model_feature))
        if duplicated:
            DUPLICATED_FEATURES.targets = duplicated
            self.warnings.append(DUPLICATED_FEATURES)
        if fragmented:
            FRAGMENTED_FEATURES.targets = fragmented
            self.warnings.append(FRAGMENTED_FEATURES)
        if partial:
            PARTIAL_FEATURES.targets = partial
            self.warnings.append(PARTIAL_FEATURES)
        if missing:
            MISSING_FEATURES.targets = missing
            self.warnings.append(MISSING_FEATURES)
        self.report.update({
            "intact_cds": intact,
            "duplicated_cds": duplicated,
            "fragmented_cds": fragmented,
            "partial_cds": partial,
            "missing_cds": missing,
            # "warnings": [warning.to_tuple() for warning in self.warnings]
        })
        
        report_json_string = json.dumps(self.report, indent=4)
        logger.info(f"vadr2gbk report:\n{report_json_string}\n")
        with open(output_file, "w") as f:
            f.write(report_json_string)
        # check duplication or fragmentation
        # check intact


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


def get_vadr_models(vadr_dir):
    vadr_prefix = os.path.basename(vadr_dir)
    mdl_file = os.path.join(vadr_dir, vadr_prefix + ".vadr.mdl")
    models = []
    for line in open(mdl_file):
        cols = line.strip().split()
        if cols[0].isnumeric():
            models.append(cols[1])
    return models

def convert_vadr_to_gbk(vadr_result_fasta, vadr_work_dir, output_gbk, isolate):
    v2d = VADR2DDBJ(vadr_result_fasta, vadr_work_dir, isolate=isolate)
    v2d.to_gbk(output_gbk)
    report_file = os.path.join(vadr_work_dir, "vadr.report.json")
    v2d.make_report(report_file)
    return output_gbk, {"vadr": v2d.report}, v2d.warnings




if __name__ == '__main__':
    
    ftr_file = "OUT/vadr/vadr.vadr.ftr"
    # input_fasta = "../meta_vry_result_rc.fa"
    input_fasta = "../keio/hcov-19_japan_donner476_2021.fasta"
    output_gbk = "OUT/test.vadr.gbk"
    vadr_dir = "OUT/vadr"
    report_file = "OUT/test.report.json"
    # ftr_file = sys.argv[1]
    # √ = sys.argv[2]
    # output_gbk = sys.argv[3]



    v2d = VADR2DDBJ(input_fasta, vadr_dir)
    v2d.to_gbk(output_gbk)
    v2d.make_report(report_file)
    # print(vsd.warnings)
"""
python /Users/tanizawa/dfast_web/app/scripts/genbank2mss.py
python genbank2mss.py <genbankFileName> <metadataFileName> <outDir> <filePrefix>
"""

