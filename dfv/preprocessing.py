import os
import sys
import io
import json
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
if __name__ == '__main__':
    from dfv_warnings import INFO_QUERY_MODIFIED, INFO_SCAFFOLDING_ENABLED, TRIM_TERMINAL_N, SEQUENCE_RENAMED, TRIM_SLASH
else:
    from .dfv_warnings import INFO_QUERY_MODIFIED, INFO_SCAFFOLDING_ENABLED, TRIM_TERMINAL_N, SEQUENCE_RENAMED, TRIM_SLASH


import logging
logger = logging.getLogger(__name__)


class Preprocessing:
    def __init__(self, query_fasta, subject_fasta, work_dir, modify_query_fasta=False):
        self.query_fasta = query_fasta
        self.subject_fasta = subject_fasta
        self.work_dir = work_dir
        self.round = 1
        self.hit_history = []
        self.report = {}
        self.warnings = []
        self.converged = False
        self.modify_query_fasta = modify_query_fasta
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir, exist_ok=True)
        logger.info("Quality check and preprocessing started.")
        logger.info(f"Query: {os.path.basename(self.query_fasta)}")
        logger.info(f"Reference: {os.path.basename(self.subject_fasta)}")

    def trim_slash_from_query(self):
        out_fasta = os.path.join(self.work_dir, "query_slash_trimmed.fa")
        R = list(SeqIO.parse(self.query_fasta, "fasta"))
        len_trimmed = 0
        len_original = 0
        str_output = ""
        flag_slash_trimmed = False
        seq_id_original = []
        for r in R:
            seq_original = str(r.seq).upper()
            seq_trimmed = seq_original.strip("/")
            len_original += len(seq_original)
            len_trimmed += len(seq_trimmed)
            str_output += f">{r.id}\n{seq_trimmed}\n"
        if len_trimmed != len_original:
            logger.warning(f"Trimmed trailing '//' in the query FASTA. Length: {len_original} --> {len_trimmed}")
            TRIM_SLASH.add(f"Length: {len_original} --> {len_trimmed}")
            self.warnings.append(TRIM_SLASH)
            flag_slash_trimmed = True
        if flag_slash_trimmed:
            with open(out_fasta, "w") as f:
                f.write(str_output)
            self.query_fasta = out_fasta 

    def trim_N_from_query(self):
        out_fasta = os.path.join(self.work_dir, "query_N_trimmed.fa")
        R = list(SeqIO.parse(self.query_fasta, "fasta"))
        len_trimmed = 0
        len_original = 0
        str_output = ""
        flag_seq_trimmed, flag_seq_renamed = False, False
        seq_id_original = []
        for r in R:
            seq_id = r.id
            seq_original = str(r.seq).upper()
            seq_trimmed = seq_original.strip("N")
            len_original += len(seq_original)
            len_trimmed += len(seq_trimmed)
            if len(seq_id) > 50:
                seq_id_original.append(seq_id)
                seq_id = seq_id[:50]
                flag_seq_renamed = True
            str_output += f">{seq_id}\n{seq_trimmed}\n"
        if len_trimmed != len_original:
            logger.warning(f"Trimmed leading/trailing Ns in the query FASTA. Length: {len_original} --> {len_trimmed}")
            TRIM_TERMINAL_N.add(f"Length: {len_original} --> {len_trimmed}")
            self.warnings.append(TRIM_TERMINAL_N)
            flag_seq_trimmed = True
        if flag_seq_renamed:
            logger.warning(f"Sequence name is too long. Truncated to 50 letters.")
            SEQUENCE_RENAMED.add(", ".join(seq_id_original))            
            self.warnings.append(SEQUENCE_RENAMED)
        if flag_seq_trimmed or flag_seq_renamed:
            with open(out_fasta, "w") as f:
                f.write(str_output)
            self.query_fasta = out_fasta


    def get_fasta_files(self, round=None):
        if round is None:
            round = self.round
        current_query = os.path.join(self.work_dir, f"query_{round}.fa")
        current_subject = os.path.join(self.work_dir, f"subject_{round}.fa")
        return current_query, current_subject

    def get_blast_out(self):
        return os.path.join(self.work_dir, f"blast_{self.round}.xml")

    def make_query_and_subject(self):
        def _mask_sequence(seq_dict, hit_history, target):
            for query_id, hit_id, hsp in hit_history:
                if target == "query":
                    start, end = min(hsp.query_start, hsp.query_end) - 1, max(hsp.query_start, hsp.query_end)  # HSP is 1-based
                    seq = seq_dict[query_id]
                    seq = seq[:start] + "N" * (end - start) + seq[end:]
                    seq_dict[query_id] = seq
                else:
                    start, end = min(hsp.sbjct_start, hsp.sbjct_end) - 1, max(hsp.sbjct_start, hsp.sbjct_end)  # HSP is 1-based
                    seq = seq_dict[hit_id]
                    seq = seq[:start] + "N" * (end - start) + seq[end:]
                    seq_dict[hit_id] = seq

        current_query, current_subject = self.get_fasta_files()
        with open(current_query, "w") as f:
            seq_dict = {r.id: str(r.seq) for r in SeqIO.parse(self.query_fasta, "fasta")}
            _mask_sequence(seq_dict, self.hit_history, target="query")
            for seq_id, seq in seq_dict.items():
                f.write(f">{seq_id}\n{seq}\n")
        with open(current_subject, "w") as f:
            seq_dict = {r.id: str(r.seq) for r in SeqIO.parse(self.subject_fasta, "fasta")}
            _mask_sequence(seq_dict, self.hit_history, target="subject")
            for seq_id, seq in seq_dict.items():
                f.write(f">{seq_id}\n{seq}\n")

    def run_blast(self):
        current_query, current_subject = self.get_fasta_files()
        blast_out = self.get_blast_out()
        out, err = NcbiblastnCommandline(query=current_query, subject=current_subject, out=blast_out, evalue=1e-10, outfmt=5, task="blastn")()

    def parse_blast_result(self):
        blast_out = self.get_blast_out()
        hsps = []
        blast_records = NCBIXML.parse(open(blast_out))
        for r in blast_records:
            for aln in r.alignments:
                for hsp in aln.hsps:
                    hsps.append((r.query, aln.hit_id, hsp))
        hsps = sorted(hsps, key=lambda hit: hit[2].score, reverse=True)
        if hsps:
            query_id, hit_id, top_hit = hsps[0]
            logger.info(f"-----   Round {self.round}   -----")
            logger.info(f"Query: {query_id}, Subject: {hit_id}")
            logger.info(f"Identity={top_hit.identities}/{top_hit.align_length}={100*top_hit.identities/top_hit.align_length:.2f}%")
            # logger.debug(f"qstart, qend, sstart, send, strand : {top_hit.query_start}, {top_hit.query_end}, {top_hit.sbjct_start},  {top_hit.sbjct_end}, {top_hit.strand}")
            logger.info(f"\n{top_hit}\n")
            assert top_hit.query_start < top_hit.query_end
            self.hit_history.append(hsps[0])
            self.round += 1
        else:
            logger.info(f"Preprocessing completed at Round {self.round}.")
            self.converged = True


    def write_output(self, output_fasta, scaffolding=False):
        def _sort_func(hit):
            query_id, hit_id, hsp = hit
            return hsp.sbjct_start if hsp.strand == ("Plus", "Plus") else hsp.sbjct_end

        if scaffolding and len(self.hit_history) > 1:
            logger.warning("Scaffolding is enabled. Contigs will be concatenated with runs of Ns of estimated length.")
            self.warnings.append(INFO_SCAFFOLDING_ENABLED)
        hsps_sorted = sorted(self.hit_history, key=lambda hit: hit[2].sbjct_start)  # hit = tuple of (query_id, hit_id, hsp)
        sequences = []
        seq_dict = {r.id: r.seq for r in SeqIO.parse(self.query_fasta, "fasta")}
        for idx, (query_id, hit_id, hsp) in enumerate(hsps_sorted, 1):
            seq = seq_dict[query_id]
            if hsp.strand == ("Plus", "Plus"):
                sequences.append(str(seq[hsp.query_start - 1:hsp.query_end]).upper())
            else:
                sequences.append(str(seq[hsp.query_start - 1:hsp.query_end].reverse_complement().upper()))
            if idx < len(hsps_sorted) and scaffolding:
                _, __, next_hsp = hsps_sorted[idx]
                start = max(hsp.sbjct_end, hsp.sbjct_start)
                end = min(next_hsp.sbjct_start, next_hsp.sbjct_end) - 1
                length = end - start

                # length = abs(next_hsp.sbjct_start - hsp.sbjct_end) - 1
                logger.warning(f'Contig{idx} and contig{idx+1} has been concatenated with a gap of {length} nt length.')
                sequences.append("N" * length)
        with open(output_fasta, "w") as f:
            logger.info(f"Writing preprocessed FASTA file to {output_fasta}")
            if scaffolding:
                seq = "".join(sequences)
                f.write(f">sequence\n{seq}\n")
            else:
                for idx, seq in enumerate(sequences, 1):
                    f.write(f">sequence{idx}\n{seq}\n")

    def write_fasta_for_msa(self, output_fasta):
        def _sort_func(hit):
            query_id, hit_id, hsp = hit
            return hsp.sbjct_start if hsp.strand == ("Plus", "Plus") else hsp.sbjct_end

        hsps_sorted = sorted(self.hit_history, key=lambda hit: hit[2].sbjct_start)  # hit = tuple of (query_id, hit_id, hsp)
        sequences = []
        seq_dict = {r.id: r.seq for r in SeqIO.parse(self.query_fasta, "fasta")}
        seq_list_sbjct = [str(r.seq) for r in SeqIO.parse(self.subject_fasta, "fasta")]
        assert len(seq_list_sbjct) == 1
        sbcjt_seq = seq_list_sbjct[0]
        # add 5'-end
        end = min(hsps_sorted[0][2].sbjct_start, hsps_sorted[0][2].sbjct_end) - 1
        if end > 0:
            sequences.append(sbcjt_seq[:end])
        for idx, (query_id, hit_id, hsp) in enumerate(hsps_sorted, 1):
            # print(idx, query_id, hit_id)
            # print(hsp)
            seq = seq_dict[query_id]
            if hsp.strand == ("Plus", "Plus"):
                sequences.append(str(seq[hsp.query_start - 1:hsp.query_end]).upper())
            else:
                sequences.append(str(seq[hsp.query_start - 1:hsp.query_end].reverse_complement().upper()))
            if idx < len(hsps_sorted):
                _, __, next_hsp = hsps_sorted[idx]
                start = max(hsp.sbjct_end, hsp.sbjct_start)
                end = min(next_hsp.sbjct_start, next_hsp.sbjct_end) - 1
                sequences.append(str(sbcjt_seq[start:end]))
                # length = end - start
        # add 3'-end
        start = max(hsps_sorted[-1][2].sbjct_start, hsps_sorted[-1][2].sbjct_end)
        if start < len(sbcjt_seq):
            sequences.append(sbcjt_seq[start:])
        with open(output_fasta, "w") as f:
            logger.info(f"Writing FASTA file for MSA to {output_fasta}")
            seq = "".join(sequences)
            f.write(f">query\n{seq}\n")
            f.write(f">sbjct\n{sbcjt_seq}\n")


    def set_report(self):
        def _get_length(fasta_file):
            R = list(SeqIO.parse(fasta_file, "fasta"))
            length = sum([len(r) for r in R])
            num_seq = len(R)
            return length, num_seq

        def _get_aligned_length(hsp):
            query_aligned = hsp.query_end - hsp.query_start + 1
            sbjct_aligned = abs(hsp.sbjct_end - hsp.sbjct_start) + 1
            return query_aligned, sbjct_aligned

        current_query, current_subject = self.get_fasta_files()
        query_length, query_num = _get_length(current_query)        
        sbjct_length, subject_num = _get_length(current_subject)

        sum_aligned_query, sum_aligned_sbjct, sum_identities, sum_aligned = 0, 0, 0, 0
        for query_id, hit_id, hsp in self.hit_history:
            sum_aligned += hsp.align_length
            query_aligned, sbjct_aligned = _get_aligned_length(hsp)
            sum_aligned_query += query_aligned
            sum_aligned_sbjct += sbjct_aligned
            sum_identities += hsp.identities

        if sum_aligned_query != query_length or query_num != len(self.hit_history):
            INFO_QUERY_MODIFIED.add(f"Length: {query_length} --> {sum_aligned_query}")
            if self.modify_query_fasta:
                self.warnings.append(INFO_QUERY_MODIFIED)
        self.report = {
            "query_length": query_length,
            "query_num_sequence": query_num,
            "sbjct_length": sbjct_length,
            "query_aligned_length": sum_aligned_query,
            "sbjct_aligned_length": sum_aligned_sbjct,
            "matched_nucleotides": sum_identities,
            "matched_fragments": len(self.hit_history),
            "average_identity": sum_identities / sum_aligned * 100,
            "query_coverage": sum_aligned_query / query_length * 100,
            "sbjct_coverage": sum_aligned_sbjct / sbjct_length * 100,
        }

    def write_report(self, output_file=None, format="json"):
        if format == "json":
            if output_file is None:
                output_file = os.path.join(self.work_dir, "preprocessing_report.json")
            with open(output_file, "w") as f:
                logger.info(f"Writing quality control and preprocessing report file to {output_file}")
                json.dump({"preprocessing": self.report}, f, indent=4)
        else:
            raise NotImplementedError


def preprocess_contigs(input_fasta, work_dir, output_fasta=None, reference_fasta=None, modify_query_fasta=False, disable_scaffolding=False):
    if reference_fasta is None:
        reference_fasta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../refs/NC_045512.2.fasta")
    if output_fasta is None:
        output_fasta = os.path.join(work_dir, "preprocessed.fa")

    scaffolding = not disable_scaffolding

    pp = Preprocessing(input_fasta, reference_fasta, work_dir, modify_query_fasta)
    pp.trim_slash_from_query()
    pp.trim_N_from_query()
    while not pp.converged:
        pp.make_query_and_subject()
        pp.run_blast()
        pp.parse_blast_result()
    pp.set_report()
    pp.write_report()
    logger.info(json.dumps(pp.report, indent=4))
    if pp.modify_query_fasta:
        logger.warning("### Since 'modify_query_fasta' is enabled, original input FASTA will be modified based on the mapping result to the reference. ###")
        pp.write_output(output_fasta, scaffolding=scaffolding)
        result_fasta = output_fasta
    else:
        pp.write_output(output_fasta, scaffolding=False)  # when modify_query_fasta is enabled, scaffolding is always disabled.
        # logger.warning("### Since preprocessing is disabled, original input FASTA will be used for the downstream steps. ###")
        result_fasta = pp.query_fasta  # query fasta (or cleaned fasta) is used in the downloading process
    pp.write_fasta_for_msa(output_fasta=os.path.join(work_dir, "msa_input.fasta"))

    return result_fasta, {"preprocessing": pp.report}, pp.warnings

# record attrs = todo
# alignment attrs = ['accession', 'hit_def', 'hit_id', 'hsps', 'length', 'title']
# hsp attrs = ['align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    # logger.setLevel(logging.DEBUG)
    input_fasta = sys.argv[1]
    output_dir = sys.argv[2]
    output_fasta, report, warnings = preprocess_contigs(input_fasta, output_dir, modify_query_fasta=True, disable_scaffolding=False)
    print("output_fasta", output_fasta)
    print(report)