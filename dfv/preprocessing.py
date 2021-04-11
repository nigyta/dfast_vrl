import os
import sys
import io
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

from logging import getLogger
logger = getLogger(__name__)


class Preprocessing:
    def __init__(self, query_fasta, subject_fasta, work_dir):
        self.query_fasta = query_fasta
        self.subject_fasta = subject_fasta
        self.work_dir = work_dir
        self.round = 0
        self.hit_history = []
        self.converged = False
        if not os.path.exists(self.work_dir):
            os.makedirs(self.work_dir, exist_ok=True)

    def get_current_files(self):
        current_query = os.path.join(self.work_dir, f"query_{self.round}.fa")
        current_subject = os.path.join(self.work_dir, f"subject_{self.round}.fa")
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

        current_query, current_subject = self.get_current_files()
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
        current_query, current_subject = self.get_current_files()
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
            print(f"---- Round {self.round}")
            print(query_id, hit_id)
            print(top_hit)
            print(f"Identity={top_hit.identities}/{top_hit.align_length}={100*top_hit.identities/top_hit.align_length:.2f}%")
            # print("qstart,qend,sstart,send", top_hit.query_start, top_hit.query_end, top_hit.sbjct_start, top_hit.strand)
            assert top_hit.query_start < top_hit.query_end
            self.hit_history.append(hsps[0])
            self.round += 1
        else:
            print(f"Converged at round-{self.round}!")
            self.converged = True


    def write_output(self, output_fasta, scaffolding=False):
        def _sort_func(hit):
            query_id, hit_id, hsp = hit
            return hsp.sbjct_start if hsp.strand == ("Plus", "Plus") else hsp.sbjct_end

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
                length = abs(next_hsp.sbjct_start - hsp.sbjct_end) - 1
                sequences.append("N" * length)
        with open(output_fasta, "w") as f:
            logger.info(f"Writing preprocessed FASTA file to {output_fasta}")
            print(f"Writing preprocessed FASTA file to {output_fasta}")
            if scaffolding:
                seq = "".join(sequences)
                f.write(f">sequence\n{seq}")
            else:
                for idx, seq in enumerate(sequences, 1):
                    f.write(f">sequence{idx}\n{seq}\n")

def preprocess_contigs(input_fasta, work_dir, output_fasta=None, reference_fasta=None, scaffolding=False):
    if reference_fasta is None:
        reference_fasta = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../refs/NC_045512.2.fasta")
    if output_fasta is None:
        output_fasta = os.path.join(work_dir, "preprocessed.fa")

    pp = Preprocessing(input_fasta, reference_fasta, work_dir)
    while not pp.converged:
        pp.make_query_and_subject()
        pp.run_blast()
        pp.parse_blast_result()
    pp.write_output(output_fasta, scaffolding=scaffolding)
    return output_fasta
# record attrs = todo
# alignment attrs = ['accession', 'hit_def', 'hit_id', 'hsps', 'length', 'title']
# hsp attrs = ['align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives', 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']

if __name__ == '__main__':
    input_fasta = sys.argv[1]
    preprocess_contigs(input_fasta)