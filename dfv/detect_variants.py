import os
import sys
import json
import subprocess
from dataclasses import dataclass
import re

from Bio import SeqIO
if __name__ == '__main__':
    from dfv_warnings import INFO_QUERY_MODIFIED, INFO_SCAFFOLDING_ENABLED, TRIM_TERMINAL_N
else:
    from .dfv_warnings import INFO_QUERY_MODIFIED, INFO_SCAFFOLDING_ENABLED, TRIM_TERMINAL_N

app_root = os.path.join(os.path.dirname(os.path.dirname(__file__)))
snpeff_conf_file = os.path.join(app_root, "refs/snpeff/snpEff.config")

import logging

logger = logging.getLogger(__name__)

pat_eff = re.compile(r"EFF=(.+?)(($)|(;)|(\t))")

@dataclass
class Var:
    idx: int
    r_pos: int
    q_pos: int
    ref: str
    alt: str

    def to_vcf(self):
        return "\t".join(["NC_045512.2", str(self.r_pos), ".", self.ref.replace("-", ""), self.alt.replace("-", ""), ".", ".", "."])


def msa2vcf(in_msa, out_vcf):
    out = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    R = [str(r.seq).upper() for r in SeqIO.parse(in_msa, "fasta")]
    query = R[0]
    reference = R[1]
    assert len(query) == len(reference)
    r_pos, q_pos = 0, 0
    variants = []
    prev_var = None
    for idx, (ref, alt) in enumerate(zip(reference, query)):
        if ref != "-":
            r_pos += 1
        if alt != "-":
            q_pos += 1
        # print(idx, r_pos, q_pos, ref, alt)
        if ref != alt:
            if alt == "N":
                continue
            var = Var(idx, r_pos, q_pos, ref, alt)
            # print(var)
            if prev_var and (prev_var.idx + len(prev_var.ref) == var.idx):
                prev_var.ref += var.ref
                prev_var.alt += var.alt
            else:
                variants.append(var)
                prev_var = var
    for var in variants:    
        if var.ref.startswith("-"):
            var.idx -= 1
            var.q_pos -= 1
            prev_pos_r = reference[var.idx]
            prev_pos_q = query[var.idx]
            assert prev_pos_r == prev_pos_q
            var.ref = prev_pos_r + var.ref
            var.alt = prev_pos_q + var.alt
        elif var.alt.startswith("-"):
            var.idx -= 1
            var.r_pos -= 1
            prev_pos_r = reference[var.idx]
            prev_pos_q = query[var.idx]
            assert prev_pos_r == prev_pos_q
            var.ref = prev_pos_r + var.ref
            var.alt = prev_pos_q + var.alt
        # print(var)
        out += var.to_vcf() + "\n"
    with open(out_vcf, "w") as f:
        f.write(out)

def parse_snpeff_vcf(in_vcf, out_json):
    def parse_anno(anno):
        effect_type = anno.split("(")[0]
        anno = anno.replace(effect_type, "").strip("()")
        values = anno.split("|")
        aa_change = values[3]
        gene_name = values[5]
        return gene_name, aa_change, effect_type

    ret = []
    for line in open(in_vcf):
        if line.startswith("#"):
            continue
        else:
            cols = line.strip("\n").split("\t")
            pos, ref, alt, info = cols[1], cols[3], cols[4], cols[7]
            nt_change = f"{ref}{pos}{alt}"
            m = pat_eff.search(info)
            annos = []
            if m:
                annos = m.group(1)
                annos = annos.split(",")
                # print(annos)
            for anno in annos:
                gene_name, aa_change, effect_type = parse_anno(anno)
                if effect_type == "SYNONYMOUS_CODING":
                    aa_change = ""
                ret.append({"gene": gene_name, "aa_change": aa_change, "nt_change": nt_change, "type": effect_type})
    return ret

def associate_variants_to_gene(variants):
    D = {}
    for dict_var in variants:  # example of dict_var = {"gene": gene_name, "aa_change": aa_change, "nt_change": nt_change, "type": effect_type}
        gene = dict_var["gene"]
        aa_change = dict_var["aa_change"]
        effect_type = dict_var["type"]
        nt_change = dict_var["nt_change"]
        if gene and aa_change:
            D.setdefault(gene, []).append((aa_change, nt_change, effect_type))
    return D

def detect_variants(input_fasta, work_dir):
    logger.info(f"Detecting variants")
    if not os.path.exists(work_dir):
        os.makedirs(work_dir, exist_ok=True)
    out_msa = os.path.join(work_dir, "msa.fasta")
    out_vcf = os.path.join(work_dir, "var.vcf")
    out_snpeff = os.path.join(work_dir, "var.anno.vcf")
    out_variant_json = os.path.join(work_dir, "var.anno.json")

    mafft_cmd = f"mafft {input_fasta} > {out_msa}"
    logger.debug(f"MAFFT command: {mafft_cmd}")
    p = subprocess.run(mafft_cmd, shell=True, encoding="UTF-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # if p.stderr:
        # logger.info(f"{'='*30} MAFFT stderr {'='*30}\n{p.stderr}\n{'='*80}\n")
    my_out_vcf = os.path.join(work_dir, "my.var.vcf")
    msa2vcf(out_msa, out_vcf)
    snpeff_cmd = f"snpeff -c {snpeff_conf_file} -noStats -no-downstream -no-upstream -no-utr -classic -formatEff nigvrl {out_vcf} > {out_snpeff}"
    logger.debug(f"SNPeff command: {snpeff_cmd}")
    p = subprocess.run(snpeff_cmd, shell=True, encoding="UTF-8", stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # if p.stderr:
        # logger.info(f"{'='*30} SNPeff stderr {'='*30}\n{p.stderr}\n{'='*80}\n")
    variants = parse_snpeff_vcf(out_snpeff, "out_json")
    variants_by_gene = associate_variants_to_gene(variants)
    with open(out_variant_json, "w") as f:
        json.dump(variants, f, indent=4)

    return variants_by_gene

if __name__ == "__main__":
    detect_variants("ppout/pp.msa.fa", "ppout")
