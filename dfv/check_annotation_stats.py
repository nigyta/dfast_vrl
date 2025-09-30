import os
import glob
import logging
from dataclasses import dataclass, field
import subprocess
from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition

from .reference_models import parse_minfo_file

logger = logging.getLogger(__name__)


def get_minfo_file(model):
    minfo_path_with_envvar = model.minfo_file
    cmd = f"ls -1 {minfo_path_with_envvar}"
    minfo_path = subprocess.run(cmd, shell=True, encoding="UTF-8", capture_output=True, text=True).stdout.strip()
    if not os.path.exists(minfo_path):
        logger.error("Minfo file not found. %s", model.minfo_file)
        exit(1)
    return minfo_path

def get_model_accession(mdl_file):
    model_accession = None
    for line in open(mdl_file):
        if line.startswith("1"):
            model_accession = line.split()[1]
        elif line.startswith("#") or line.startswith("-"):
            pass  # ok
        else:
            model_accession = "-"
            # logger.error("Unexpected .mdl file. %s", line)
            # exit(1)
    assert model_accession
    return model_accession

def get_model_info(out_dir, minfo_file):
    prefix = os.path.basename(out_dir)
    mdl_files = glob.glob(os.path.join(out_dir, prefix + ".vadr.mdl"))
    if len(mdl_files) == 0:
        logger.error("No .vadr.mdl file found. VADR may have failed. %s", mdl_files)
        exit(1)
    elif len(mdl_files) > 1:
        logger.error("Multiple .vadr.mdl file found. VADR may have failed. Not expected. Aborted. %s", mdl_files)
        exit(1)
    mdl_file = mdl_files[0]                          
    model_accession = get_model_accession(mdl_file)
    # print(model_accession)
    dict_model_info = parse_minfo_file(minfo_file)
    model_info = dict_model_info.get(model_accession)
    if model_info:
        # print(model_info)
        model_length = model_info[0].attrs.get("length")
        # print(model_length)
        model_num_cds = len([x for x in model_info if x.type == "CDS"])
        # print(num_cds)
    else:
        # not available for models with multi seqments
        model_length = None
        model_num_cds = None
    return model_accession, model_length, model_num_cds

def get_query_info(out_dir):
    gbk_file = os.path.join(out_dir, "annotation.gbk")
    if not os.path.exists(gbk_file):
        logger.error("GBK file %s not found. VADR may have failed.", gbk_file)
        exit(1)
    # print(gbk_file)
    R = list(SeqIO.parse(gbk_file, "genbank"))
    num_sequence = len(R)
    query_total_length = sum([len(r.seq) for r in R])
    gap_length = sum([str(r.seq).upper().count("N") for r in R])
    num_cds_intact, num_cds_partial = 0, 0
    for r in R:
        for f in r.features:
            if f.type != "CDS":
                continue
            if isinstance(f.location.start, BeforePosition) or isinstance(f.location.end, AfterPosition):
                num_cds_partial += 1
            else:
                num_cds_intact += 1
    return num_sequence, query_total_length, gap_length, num_cds_intact, num_cds_partial

def check_annotation_stats(out_dir, vadr_dir, model):
    number_of_sequence, query_total_length, gap_length, num_cds_intact, num_cds_partial = get_query_info(out_dir)
    minfo_file = get_minfo_file(model)
    model_accession, model_length, model_num_cds = get_model_info(vadr_dir, minfo_file)
    if model_length:
        model_length = int(model_length)
        query_coverage = query_total_length / model_length
        query_coverage = f"{query_coverage:.2%}"
        if query_coverage > 0.98 and gap_length == 0:
            status = "complete"
        elif query_coverage > 0.9:
            status = "nearly complete"
        else:
            status = "draft"

    else:
        query_coverage = "n.a."
        status = "n.a."
    if model_num_cds:
        num_cds_missing = model_num_cds - num_cds_intact - num_cds_partial
    else:
        num_cds_missing = "n.a."
    # status:
    # complete: coverage > 98%, qap_length = 0
    # nearly complete: coverage > 90%
    # draft: other
    # if query_coverage > 0.98 and gap_length == 0:
    #     status = "complete"
    # elif query_coverage > 0.9:
    #     status = "nearly complete"
    # else:
    #     status = "draft"

    ret = {
        "status": status,
        "number_of_sequence": number_of_sequence,
        "total_length": query_total_length,
        "model_length": model_length,
        "query_coverage": query_coverage,
        "qap_length": gap_length,
        # "cds_completeness": f"{num_cds_intact} / {num_cds_partial} / {num_cds_missing} / {model_num_cds} [intact/partial/missing/total]"
        "cds_completeness": f"{num_cds_intact} / {num_cds_partial} / {model_num_cds} [intact/partial/expected]"
    }
    return ret

if __name__ == '__main__':
    minfo_file = "/opt/vadr/vadr-models-flavi/dengue.minfo"
    out_dir = "/data/out_dengu_OM281574_mss2vadr"
    # print(get_query_info(out_dir))
    from .vadr2mss_config import models
    model = models["Dengue"]
    ret = check_annotation_stats(out_dir, model)
    print(ret)
""" mdl file
#                                   num   num   num
#idx  model      group   subgroup  seqs  pass  fail
#---  ---------  ------  --------  ----  ----  ----
1     NC_001477  Dengue  1            1     1     0
#---  ---------  ------  --------  ----  ----  ----
-     *all*      -       -            1     1     0
-     *none*     -       -            0     0     0
#---  ---------  ------  --------  ----  ----  ----
"""