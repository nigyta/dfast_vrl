import os
import logging

from dfv.fix_mss_for_ful import parse_ftr_file, read_mss_file

logger = logging.getLogger(__name__)


def fix_hmpv_mss(work_dir, mss_file_prefix):
    """
    Fix MSS file for hMPV: set mol_type to "viral cRNA" (the sequence is in
    antigenome sense, as for RSV), and log the subgroup (A1/A2a/A2b1/A2b2/B1/B2)
    of the VADR model each sequence was assigned to.
    """
    mss_file = os.path.join(work_dir, f"{mss_file_prefix}.annt.tsv")
    features = list(read_mss_file(mss_file))

    seq_model = parse_ftr_file(work_dir) or {}  # {seq_id: model accession}
    subgroups = _parse_mdl_subgroups(work_dir)  # {model accession: subgroup}

    for feature in features:
        if feature[0][1] != "source":
            continue
        for qualifier in feature:
            if qualifier[3] == "mol_type":
                qualifier[4] = "viral cRNA"
        seq_id = feature[0][0]
        model = seq_model.get(seq_id)
        logger.info(f"hMPV sequence {seq_id}: model={model}, subgroup={subgroups.get(model)}")

    with open(mss_file, "w") as f:
        f.write("\n".join("\t".join(row) for feature in features for row in feature) + "\n")
    logger.info(f"Fixed MSS file for hMPV: {mss_file}")


def _parse_mdl_subgroups(work_dir):
    # rows: idx model group subgroup seqs pass fail
    subgroups = {}
    for line in open(os.path.join(work_dir, "vadr", "vadr.vadr.mdl")):
        cols = line.split()
        if len(cols) >= 7 and cols[0].isdigit():
            subgroups[cols[1]] = cols[3]
    return subgroups
