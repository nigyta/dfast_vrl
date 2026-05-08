import os
import logging

from dfv.fix_mss_for_ful import parse_ftr_file, read_mss_file

logger = logging.getLogger(__name__)

ORGANISM_BY_SUBTYPE = {
    "A": "Human respiratory syncytial virus A",
    "B": "Human respiratory syncytial virus B",
}


def fix_rsv_mss(work_dir, mss_file_prefix, out_mss_file=None):
    """
    Fix MSS file for RSV (overwrite organism with subtype-specific name and mol_type).
    """
    mss_file = os.path.join(work_dir, f"{mss_file_prefix}.annt.tsv")
    if not os.path.exists(mss_file):
        logger.error(f"MSS file not found: {mss_file}")
        return

    dict_model = parse_mdl_file_rsv(work_dir)
    logger.debug(f"Parsed model info: {dict_model}")

    dict_seq_model = parse_ftr_file(work_dir)
    logger.debug(f"Parsed sequence to model mapping: {dict_seq_model}")

    out_buffer = []
    for feature in read_mss_file(mss_file):
        feature_type = feature[0][1]
        if feature_type == "source":
            seq_id = feature[0][0]
            if seq_id in dict_seq_model:
                model = dict_seq_model.get(seq_id)
                if model is None:
                    logger.error(f"Model not found for sequence: {seq_id}")
                    exit(1)
                subtype = dict_model.get(model)
                if subtype is None:
                    logger.error(f"Model info not found for model: {model}")
                    exit(1)

                organism = ORGANISM_BY_SUBTYPE.get(subtype)
                if organism is None:
                    logger.warning(f"Unknown RSV subtype '{subtype}' for sequence {seq_id}; organism will not be overwritten.")

                for qualifier in feature:
                    q_key = qualifier[3]
                    if q_key == "organism" and organism:
                        qualifier[4] = organism
                    if q_key == "mol_type":
                        qualifier[4] = "viral cRNA"
                logger.info(f"Updated source for sequence {seq_id}: subtype={subtype}, organism={organism}")
        out_buffer.append(feature)

    if out_mss_file is None:
        out_mss_file = mss_file
        logger.info(f"Overwriting original MSS file: {mss_file}")
    with open(out_mss_file, "w") as f:
        out_str = "\n".join(["\t".join(row) for feature in out_buffer for row in feature])
        f.write(out_str + "\n")
    logger.info(f"Fixed MSS file for RSV: {out_mss_file}")


def parse_mdl_file_rsv(work_dir):
    """
    Parse MDL file to map model id to RSV subtype.
    The VADR RSV model MDL has group="RSV" and subgroup="A" or "B".
    Returns: {model_id: subtype} where subtype is "A" or "B".

    Example MDL content:
    #idx  model     group  subgroup  seqs  pass  fail
    1     MZ516105  RSV    B            1     1     0
    """
    mdl_file = os.path.join(work_dir, "vadr", "vadr.vadr.mdl")
    if not os.path.exists(mdl_file):
        logger.error(f"MDL file not found: {mdl_file}")
        return None

    dict_model = {}
    with open(mdl_file) as f:
        for line in f:
            if line.startswith("#") or line.startswith("-"):
                continue
            cols = line.strip().split()
            if len(cols) < 6:
                continue
            model = cols[1]
            subgroup = cols[3]
            if subgroup in ("A", "B"):
                subtype = subgroup
            else:
                logger.warning(f"Could not determine RSV subtype from subgroup: {subgroup}")
                subtype = "unknown_subtype"
            dict_model[model] = subtype
    return dict_model


if __name__ == "__main__":
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())
    fix_rsv_mss("OUT_rsv", "DDBJ")
    print("done")
