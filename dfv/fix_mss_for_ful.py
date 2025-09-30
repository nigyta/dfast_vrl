import os
import logging

logger = logging.getLogger(__name__)

def fix_flu_mss(work_dir, mss_file_prefix, out_mss_file=None):
    """
    Fix MSS file for influenza virus (adding segment information, etc.)
    """

    mss_file = os.path.join(work_dir, f"{mss_file_prefix}.annt.tsv")
    if not os.path.exists(mss_file):
        logger.error(f"MSS file not found: {mss_file}")
        return

    
    dict_model = parse_mdl_file(work_dir)
    logger.debug(f"Parsed model info: {dict_model}")

    dict_seq_model = parse_ftr_file(work_dir)
    logger.debug(f"Parsed sequence to model mapping: {dict_seq_model}")

    subgroup_h = [subgroup for flu_type, segment, subgroup in dict_model.values() if subgroup.startswith("H")]
    subgroup_n = [subgroup for flu_type, segment, subgroup in dict_model.values() if subgroup.startswith("N")]
    if len(subgroup_h) == 1 and len(subgroup_n) == 1:
        serotype = f"{subgroup_h[0]}{subgroup_n[0]}"
        logger.info(f"Determined serotype by VADR: {serotype}")
    elif len(subgroup_h) > 1 or len(subgroup_n) > 1:
        logger.warning(f"Multiple H or N subgroups found: H={subgroup_h}, N={subgroup_n}. Will not set serotype.")
        serotype = ""
    else:
        logger.warning(f"Could not determine serotype from subgroups: H={subgroup_h}, N={subgroup_n}. Will not set serotype.")
        serotype = ""


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
                model_info = dict_model.get(model)
                if model_info is None:
                    logger.error(f"Model info not found for model: {model}")
                    exit(1)
                flu_type, segment, subgroup = model_info

                serotype_qualifier = [q[4] for q in feature if q[3] == "serotype"]
                if serotype_qualifier:
                    serotype_qualifier = serotype_qualifier[0]
                else:
                    serotype_qualifier = ""

                if serotype_qualifier and serotype != serotype_qualifier:
                    logger.error(f"serotype in MSS file ({serotype_qualifier}) is different from determined serotype ({serotype})")
                    exit(1)
                if serotype != "" and serotype_qualifier == "":
                    serotype_qualifier = serotype
                    feature.append(["", "", "", "serotype", serotype_qualifier])
                    logger.info(f"Added serotype qualifier infered by VADR: {serotype_qualifier} to sequence {seq_id}")

                for qualifier in feature:
                    q_key = qualifier[3]
                    q_value = qualifier[4]
                    if q_key == "organism":
                        # print(f"serotype_qualifier: {serotype_qualifier}, serotype: {serotype}, flu_type: {flu_type}, organism: {q_value}")
                        if flu_type == "A":
                            if serotype_qualifier:
                                q_value = f"Influenza A virus {serotype_qualifier}"
                            else:
                                q_value = "Influenza A virus"
                        elif flu_type == "B":
                            q_value = "Influenza B virus"
                        else:
                            q_value = f"Influenza {flu_type} virus"
                        qualifier[4] = q_value
                    if q_key == "mol_type":
                        qualifier[4] = "viral cRNA"
                    if q_key == "ff_definition":
                        if "isolate" in q_value:
                            qualifier[4] = q_value.replace("@@[organism]@@ @@[isolate]@@ RNA,", f"@@[organism]@@ @@[isolate]@@ RNA, segment @@[segment]@@,").replace("genome", "sequence")
                        elif "strain" in q_value:
                            qualifier[4] = q_value.replace("@@[organism]@@ @@[strain]@@ RNA,", f"@@[organism]@@ @@[strain]@@ RNA, segment @@[segment]@@,").replace("genome", "sequence")
                        
                feature.append(["", "", "", "segment", segment])
        out_buffer.append(feature)

    if out_mss_file is None:
        out_mss_file =  mss_file
        logger.info(f"Overwriting original MSS file: {mss_file}")           
    with open(out_mss_file, "w") as f:
        out_str = "\n".join(["\t".join(row) for feature in out_buffer for row in feature])
        f.write(out_str + "\n")
    logger.info(f"Fixed MSS file for influenza virus: {out_mss_file}")

def parse_mdl_file(work_dir):
    """
    Parse MDL file to get segment names and Flu type and subgroup as dictionary
    key: model, value: (flu_type, segment, subgroup)
    example: AF387493 -> (B, 4, -)

    #                                     num   num   num
    #idx  model     group      subgroup  seqs  pass  fail
    #---  --------  ---------  --------  ----  ----  ----
    #---  --------  ---------  --------  ----  ----  ----
    1     AF387493  fluB-seg4  -            1     1     0
    2     AY191501  fluB-seg6  -            1     1     0
    3     AY504599  fluB-seg2  -            1     1     0
    4     AY504605  fluB-seg7  -            1     1     0
    5     AY504614  fluB-seg8  -            1     1     0
    6     EF626631  fluB-seg5  -            1     1     0
    7     EF626633  fluB-seg3  -            1     1     0
    8     EF626642  fluB-seg1  -            1     1     0
    #---  --------  ---------  --------  ----  ----  ----
    -     *all*     -          -            8     8     0
    -     *none*    -          -            0     0     0
    #---  --------  ---------  --------  ----  ----  ----
    """

    mdl_file = os.path.join(work_dir, "vadr", "vadr.vadr.mdl")
    if not os.path.exists(mdl_file):
        logger.error(f"MDL file not found: {mdl_file}")
        return None

    dict_model = {}
    with open(mdl_file) as f:
        in_table = False
        for line in f:
            if line.startswith("#") or line.startswith("-"):
                continue
            cols = line.strip().split()
            if len(cols) < 6:
                continue
            model = cols[1]
            group = cols[2]
            subgroup = cols[3]
            flu_type = group.split("-")[0] if "-" in group else None
            if flu_type:
                flu_type = flu_type.replace("flu", "").upper()  # "fluA" -> "A", "fluB" -> "B"
            else:
                logger.warning(f"Could not determine flu type from group: {group}")
                flu_type = "unknown_type"

            segment = group.split("-")[1] if "-" in group else None
            if segment:
                segment = segment.replace("seg", "")  # "seg1" -> "1"
            else:
                logger.warning(f"Could not determine segment from group: {group}")
                segment = "unknown_segment"            
            dict_model[model] = (flu_type, segment, subgroup)
    return dict_model

def parse_ftr_file(work_dir):
    """
    Parse vadr.ftr file to get seq_name and model name as dictionary
    Also check consistency of the model (confirm only one model is used for each segment)

    """
    ftr_file = os.path.join(work_dir, "vadr", "vadr.vadr.ftr")
    if not os.path.exists(ftr_file):
        logger.error(f"FTR file not found: {ftr_file}")
        return None
    dict_seq_model = {}
    with open(ftr_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split()
            seq_name = cols[1]
            model = cols[4]
            if seq_name in dict_seq_model:
                if dict_seq_model[seq_name] != model:
                    logger.error(f"Inconsistent model for sequence {seq_name}: {dict_seq_model[seq_name]} and {model}")
            dict_seq_model[seq_name] = model
    return dict_seq_model

def read_mss_file(mss_file):
    """
    Read MSS file and return list of lines
    """
    def iterEntries(fileName):
        Buffer = []
        for line in open(fileName):
            rows = line.strip("\n").split("\t")
            if rows[0] != "" and len(Buffer) > 0:
                yield Buffer
                Buffer = []
            Buffer.append(rows)
        yield Buffer

    def iterFeatures(entry):
        Buffer = []
        for rows in entry:
            if rows[1] != "" and len(Buffer) > 0:
                yield Buffer
                Buffer = []
            Buffer.append(rows)
        yield Buffer
    for entry in iterEntries(mss_file):
        for feature in iterFeatures(entry):
            yield feature

if __name__ == "__main__":
    # test
    logger.setLevel(logging.DEBUG)
    logger.addHandler(logging.StreamHandler())
    # fix_flu_mss("flu/test/", "DDBJ")
    fix_flu_mss("flu/test/", "A_Japan_NIG-xxxxx_2023", out_mss_file="flu/test/A_Japan_NIG-xxxxx_2023.fixed.annt.tsv")
    print("done")