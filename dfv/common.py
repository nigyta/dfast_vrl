import os
import shutil
import logging
import sys


def update_metadata_file(metadata_file, seq_status, number_of_sequence, mol_type="RNA"):
    """
    Update metadata file depending on sequence status and number of sequence
    """

    lines = open(metadata_file).readlines()
    ret = ""
    for line in lines:
        key, value = line.strip("\n").split("\t", 1)
        if not (key in ["ff_definition", "seqStatus", "mol_type"]):
            ret += line
    
    ret += f"seqStatus\t{seq_status}\n"
    ret += f"mol_type\tgenomic {mol_type}\n"
    if seq_status == "complete":
        ret += f"ff_definition\t@@[organism]@@ @@[isolate]@@ {mol_type}, complete genome\n"
    elif seq_status == "nearly complete":
        ret += f"ff_definition\t@@[organism]@@ @@[isolate]@@ {mol_type}, nearly complete genome\n"
    elif seq_status == "draft":
        if number_of_sequence == 1:
            ret += f"ff_definition\t@@[organism]@@ @@[isolate]@@ {mol_type}, partial genome\n"
        else:
            ret += f"ff_definition\t@@[organism]@@ @@[isolate]@@ {mol_type}, draft genome, @@[entry]@@\n"
    else:
        raise AssertionError
    with open(metadata_file, "w") as f:
        f.write(ret)


def get_isolate(metadata_file, args):
    """
    Return isolate name, prefix for the DDBJ submission file
    priority:
    1. 'isolate' specified by command-line arg 
    2. isolate name described in metadata file
    3. If neither of above available, 'unspecified_isolate' is used
    """
    metadata = {}
    for line in open(metadata_file):
        key, value = line.strip("\n").split("\t")
        metadata[key] = value
    isolate_from_metadata_file = metadata.get("isolate")
    if args.isolate:
        isolate = args.isolate
        mss_file_prefix = isolate.replace("/", "_")
        metadata["isolate"] = isolate
        with open(metadata_file, "w") as f:
            for key, value in metadata.items():
                f.write(f"{key}\t{value}\n")
    elif isolate_from_metadata_file:
        isolate = isolate_from_metadata_file
        mss_file_prefix = isolate.replace("/", "_")
    else:
        isolate = "unspecified_isolate"
        mss_file_prefix = "DDBJ"
    return isolate, mss_file_prefix


def copy_or_create_metadata_file(work_dir, args):
    metadata_file_copy = os.path.join(work_dir, "metadata.txt")
    if args.metadata_file is None:
        # Create dummy metadata file
        with open(metadata_file_copy, "w") as f:
            f.write("projectType\tvrl\n")
    else:
        shutil.copy(args.metadata_file, metadata_file_copy)
    return metadata_file_copy


def get_logger(name=None, debug=False, work_dir="."):
    if debug:
        log_level = logging.DEBUG
    else:
        log_level = logging.INFO

    logger = logging.getLogger(__name__)
    sh = logging.StreamHandler(stream=sys.stdout)
    log_file = os.path.join(work_dir, "dfast_vrl.log")
    fh = logging.FileHandler(log_file, mode="w", encoding="utf-8", delay=True)

    logging.basicConfig(
        format="[%(asctime)s] [%(levelname)s] %(message)s",
        level=log_level,
        handlers=[sh, fh]) 
    logger = logging.getLogger(__name__)
    return logger