import os
import shutil
# from logging import getLogger, StreamHandler, FileHandler, Formatter, DEBUG, INFO
import sys


def update_metadata_file(metadata_file, seq_status, number_of_sequence, mol_type="RNA", organism=None):
    """
    Update metadata file depending on sequence status and number of sequence
    """

    lines = open(metadata_file).readlines()
    D = {}  # for metadata
    for line in lines:
        key, value = line.strip("\n").split("\t", 1)
        D[key] = value

    # Add organism if not exists
    if "organism" not in D and organism:
        D["organism"] = organism

    ret = ""
    for key, value in D.items():
        if not (key in ["ff_definition", "seqStatus", "mol_type"]):
            ret += f"{key}\t{value}\n"

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


# def get_logger(name=None, debug=False, work_dir="."):
#     if debug:
#         log_level = DEBUG
#     else:
#         log_level = INFO
#     logger = getLogger(name)
#     logger.setLevel(log_level)
#     log_format = "[%(asctime)s] [%(levelname)s] %(message)s"
#     sh = StreamHandler(stream=sys.stdout)
#     sh.setFormatter(Formatter(log_format))
#     sh.setLevel(log_level)
#     log_file = os.path.join(work_dir, "dfast_vrl.log")
#     fh = FileHandler(filename=log_file, mode="w", encoding="utf-8", delay=True)
#     fh.setFormatter(Formatter(log_format))
#     fh.setLevel(log_level)
#     logger.addHandler(fh)
#     logger.addHandler(sh)
#     logger.debug("Creating log file: %s, %s", log_file, name)


#     return logger

