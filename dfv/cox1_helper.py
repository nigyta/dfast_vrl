#!/usr/bin/env python

import os
import sys
import argparse
import logging
import json
import re
from io import StringIO

logger = logging.getLogger(__name__)

def create_metadata_file(work_dir, args):
    metadata_file = args.metadata_file
    output_metadata_file = os.path.join(work_dir, "metadata.txt")
    metadata = {"projectType": "cox1"}

    common_metadata, specific_metadata = parse_common_and_specific(metadata_file)

    for key, value in common_metadata.items():
        # 1st layer: SUBMITTER, REFERENCE, COMMENT
        for key2, value2 in value.items():
            metadata[key2] = value2
    with open(output_metadata_file, "w", encoding="utf-8") as f:
        for key, value in metadata.items():
            f.write(f"{key}\t{value}\n")
    return output_metadata_file, common_metadata, specific_metadata

def parse_common_and_specific(filepath):
    """
    Parse a text file and extract the COMMON JSON block
    and the SPECIFIC TSV block.
    Returns a dictionary containing both blocks.
    """
    with open(filepath, encoding="utf-8") as f:
        text = f.read()

    # --- Extract COMMON block ---
    common_match = re.search(
        r"##COMMON.*?(\{.*?\})\s*##SPECIFIC",
        text,
        flags=re.DOTALL
    )
    if not common_match:
        raise ValueError("COMMON block not found")

    common_json_text = common_match.group(1)

    # Parse JSON (assumed to be valid)
    common_data = json.loads(common_json_text)
    # print(common_data)
    # logger.debug("Parsed COMMON data: %s", common_data)

    # --- Extract SPECIFIC block ---
    specific_match = re.search(
        r"##SPECIFIC.*?\n(.*)$",
        text,
        flags=re.DOTALL
    )
    if not specific_match:
        raise ValueError("SPECIFIC block not found")

    specific_text = specific_match.group(1).strip()
    # print("specific_text:", specific_text)
    # logger.debug("specific_text:", specific_text)
    # Parse TSV only if not empty

    specific_metadata = {}
    if specific_text:
        f = StringIO(specific_text)
        lines = f.readlines()
        header = lines[0].lstrip("# ").rstrip("\n").split("\t")
        for line in lines[1:]:
            fields = line.rstrip("\n").split("\t")
            # print("fields:", fields )

            D = {key: value for key, value in zip(header, fields)}
            entry = D.get("entry")
            if entry:
                specific_metadata[entry] = D
            else:
                logger.error("Missing 'entry' field in SPECIFIC metadata line: %s", line)
                exit(1)
        # print(header)
        # print("specific_metadata:", specific_metadata)
        # specific_metadata = parse_specific_tsv(specific_text)

    return common_data, specific_metadata

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Create metadata file for cox1_to_ddbj")
    parser.add_argument("-m", "--metadata_file", required=True, help="Input metadata file containing COMMON and SPECIFIC blocks")
    parser.add_argument("-o","--out_dir", required=True, help="Output directory to save the metadata file")
    args = parser.parse_args()

    work_dir = args.out_dir
    os.makedirs(work_dir, exist_ok=True)

    common_metadata, specific_metadata = create_metadata_file(work_dir, args)
    output_metadata_file = os.path.join(work_dir, "metadata.txt")

    print(f"Metadata file created at: {output_metadata_file}")
    print("COMMON metadata:")
    print(json.dumps(common_metadata, indent=4))
    print("SPECIFIC metadata:")
    for row in specific_metadata:
        print("\t".join(row))