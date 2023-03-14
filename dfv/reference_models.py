import os
import sys
import json
from copy import copy
# from dataclasses import dataclass

# @dataclass
class VadrModelFeature:

    def __init__(self, model_name, f_type, attrs):
        self.model_name = model_name
        self.type = f_type
        self.attrs = attrs
        self.hits = []

    def __repr__(self):
        if self.type is None:
            return f"<Model:{self.model_name}>"
        elif self.type == "gene":
            return f"{self.type}:{self.attrs['gene']}"
        elif "product" in self.attrs:
            return f"{self.type}:{self.attrs['product']}"
        else:
            return f"{self.type}:{self.attrs['coords']}"

VADR_MODEL_DIR = os.getenv("VADRMODELDIR")
VADR_SCOV2_MODELS_VERSION = os.getenv("VADR_SCOV2_MODELS_VERSION")

sar2_cov_2_minfo = os.path.join(VADR_MODEL_DIR, f"vadr-models-sarscov2-{VADR_SCOV2_MODELS_VERSION}", "sarscov2.minfo")



def parse_minfo_file(file_name):
    def attrs_to_dict(items):
        f_attrs = {}
        items = items.replace('" ', ('"@@@')).split('@@@')
        for item in items:
            key, value = item.split(":", 1)
            value = value.strip('"')
            f_attrs[key] = value
        # print(D)
        return f_attrs


    D = {}
    dict_count = {}
    for line in open(file_name):
        cols = line.strip().split(" ", 2)
        model_name = cols[1]
        f_attrs = attrs_to_dict(cols[2])
        # print(f_attrs)
        if cols[0] == "MODEL":
            f_type = None
            vadr_model_feature = VadrModelFeature(model_name, f_type, f_attrs)
            D[model_name] = [vadr_model_feature]
        elif cols[0] == "FEATURE":
            f_type = f_attrs.get("type")
            vadr_model_feature = VadrModelFeature(model_name, f_type, f_attrs)
            D[model_name].append(vadr_model_feature)
            model_name = cols[1]
            feature_info = attrs_to_dict(cols[2])
 
        else:
            sys.stderr.write(f"Error while reading model file {file_name}. Aborted.")
            raise AssertionError
    return D



def get_reference_model(name):
    # if name == "sars_cov_2":
    #     return copy(sars_cov_2)
    # else:
    #     raise NotImplementedError
    D = parse_minfo_file(sar2_cov_2_minfo)
    # return D["NC_045512-MW422255"]
    # print(D["NC_045512-MW422255"])
    return D

if __name__ == '__main__':
    D = parse_minfo_file(sar2_cov_2_minfo)
    # D = parse_minfo_file("/opt/vadr/vadr-models-flavi/dengue.minfo")
    for key, value in D.items():
        print(key, value)
    print("---")
    # print(D["NC_001477"])
    # print(D["NC_001477"][0].attrs["length"])
    # print(D["NC_045512-MW422255"])
    print([x for x in D["NC_045512"] if x.type =="CDS"])
    # print([x for x in D["NC_001477"] if x.type =="CDS"])