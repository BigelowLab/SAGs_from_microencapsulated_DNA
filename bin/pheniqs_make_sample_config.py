#!/usr/bin/env python

import sys
import json
from pathlib import Path
import os


sample_fastq_pth = sys.argv[1]
bc_D_pth = sys.argv[2]
bc_C_pth = sys.argv[3]
bc_B_pth = sys.argv[4]
bc_A_pth = sys.argv[5]
sample_hamming_dist = int(sys.argv[6])


full_dict = {
    "input": [sample_fastq_pth],
    "output": ["sample_demux.bam"],
    "cellular": [
        {
        "algorithm":"mdd",
        "base":"round_D",
        "codec": {},
        "distance tolerance":[sample_hamming_dist],
        "transform": {"token": ["0:0:8"]}
            },
        {"algorithm":"mdd",
        "base":"round_C",
        "codec": {},
        "distance tolerance":[sample_hamming_dist],
        "transform": {"token": ["0:12:20"]}
            },
        {"algorithm":"mdd",
        "base":"round_B",
        "codec": {},
        "distance tolerance":[sample_hamming_dist],
        "transform": {"token": ["0:24:32"]}
            },
        {"algorithm":"mdd",
        "base":"round_A",
        "codec": {},
        "distance tolerance":[sample_hamming_dist],
        "transform": {"token": ["0:36:44"]}
            }
            ],
    "template": {"transform": {"token": [ "0::"]}}
    
}

with open(bc_D_pth, "rt") as D:
    ls_D = list(line.strip() for line in D)
for idx, seq in enumerate(ls_D):
    full_dict["cellular"][0]["codec"][f"d{idx:02}"] = {"barcode": [f"{seq}"], "concentration": 1}

with open(bc_C_pth, "rt") as C:
    ls_C = list(line.strip() for line in C)
for idx, seq in enumerate(ls_C):
    full_dict["cellular"][1]["codec"][f"c{idx:02}"] = {"barcode": [f"{seq}"], "concentration": 1}

with open(bc_B_pth, "rt") as B:
    ls_B = list(line.strip() for line in B)
for idx, seq in enumerate(ls_B):
    full_dict["cellular"][2]["codec"][f"b{idx:02}"] = {"barcode": [f"{seq}"], "concentration": 1}

with open(bc_A_pth, "rt") as A:
    ls_A = list(line.strip() for line in A)
for idx, seq in enumerate(ls_A):
    full_dict["cellular"][3]["codec"][f"a{idx:02}"] = {"barcode": [f"{seq}"], "concentration": 1}


with open("1_sample_pheniqs_config.json", 'w') as fp:
    json.dump(full_dict, fp)