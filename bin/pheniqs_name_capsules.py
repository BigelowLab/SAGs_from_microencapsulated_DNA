#!/usr/bin/env python

import json
import os
import sys

PATH_filt_bc = sys.argv[1]
PATH_r1 = sys.argv[2]
PATH_r2 = sys.argv[3]
split_hamming_dist = int(sys.argv[4])
library = sys.argv[5]
PATH_bc_D = sys.argv[6]
PATH_bc_C = sys.argv[7]
PATH_bc_B = sys.argv[8]
PATH_bc_A = sys.argv[9]

DICT_D = {}
with open(PATH_bc_D, "rt") as D: LIST_D = list(line.strip() for line in D)
for idx, seq in enumerate(LIST_D): DICT_D[f"{seq}"] = f"{idx:02}"

DICT_C = {}
with open(PATH_bc_C, "rt") as C: LIST_C = list(line.strip() for line in C)
for idx, seq in enumerate(LIST_C): DICT_C[f"{seq}"] = f"{idx:02}"

DICT_B = {}
with open(PATH_bc_B, "rt") as B: LIST_B = list(line.strip() for line in B)
for idx, seq in enumerate(LIST_B): DICT_B[f"{seq}"] = f"{idx:02}"

DICT_A = {}
with open(PATH_bc_A, "rt") as A: LIST_A = list(line.strip() for line in A)
for idx, seq in enumerate(LIST_A): DICT_A[f"{seq}"] = f"{idx:02}"

#### II. Assign an ID to each barcode combination -> DICT_combo_to_id ##### __________________________________
'''
E.g. a 'GGTCATTGAAGGTGACTACAACCGGTCCGAT' gets ID'd as "19140024" because:
* barcode D positions 0-8 ("GGTCATTG") match the 19th barcode in bcD_24.txt -> 19
* barcode C pos' 8-16 match 14th in bcC_24.txt -> 14
* barcode B pos' 16-24 match 00th in BcB_24.txt -> 00
* barcode A pos' 24-32 match 24th in BcD_24.txt -> 24
'''

DICT_combo_to_id = {}

with open(PATH_filt_bc) as file:
    for line in file:
        id = DICT_D[line[0:8]]+DICT_C[line[8:16]]+DICT_B[line[16:24]]+DICT_A[line[24:32]]
        barcode = line[0:32]
        DICT_combo_to_id[barcode] = id

#### III. Write the config file for Pheniqs ##### __________________________________

DICT_config = {
    "input":[PATH_r1,PATH_r2],
    "output": [f"undetermined_r1.fastq.gz",f"undetermined_r2.fastq.gz"],

    "cellular": [
        {
        "algorithm":"mdd",
        "base":"round_A",
        "codec": {},
        "distance tolerance":[split_hamming_dist],
        "transform": {"token": ["1:0:8","1:12:20","1:24:32","1:36:44"],"knit": ["0:1:2:3"]}
            },
            ],
    "template": {"transform": {"token": [ "0::", "1::"]}}
    
}

with open(PATH_filt_bc, "rt") as A:
    ls_A = list(line.strip() for line in A)
for idx, seq in enumerate(ls_A):
    DICT_config["cellular"][0]["codec"][DICT_combo_to_id[seq]] = {"barcode": [f"{seq}"], "concentration": 1,"output":[f"{library}_{DICT_combo_to_id[seq]}_r1.fastq.gz",f"{library}_{DICT_combo_to_id[seq]}_r2.fastq.gz"]}

with open("6_split_pheniqs_config.json", 'w') as fp:
    json.dump(DICT_config, fp)
