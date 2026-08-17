#!/usr/bin/env python

import pandas as pd
import json
import sys

split_report = sys.argv[1]

with open(split_report, 'r') as file:
    json_data = json.load(file)

def extract_barcode_info(data):
    barcode_info = []
    for entry in data.get('cellular', []):
        for classified_entry in entry.get('classified', []):
            barcode_info.append({
                'barcode': classified_entry['barcode'][0],
                'count': classified_entry['count'],
                'fastq_index': classified_entry['index']
            })
    return barcode_info


barcode_data = extract_barcode_info(json_data)
df = pd.DataFrame(barcode_data)
df['fastq_index'] = df['fastq_index'].apply(lambda x: str(x - 1).zfill(4))
df.to_csv("8_pheniqs_report.csv",header=False,index=False)

