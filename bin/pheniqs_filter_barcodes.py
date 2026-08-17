#!/usr/bin/env python

import sys
import pandas as pd


observed_bc = sys.argv[1]
read_threshold = int(sys.argv[2])
cell_threshold = int(sys.argv[3])

df = pd.read_csv(observed_bc, sep='\s+', names=['number','bc'])
df1 = df[(df['number']>=read_threshold)&(~df['bc'].str.contains('='))].sort_values('number', ascending=False)
df1.head(cell_threshold).bc.str.replace('-', '').to_csv("5_filt_bc.txt",sep='\n',header=False,index=False)