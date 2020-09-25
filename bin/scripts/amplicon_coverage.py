#!/usr/bin/env python3
import pandas as pd
import subprocess
from sys import argv
import os

from pandas.core.dtypes import base


## import the primer bedfile as an argument
LOCAL, PRIMERBED = argv


df1 = pd.read_csv(PRIMERBED, sep='\t', header=None, usecols=[0, 1, 2, 3])
df1 = df1.applymap(str)
df1 = df1.rename(columns={3: "primer_pairs"})

df_regions = df1[["primer_pairs"]].copy()
df_regions['samtools_positions'] = df1[0] + ":" + df1[1] + "-" + df1[2]

### Walk over the folder with alignments and generate a list of bam files to use.
bampaths = []

for file in os.listdir("data/alignment/bam-files"):
    if file.endswith(".bam"):
        filepath = os.path.join("data/alignment/bam-files", file)
    
        bampaths.append(filepath)

bampaths.sort()

def execute_samtools(pos, infile):
    cmd = f"samtools coverage -H -r {pos} {infile}"
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True, shell=True)
    prc_output = process.stdout.readline()
    intermediate = prc_output.split('\t')
    tup = tuple(intermediate)
    avg_cov = tup[6]
    return avg_cov

for a in bampaths:
    base = os.path.basename(a)
    name = os.path.splitext(base)[0]
    df_regions[name] = df_regions['samtools_positions'].apply(execute_samtools, args=[a])

headers = ['primer_pairs']
for n in bampaths:
    b = os.path.basename(n)
    fname = os.path.splitext(b)[0]
    headers.append(fname)

df_regions.to_csv("results/fragment_coverage.tsv", sep='\t', columns = headers, index=False)