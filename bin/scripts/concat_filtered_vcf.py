#!/usr/bin/env python3
"""Concatenate the filtered vcf files.

Usage:
  concat_filtered_vcf.py <input> <output>

<input> is a "glob" of vcf files, e.g. `data/scaffolds_filtered/*_filtered.vcf`.
<output> is the output tsv (tab seperated values).

Example:
  concat_filtered_vcf.py data/scaffolds_filtered/*_filtered.vcf results/all_filtered_SNPs.tsv
"""

import pandas as pd
import glob
import os
from sys import argv

SCRIPT, VCF_FOLDER_GLOB, OUTPUT_FILE = argv

samples_list_MergedData = glob.glob(VCF_FOLDER_GLOB)

column_names_list = [
    "Contig_name",
    "Position",
    "Identifier",
    "Reference_base",
    "Alternative_base",
    "Quality",
    "Filter_status",
    "Info",
    "Sample_name",
]

MergedData_df = pd.concat(
    [
        pd.read_csv(
            f, sep="\t", comment="#", header=None, names=column_names_list
        ).assign(Sample_name=os.path.basename(f))
        for f in samples_list_MergedData
    ]
)

cols = list(MergedData_df.columns.values)
cols.pop(cols.index("Sample_name"))
MergedData_df = MergedData_df[["Sample_name"] + cols]

MergedData_df[
    ["Total_depth_of_coverage", "Allele_frequency", "Strand_bias", "DP4"]
] = MergedData_df["Info"].str.split(";", expand=True)
MergedData_df["Total_depth_of_coverage"] = (
    MergedData_df["Total_depth_of_coverage"].str.split("=").str[-1]
)
MergedData_df["Allele_frequency"] = (
    MergedData_df["Allele_frequency"].str.split("=").str[-1]
)
MergedData_df["Strand_bias"] = MergedData_df["Strand_bias"].str.split("=").str[-1]
MergedData_df["DP4"] = MergedData_df["DP4"].str.split("=").str[-1]
MergedData_df[
    [
        "DoC_forward_ref_allele",
        "DoC_reverse_ref_allele",
        "DoC_forward_non-ref_allele",
        "DoC_reverse_non-ref_allele",
    ]
] = MergedData_df["DP4"].str.split(",", expand=True)
MergedData_df.drop("Info", 1, inplace=True)
MergedData_df.drop("DP4", 1, inplace=True)
MergedData_df["Total_depth_of_coverage"] = MergedData_df[
    "Total_depth_of_coverage"
].astype("int")
MergedData_df["Allele_frequency"] = MergedData_df["Allele_frequency"].astype("float")
MergedData_df["Strand_bias"] = MergedData_df["Strand_bias"].astype("int")
MergedData_df["DoC_forward_ref_allele"] = MergedData_df[
    "DoC_forward_ref_allele"
].astype("int")
MergedData_df["DoC_reverse_ref_allele"] = MergedData_df[
    "DoC_reverse_ref_allele"
].astype("int")
MergedData_df["DoC_forward_non-ref_allele"] = MergedData_df[
    "DoC_forward_non-ref_allele"
].astype("int")
MergedData_df["DoC_reverse_non-ref_allele"] = MergedData_df[
    "DoC_reverse_non-ref_allele"
].astype("int")

MergedData_df.to_csv(OUTPUT_FILE, index=False, sep="\t")
