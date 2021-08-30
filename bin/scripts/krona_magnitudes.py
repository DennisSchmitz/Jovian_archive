#!/usr/bin/env python3
"""
Thierry Janssens 01 July 2019
Add contigs withou blastn hits and magnitude information to the Krona chart. 

Usage:
  krona_magnitudes.py <input_taxtab> <input_stats> <output>

<input_taxtab> is the file generated by `ktClassifyBLAST` or by the mgkit function taxonutils lca.
<input_stats> is the file generated by `BBtools pileup.sh` script.
<output> is output taxMagtab; `ktImportTaxonomy` will use this to generate
a Krona chart where the size of the pie-parts are scaled to the number of 
aligned reads.

Example:
  krona_magnitudes.py data/taxonomic_classification/[sample_name].taxtab ata/taxonomic_classification/[sample_name].blastn data/scaffolds_filtered/[sample_name]_perMinLenFiltScaffold.stats data/taxonomic_classification/[sample_name].taxMagtab
"""

import pandas as pd
import numpy as np
import sys
from sys import argv


SCRIPT, INPUTTAX, INPUTSTATS, OUTPUTFILE = argv

df_tax = pd.read_csv(INPUTTAX, sep="\t")
fields = ["#ID", "Plus_reads", "Minus_reads"]

df_stats = pd.read_csv(INPUTSTATS, sep="\t", usecols=fields)
df_statsID = df_stats["#ID"]
df_statsID.rename("#queryID", inplace=True)
df_nohits = pd.DataFrame(df_statsID[~df_statsID.isin(df_tax[df_tax.columns[0]])])

df_nohits.insert(loc=1, column="taxID", value=0)
df_nohits.insert(2, "Avg. log e-value", "1")
df_tax = pd.concat([df_tax, df_nohits], join="inner")


df_stats["Nr_of_reads"] = df_stats["Plus_reads"].add(df_stats["Minus_reads"])
df_stats.drop(["Plus_reads", "Minus_reads"], axis=1, inplace=True)

merged_tax_readCount = pd.merge(
    df_tax, df_stats, how="left", left_on="#queryID", right_on="#ID"
).drop("#ID", axis=1)

merged_tax_readCount.to_csv(OUTPUTFILE, index=False, sep="\t")
