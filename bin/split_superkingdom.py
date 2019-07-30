#!/usr/bin/env python3
#Thierry Janssens, 19JUL2019
#This script generates multiple fasta files of the generated sscaffold per superkingdom for 
#further downstream annotation.
#Usage:
#split_superkingdom.py sample_taxClassified.tsv sample_taxUnclassified.tsv sample_Bacteria.fasta sample_Viruse.fasta sample_Archaea.fasta
#To do: The generation of an output file for Eukaryotes is still switched off, since the downstream annotation is more complex""""

import pandas as pd
from sys import argv
from Bio import SeqIO

SCRIPT, INPUTFILE_TAX, INPUTSCAFFOLDS, PATH_TAXDUMP_RANKEDLINEAGE, OUTPUTFILE_BAC, OUTPUTFILE_VIR, OUTPUTFILE_ARCH = argv

TaxLCA = pd.read_csv(INPUTFILE_TAX, sep="\t", header=0)
TaxLCA.columns = TaxLCA.columns.str.replace('\s+','_')
colnames_rankedlineage=["tax_id","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"]
taxdump_rankedlineage = pd.read_csv(PATH_TAXDUMP_RANKEDLINEAGE, sep="|", header = None, names=colnames_rankedlineage, low_memory=False)

scaffolds_dict = {"scaffold_name":[], "scaffold_seq":[]}
for seq_record in SeqIO.parse(INPUTSCAFFOLDS, "fasta"):
    scaffolds_dict["scaffold_name"].append(seq_record.id)
    scaffolds_dict["scaffold_seq"].append(str(seq_record.seq))
scaffoldsFasta = pd.DataFrame.from_dict(scaffolds_dict)
df1 = pd.merge(TaxLCA, taxdump_rankedlineage, how = "left", left_on = "taxID", right_on = "tax_id").drop("tax_id", axis = 1)
df2 = pd.merge(df1, scaffoldsFasta, how = "left", left_on = "#queryID", right_on = "scaffold_name").drop("#queryID", axis = 1)

taxClassified=df2.loc[df2['taxID'].notnull()]
taxUnclassified = df2.loc[df2['taxID'].isnull()].drop(["taxID","tax_name","species","genus","family","order","class","phylum","kingdom","superkingdom"], axis = 1)
#is_Eukaryote = taxClassified.loc[taxClassified['superkingdom']=='Eukaryota']
is_Virus = taxClassified[taxClassified['superkingdom']=='Viruses']
is_Archaeon = taxClassified[taxClassified['superkingdom']=='Archaea']

is_Bacterium1 = taxClassified.loc[taxClassified['superkingdom'] =='Bacteria'] 
is_Bacterium2 = taxClassified.loc[taxClassified['taxID']==2]
is_Bacterium3 = taxClassified.loc[taxClassified['taxID']==131567]
is_Bacterium = is_Bacterium1.append(is_Bacterium2.append(is_Bacterium3))

bacteriafile = open(OUTPUTFILE_BAC, 'a')
for index, row in is_Bacterium.iterrows():
    print(">{}\n{}".format(row[11], row[12]), file = bacteriafile)
for index, row in taxUnclassified.iterrows():
    print(">{}\n{}".format(row[1], row[2]), file = bacteriafile)
bacteriafile.close()    

virusesfile = open(OUTPUTFILE_VIR, 'a')
for index, row in is_Virus.iterrows():
    print(">{}\n{}".format(row[11], row[12]), file = virusesfile)
virusesfile.close()  

#eukaryotafile=open("4311801221_Eukaryota_scaffolds.fasta", 'a')
#for index, row in is_Eukaryote.iterrows():
#    print(">{}\n{}".format(row[1], row[23]), file=eukaryotafile)
#eukaryotafile.close() 

archaeafile = open(OUTPUTFILE_ARCH, 'a')
for index, row in is_Archaeon.iterrows():
    print(">{}\n{}".format(row[11], row[12]), file = archaeafile)
archaeafile.close() 