#!/usr/bin/env python3
#Thierry Janssens, 19JUL2019
#This script generates multiple fasta files of the generated sscaffold per superkingdom for 
#further downstream annotation.
#Usage:
#split_superkingdom.py sample_taxClassified.tsv sample_taxUnclassified.tsv sample_Bacteria.fasta sample_Viruse.fasta sample_Archaea.fasta
#To do: The generation of an output file for Eukaryotes is still switched off, since the downstream annotation is more complex""""

import pandas as pd
from sys import argv

#SCRIPT, INPUTFILE_CLASSIFIED, INPUTFILE_UNCLASSIFIED, OUTPUTFILE_BAC, OUTPUTFILE_VIR, OUTPUTFILE_ARCH = argv

taxClassified = pd.read_csv(INPUTFILE_CLASSIFIED, sep ="\t")
taxUnclassified = pd.read_csv(INPUTFILE_UNCLASSIFIED, sep ="\t")
#is_Eukaryote = taxClassified.loc[taxClassified['superkingdom']=='Eukaryota']
is_Virus = taxClassified[taxClassified['superkingdom']=='Viruses']
is_Archaeon = taxClassified[taxClassified['superkingdom']=='Archaea']

is_Bacterium1 = taxClassified.loc[taxClassified['superkingdom'] =='Bacteria'] 
is_Bacterium2 = taxClassified.loc[taxClassified['taxID']==2]
is_Bacterium3 = taxClassified.loc[taxClassified['taxID']==131567]
is_Bacterium = is_Bacterium1.append(is_Bacterium2.append(is_Bacterium3))

bacteriafile=open(OUTPUTFILE_BAC, 'a')
for index, row in is_Bacterium.iterrows():
    print(">{}\n{}".format(row[1], row[23]), file=bacteriafile)
for index, row in taxUnclassified.iterrows():
    print(">{}\n{}".format(row[1], row[12]), file=bacteriafile)
bacteriafile.close()    

virusesfile=open(OUTPUTFILE_VIR, 'a')
for index, row in is_Virus.iterrows():
    print(">{}\n{}".format(row[1], row[23]), file=virusesfile)
virusesfile.close()  

#eukaryotafile=open("4311801221_Eukaryota_scaffolds.fasta", 'a')
#for index, row in is_Eukaryote.iterrows():
#    print(">{}\n{}".format(row[1], row[23]), file=eukaryotafile)
#eukaryotafile.close() 

archaeafile=open(OUTPUTFILE_ARCH, 'a')
for index, row in is_Archaeon.iterrows():
    print(">{}\n{}".format(row[1], row[23]), file=archaeafile)
archaeafile.close() 