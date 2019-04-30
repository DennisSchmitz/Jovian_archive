#! /usr/bin/env python

# coding: utf-8

# # Create a heatmap of taxa in all analysed samples, quantified by the number of reads per sample
# 
# _Date: 13 Nov 2018_  
# _Author: Sam Nooij_
# 
# Input: Identified taxa, number of reads (or read pairs) per taxon, total number of reads per sample
# 
# Output: Heatmap of taxa in different samples from the same run, quantified by their read numbers to look at the differences between samples
# 
# -- 14 Nov 2018 update: changed notebook into regular Python script that is run by the Snakefile
# -- 30 Apr 2019 update: start complete rework to make the script snakemake-independent and fix bugs
# 
# Required Python packages:
#  - Pandas
#  - Bokeh
#  - argparse

#IMPORT required libraries---------------------------------
import sys
import argparse
import numpy as np
import pandas as pd
from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, ColumnDataSource

#Define FUNCTIONS------------------------------------------
def parse_arguments():
    """
    Parse the arguments from the command line, i.e.:
     -c/--classified = table with taxonomic classifications
     -n/--numbers = file with (MultiQC, Trimmomatic's) read numbers
     -s/--super = output files for superkingdoms heatmap (as list)
     -v/--virus = output files for virus heatmaps (as list)
     -p/--phage = output files for phage heatmaps (as list)
     -b/--bact = output files for bacteria heatmaps (as list)
     -sq/--super-quantities = output file for superkingdom quantities
     -st/--stats = output file with taxonomic rank statistics
     -vs/--vir-stats = ouput file for viral statistis
     -ps/--phage-stats = output file for phage statistics
     -bs/--bact-stats = output file for bacterial statistics
     -col/--colour = heatmap colour
     -h/--help = show help
    """
    parser = argparse.ArgumentParser(prog="draw heatmaps",
             description="Draw heatmps for the Jovian taxonomic output",
             usage="draw_heatmaps.py -c -n -s -v -p -b -sq -st -vs -ps -bs"
             "[-h / --help]",
             add_help=False)


    required = parser.add_argument_group("Required arguments")

    required.add_argument('-c',
                          '--classified',
                          dest="classified",
                          metavar='',
                          required=True,
                          type=str,
                          help="Table with taxonomic classifications.")

    required.add_argument('-n',
                          '--numbers',
                          dest="numbers",
                          metavar='',
                          required=True,
                          type=str,
                          help="Multiqc Trimmomatic file with read numbers")

    optional = parser.add_argument_group("Optional arguments")

    optional.add_argument('-col',
                          "--colour",
                          dest="colour",
                          metavar='',
                          required=False,
                          type=str,
                          nargs='+',
                          default=[ "#000000" ],
                          help="Colour of the heatmap tiles")

    (args, extra_args) = parser.parse_known_args()

    return(args)

def read_numbers(infile):
    """
    Input: Tabular text file (.tsv) with number of reads/read pairs per sample
    Output: Pandas Dataframe with sample names and numbers in columns
    """
    #Read the number of read pairs per read set/sample
    numbers_df = pd.read_csv(infile, delimiter='\t')
    numbers_df = numbers_df[[ "Sample", "input_read_pairs" ]]
    numbers_df = numbers_df.rename(columns={"input_read_pairs" : "read_pairs"})

    numbers_df["Sample"] = numbers_df.Sample.apply(lambda x: x[:x.rfind("_R1")]) # On every value in column named "Sample" perform function that chops off "_R1" and any character after it
    
    return(numbers_df)

def read_classifications(infile):
    """
    Input: Tabulers text file (.tsv) with output from PZN analysis: 
      classifications for scaffolds and quantitative information of mapped-back reads,
      for _all samples analysed in the same run_
    Output: Pandas Dataframe with the information of the classified scaffolds
    """
    #Initialise the dataframe with taxonomic classifications
    # and numbers of reads mapped to the scaffolds (i.e.
    # the result/output of the pipeline).
    classifications_df = pd.read_csv(infile, delimiter='\t')

    #Check column names for debugging:
    #print(classifications_df.columns)

    #Select only relevant columns:
    classifications_df = classifications_df[[ "Sample_name", "#ID", "taxID", 
                                             "tax_name", "superkingdom", 
                                             "kingdom", "phylum", "class", 
                                             "order", "family", "genus",
                                             "species", "Plus_reads", 
                                             "Minus_reads", "Avg_fold", "Length" 
                                            ]]

    #Calculate the number of read pairs matched to each scaffold
    # by averaging the plus and minus reads.
    #N.B. This is an imperfect approximation.
    classifications_df["reads"] = round((classifications_df.Plus_reads +
                                       classifications_df.Minus_reads) / 2 )
    
    return(classifications_df)

def main():
    """
    Main execution of the script
    """
    arguments = parse_arguments()

    message =  ("\n"
                "These are the arguments you have provided:\n"
                "  INPUT:\n"
                "classified = {0},\n"
                "numbers = {1}\n"
                "  OPTIONAL PARAMETERS:\n"
                "colour = {2}\n".format(arguments.classified,
                             arguments.numbers,
                             arguments.colour))

    print(message)
    
    numbers_df = read_numbers(arguments.numbers)
    classifications_df = read_classifications(arguments.classified)

    merged_df = classifications_df.merge(numbers_df, left_on="Sample_name", right_on="Sample")
    merged_df["Percentage"] = merged_df.reads / merged_df.read_pairs * 100

    print(merged_df.head()[["#ID", "reads", "read_pairs", "Percentage"]])
    print(merged_df.info())
    

#EXECUTE script--------------------------------------------
if __name__ == "__main__":
    main()
