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

## Define FUNCTIONS----------------------------------------
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

    (args, extra_args) = parser.parse_known_args()

    return(args)


def main():
    """
    Main execution of the script
    """
    arguments = parse_arguments()

    message =  ("\n"
                "These are the arguments you have provided:\n"
                "classified = {0},\n"
                "numbers = {1}".format(arguments.classified,
                             arguments.numbers))

    print(message)

if __name__ == "__main__":
    main()