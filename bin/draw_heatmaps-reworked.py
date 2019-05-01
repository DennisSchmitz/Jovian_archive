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

# IMPORT required libraries--------------------------------
import argparse
import numpy as np
import pandas as pd
from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, ColumnDataSource


# Set global VARIABLES-------------------------------------
RANKS = ["superkingdom", "phylum", "class", "order",
          "family", "genus", "species"]

PHAGE_FAMILY_LIST = ["Myoviridae", "Siphoviridae", "Podoviridae", "Lipothrixviridae", 
              "Rudiviridae", "Ampullaviridae", "Bicaudaviridae", "Clavaviridae", 
              "Corticoviridae", "Cystoviridae", "Fuselloviridae", "Globuloviridae", 
              "Guttaviridae", "Inoviridae", "Leviviridae", "Microviridae", 
              "Plasmaviridae", "Tectiviridae"]


# Define FUNCTIONS-----------------------------------------
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
             usage="draw_heatmaps.py -c -n -s -v -p -b -sq -st -vs -ps -bs -col"
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

    required.add_argument('-sq',
                          '--super-quantities',
                          dest="super_quantities",
                          metavar='',
                          required=True,
                          type=str,
                          help="Table with superkingdom quantities per sample")

    required.add_argument('-st',
                          '--stats',
                          dest="stats",
                          metavar='',
                          required=True,
                          type=str,
                          help="Table with taxonomic rank statistics")

    required.add_argument('-vs',
                          '--vir-stats',
                          dest="vir_stats",
                          metavar='',
                          required=True,
                          type=str,
                          help="Table with virual taxonomic rank statistics")

    required.add_argument('-ps',
                          '--phage-stats',
                          dest="phage_stats",
                          metavar='',
                          required=True,
                          type=str,
                          help="Table with phage taxonomic rank statistics")

    required.add_argument('-bs',
                          '--bact-stats',
                          dest="bact_stats",
                          metavar='',
                          required=True,
                          type=str,
                          help="Table with bacterial taxonomic rank statistics")

    optional = parser.add_argument_group("Optional arguments")

    optional.add_argument('-col',
                          '--colour',
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
    # Read the number of read pairs per read set/sample
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
    # Initialise the dataframe with taxonomic classifications
    # and numbers of reads mapped to the scaffolds (i.e.
    # the result/output of the pipeline).
    classifications_df = pd.read_csv(infile, delimiter='\t')

    # Check column names for debugging:
    #print(classifications_df.columns)

    # Select only relevant columns:
    classifications_df = classifications_df[[ "Sample_name", "#ID", "taxID", 
                                             "tax_name", "superkingdom", 
                                             "kingdom", "phylum", "class", 
                                             "order", "family", "genus",
                                             "species", "Plus_reads", 
                                             "Minus_reads", "Avg_fold", "Length" 
                                            ]]

    # Calculate the number of read pairs matched to each scaffold
    # by averaging the plus and minus reads.
    # N.B. This is an imperfect approximation.
    classifications_df["reads"] = round((classifications_df.Plus_reads +
                                       classifications_df.Minus_reads) / 2 )
    
    return(classifications_df)


def filter_taxa(df, taxon, rank):
    """
    Filter taxa of interest of a certain rank from
    a dataframe.
    (taxon may be a single taxon as string, or a list of taxa)
    """
    if isinstance(taxon, str):
        # If a string is provided, continue as intended
        subset_df = df[df[rank] == taxon]
    elif isinstance(taxon, list) and len(taxon) == 1:
        # If a single-entry list is provided, use taxon as string
        taxon = taxon[0]
        subset_df = df[df[rank] == taxon]
    else:
        # If a list is provided, filter all given taxa
        taxa_list = taxon
        subset_df = df[df[rank].isin(taxa_list)]

    return(subset_df)


def remove_taxa(df, taxon, rank):
    """
    Negative filter of taxa of a certain rank from
    a dataframe: remove them and keep the rest.
    (taxon may be a single taxon as string, or a list of taxa)
    """
    if isinstance(taxon, str):
        # If a string is provided, continue as intended
        subset_df = df[~df[rank] == taxon]
    elif isinstance(taxon, list) and len(taxon) == 1:
        # If a single-entry list is provided, use taxon as string
        taxon = taxon[0]
        subset_df = df[~df[rank] == taxon]
    else:
        # If a list is provided, filter all given taxa
        taxa_list = taxon
        subset_df = df[~df[rank].isin(taxa_list)]

    return(subset_df)

def report_taxonomic_statistics(df, outfile):
    """
    Input: dataframe with classifications of scaffolds, a name for an output file (txt)
    Output: a list of statistics in a text file, like:
        superkingdom 4
        phylum 50
        class 99
        order 220
        family 373
        genus 649
        species 337
    """
    header="taxonomic_level\tnumber_found\n"
    with open(outfile, 'w') as f:
        f.write(header)
        # Count how many taxa have been reported
        for t in [ "superkingdom", "phylum", "class", "order", "family", "genus", "species" ]:
            f.write("%s\t%i\n" % (t, df[t].nunique()))
    
    return(None)

def main():
    """
    Main execution of the script
    """
    #1. Parse and show arguments
    arguments = parse_arguments()

    message =  ("\n"
                "These are the arguments you have provided:\n"
                "  INPUT:\n"
                "classified = {0},\n"
                "numbers = {1}\n"
                "  OUTPUT:\n"
                "super_quantities = {2}\n"
                "stats = {3}\n"
                "vir_stats = {4}\n"
                "phage_stats = {5}\n"
                "bact_stats = {6}\n"
                "  OPTIONAL PARAMETERS:\n"
                "colour = {7}\n".format(arguments.classified,
                                        arguments.numbers,
                                        arguments.super_quantities,
                                        arguments.stats,
                                        arguments.vir_stats,
                                        arguments.phage_stats,
                                        arguments.bact_stats,
                                        arguments.colour))

    print(message)
    
    #2. Read input files and make dataframes
    numbers_df = read_numbers(arguments.numbers)
    classifications_df = read_classifications(arguments.classified)

    merged_df = classifications_df.merge(numbers_df, left_on="Sample_name", right_on="Sample")
    merged_df["Percentage"] = merged_df.reads / merged_df.read_pairs * 100

    print(merged_df.head()[["#ID", "reads", "read_pairs", "Percentage"]])
    print(merged_df.info())

    #3. Create chunks of information required for the heatmaps
    #3.1. Aggregate superkingdom-rank information
    # Count the percentages of Archaea, Bacteria, Eukaryota and Viruses per sample:
    superkingdom_sums = pd.DataFrame(merged_df.groupby(
                        [ "Sample_name", "superkingdom" ]).sum()
                        [[ "reads", "Percentage" ]])
    superkingdom_sums.reset_index(inplace=True) #to use MultiIndex "Sample_name" and "superkingdom" as columns
    
    missing_superkingdoms = { "Sample_name": [],
                            "superkingdom": [],
                            "reads": [],
                            "Percentage": []}
    
    # Check for missing taxa:
    for sample in set(superkingdom_sums["Sample_name"]):
        subset = superkingdom_sums.loc[superkingdom_sums.Sample_name == sample, ["superkingdom"]]
        for taxon in [ "Archaea", "Bacteria", "Eukaryota", "Viruses" ]:
            if taxon not in subset.values:
                missing_superkingdoms["Sample_name"].append(sample)
                missing_superkingdoms["superkingdom"].append(taxon)
                missing_superkingdoms["reads"].append(0)
                missing_superkingdoms["Percentage"].append(0)

    complete_superkingdoms = pd.concat([superkingdom_sums, pd.DataFrame(missing_superkingdoms)])
    complete_superkingdoms.sort_values(by=["Sample_name", "superkingdom"], inplace=True)
    complete_superkingdoms.reset_index(inplace=True)
    complete_superkingdoms["reads"] = complete_superkingdoms["reads"].astype(int)

    complete_superkingdoms.to_csv(arguments.super_quantities, index=False)

    print(complete_superkingdoms)
    
    #3.2. Filter viruses from the table
    virus_df = filter_taxa(df=merged_df, taxon="Viruses", rank="superkingdom")
    # Remove the phages from the virus df to make less cluttered heatmaps
    virus_df = remove_taxa(df=virus_df, taxon=PHAGE_FAMILY_LIST, rank="family")
    
    #3.3. Filter phages
    phage_df = filter_taxa(df=merged_df, taxon=PHAGE_FAMILY_LIST, rank="family")

    #3.4. Filter bacteria
    bacterium_df = filter_taxa(df=merged_df, taxon="Bacteria", rank="superkingdom")

    print(virus_df.head())
    print(phage_df.head())
    print(bacterium_df.head())

    #4. Write taxonomic rank statistics to a file, for each chunk
    #4.1. All taxa
    report_taxonomic_statistics(df = merged_df, outfile = arguments.stats)
    #4.2. Viruses
    report_taxonomic_statistics(df = virus_df, outfile = arguments.vir_stats)
    #4.3. Phages
    report_taxonomic_statistics(df = phage_df, outfile = arguments.phage_stats)
    #4.4. Bacteria
    report_taxonomic_statistics(df = bacterium_df, outfile = arguments.bact_stats)


#EXECUTE script--------------------------------------------
if __name__ == "__main__":
    main()
