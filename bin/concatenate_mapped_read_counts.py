#! /usr/bin/env python3

"""
Author: Sam Nooij
Date: 8 November 2019

Concatenate mapped read counts per sample into an overall
read count table. Also checks for duplicated scaffold 
names.

Requires input files with a file name that ends with:
-{sample}.tsv
where {sample} is the name of the sample.

Example use:
$ python3 concatenate_mapped_read_counts.py mapped_read_counts-1.tsv mapped_read_counts-2.tsv mapped_read_counts-3.tsv mapped_read_counts.tsv
"""

# IMPORT required libraries--------------------------------
import pandas as pd
import argparse
import os

# Define FUNCTIONS-----------------------------------------

def parse_arguments():
    """
    Parse the arguments from the command line, i.e.:
     -i/--input = list of input files (tab-separated tables)
     -o/--output = output file (tab-separated table)
     -h/--help = show help
    """
    parser = argparse.ArgumentParser(prog="concatenate mapped read counts",
             description="Concatenate mapped read count tables",
             usage="concatenate_mapped_read_counts.py -i [input] -o [output]"
             " [-h / --help]",
             add_help=False)

    required = parser.add_argument_group("Required arguments")

    required.add_argument('-i',
                          '--input',
                          dest="input",
                          metavar='',
                          required=True,
                          type=str,
                          nargs='+',
                          help="List of input files (counts per sample).")
    
    required.add_argument('-o',
                          '--output',
                          dest="output",
                          metavar='',
                          required=True,
                          type=str,
                          help="Output file name (and directory).")
    
    (args, extra_args) = parser.parse_known_args()

    return(args)

def main():
    """
    Main execution of the script
    """
    #1. Parse and show arguments
    arguments = parse_arguments()

    message =  ("\n"
                "These are the arguments you have provided:\n"
                "  INPUT:\n"
                "{0},\n"
                "  OUTPUT:\n"
                "{1}\n".format(arguments.input,
                               arguments.output))

    print(message)
    
    #2. Read input files and make into one dataframe
    concat_df = pd.DataFrame()
    
    for file in arguments.input:
        if file.count('-') > 1:
            #If there are multiple dashes in the file name, the
            # sample name contains one or more dashes.
            sample = os.path.splitext('-'.join(file.split('-')[1:]))[0]
            #So sample should include all parts separated by dashes.
        else:
            #Otherwise, the sample name has no dashes and splitting
            # the filename in two will work.
            sample = os.path.splitext(file.split('-')[-1])[0]
        df = pd.read_csv(file, sep='\t')
        df["Sample_name"] = sample
        
        concat_df = pd.concat([concat_df, df])
        
    #3. Check for duplicated scaffold names
    if len(concat_df["scaffold_name"]) > len(set(concat_df["scaffold_name"])):
        #If the whole list is longer than the deduplicated 'set',
        # there must be a duplicate in the list.
        print("Warning, duplicated scaffold name suspected!")
        print("Total number of scaffolds: %i" % len(concat_df["scaffold_name"]))
        print("Number of unique scaffold names: %i" % len(set(concat_df["scaffold_name"])))
        print("-------------------------------------")
        duplicates = concat_df[concat_df.duplicated(["scaffold_name"], keep = False)]
        print("These are the duplicates:\n\n", duplicates)
        print("\nAs contig names are unique _per sample_, "
            "you should also filter by *Sample_name* to "
            "separate scaffolds and count correctly.")
        print("(This is done automatically with "
            "'quantify_output.py')")
        
    else:
        #Otherwise, there should be no duplicates and everything
        # is fine.
        print("No duplicate scaffold names found. Proceed normally.")
        
    #4. Write table to a tsv file
    concat_df.to_csv(arguments.output, sep = '\t', 
                     index = False)

        

#EXECUTE script--------------------------------------------
if __name__ == "__main__":
    main()
