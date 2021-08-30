#! /usr/bin/env python
# coding: utf-8

"""
Author: Sam Nooij
Date: 4 December 2018

-- 21 May 2019 update: start complete rework to make the script snakemake-independent and fix bugs
        N.B. For now, the number of threads is still passed by snakemake!
          It can optionally be provided as command-line argument.
-- 19 Sep 2019 update: count mapped reads with separate script (with samtools view)
        merge these reads to all_tax[Un]classified.tsv for correct quantifications

Script that gathers numbers from the PZN analysis, i.e.:
 - number of reads per sample (raw)
 - number of reads removed by trimmomatic
 - number of human reads identified by Bowtie2 (and removed)
 - number of reads mapped to scaffolds that have been classified as:
   - Archaea
   - Bacteria
   - Eukaryota
   - Viruses
 - number of reads that map to unclassified scaffolds ("dark matter")
(If there is a gap between these groups and the total number of raw reads,
 these were unable to assemble into contigs > 500 nt.)
   
These 7 numbers are collected for each sample and written to a table (.csv file)
and a stacked bargraph is made to visualise these data.

Input:
    Requires output files from FastQC, Trimmomatic, Bowtie, and the PZN taxonomic report:
     - results/multiqc_data/multiqc_fastqc.txt
     - results/multiqc_data/multiqc_trimmomatic.txt
     - data/HuGo_removal/*.fq
     - results/all_taxClassified.tsv
     - results/all_taxUnclassified.tsv
     
Output:
    Creates two tables (number of reads and percentages), and stacked bar charts:
     - results/read_counts.csv
     - results/profile_percentages.csv
     - results/Sample_composition_graph.html
"""

### Import required libraries ------------------------------------
import argparse
import datetime
import os
import sys
import re
import pandas as pd
from bokeh.core.properties import value
from bokeh.io import save, output_file
from bokeh.plotting import figure
from bokeh.models.widgets import Tabs, Panel
from bokeh.models import ColumnDataSource
import concurrent.futures

### Define global variables --------------------------------------
FASTQ_COUNT_WARINING_MESSAGE = """
Now counting the number of lines in the fastq files from which
reads have been discarded that were mapped to the human genome 
by Bowtie2 to determine the number of non-human reads.

This may take a while...
"""

COLOURS = [
    "#FFDBE5",
    "#7A4900",
    "#0000A6",
    "#63FFAC",
    "#B79762",
    "#004D43",
    "#8FB0FF",
    "#997D87",
]
# Colour scheme with distinct colours thanks to Tatarize
# (https://godsnotwheregodsnot.blogspot.com/2013/11/kmeans-color-quantization-seeding.html)
# And Alexey Popkov (https://graphicdesign.stackexchange.com/revisions/3815/8)

### Create functions that do all the work ------------------------
def parse_arguments():
    """
    Parse the arguments from the command line, i.e.:
     -f/--fastqc = file with (multiqc) fastqc output
     -t/--trimmomatic = file with (multiqc) trimmomatic output
     -hg/--hugo = fastq files of extracted human reads
     -c/--classified = table with taxonomic classifications
     -u/--unclassified = table with unclassified contigs
     -m/--mapped_reads = table with number of mapped reads
     -co/--counts = output file (table) with counts
     -pg/--percentages = output file (table) with percentages
     -cpu/--cpu-cores = number of cores (threads) to use
     -col/--colours = colours to use in figure/barchart (8 colours)
     -h/--help = show help
    """
    parser = argparse.ArgumentParser(
        prog="draw heatmaps",
        description="Draw heatamps for the Jovian taxonomic output",
        usage="quantify_profiles.py -f -t -hg -c -u" " [-h / --help]",
        add_help=False,
    )

    required = parser.add_argument_group("Required arguments")

    required.add_argument(
        "-f",
        "--fastqc",
        dest="fastqc",
        metavar="",
        required=True,
        type=str,
        help="MultiQC-FastQC output file",
    )

    required.add_argument(
        "-t",
        "--trimmomatic",
        dest="trimmomatic",
        metavar="",
        required=True,
        type=str,
        help="MultiQC-Trimmomatic output file",
    )

    required.add_argument(
        "-hg",
        "--hugo",
        dest="hugo",
        metavar="",
        required=True,
        type=str,
        nargs="+",
        help="Fastq files of extracted human reads",
    )

    required.add_argument(
        "-c",
        "--classified",
        dest="classified",
        metavar="",
        required=True,
        type=str,
        help="Table with taxonomic classifications",
    )

    required.add_argument(
        "-u",
        "--unclassified",
        dest="unclassified",
        metavar="",
        required=True,
        type=str,
        help="Table with unclassified contigs",
    )

    required.add_argument(
        "-m",
        "--mapped-reads",
        dest="mapped_reads",
        metavar="",
        required=True,
        type=str,
        help="Table with mapped read counts",
    )

    required.add_argument(
        "-co",
        "--counts",
        dest="counts",
        metavar="",
        required=True,
        type=str,
        help="Table of read counts",
    )

    required.add_argument(
        "-p",
        "--percentages",
        dest="percentages",
        metavar="",
        required=True,
        type=str,
        help="Table of read counts as percentages",
    )

    required.add_argument(
        "-g",
        "--graph",
        dest="graph",
        metavar="",
        required=True,
        type=str,
        help="Graph of sample compositions",
    )

    optional = parser.add_argument_group("Optional arguments")

    optional.add_argument(
        "-l",
        "--log",
        dest="log",
        metavar="",
        required=False,
        default=False,
        type=str,
        help="Log file to write warnings to",
    )

    optional.add_argument(
        "-cpu",
        "--cpu-cores",
        dest="cores",
        metavar="",
        required=False,
        type=int,
        default=4,
        help="Number of threads to read fastq files",
    )

    optional.add_argument(
        "-col",
        "--colours",
        dest="colours",
        metavar="",
        required=False,
        type=str,
        nargs="+",
        default=[
            "#FFDBE5",
            "#7A4900",
            "#0000A6",
            "#63FFAC",
            "#B79762",
            "#004D43",
            "#8FB0FF",
            "#997D87",
        ],
        help="Colours (8) of the barchart (stacks)",
    )
    # Colour scheme with distinct colours thanks to Tatarize
    # (https://godsnotwheregodsnot.blogspot.com/2013/11/kmeans-color-quantization-seeding.html)
    # And Alexey Popkov (https://graphicdesign.stackexchange.com/revisions/3815/8)

    optional.add_argument(
        "-h", "--help", action="help", help="Show this message and exit."
    )

    (args, extra_args) = parser.parse_known_args()

    return args


def count_sequences_in_fastq(infile):
    """
    Input: fastq file
    Output: line number / 4 (number of sequences)
    """
    with open(infile, "r") as f:
        for i, l in enumerate(f):
            pass
    try:
        lines = i + 1
        return lines / 4
    except UnboundLocalError:  # happens when file is empty, there is no 'i'
        return 0


def progress(count, total, status=""):
    """
    Copyright (c) Vladimir Ignatev (MIT Licence)
    https://gist.github.com/vladignatyev/06860ec2040cb497f0f3
    """
    bar_len = 60
    filled_len = int(round(bar_len * count / float(total)))

    percents = round(100.0 * count / float(total), 1)
    bar = "=" * filled_len + "-" * (bar_len - filled_len)

    sys.stdout.write(
        "[%s] %s%s - %s  [ %i / %i ]\r" % (bar, percents, "%", status, count, total)
    )
    # Adapted above line to also show how many files are counted - Sam
    sys.stdout.flush()  # As suggested by Rom Ruben (see: http://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console/27871113#comment50529068_27871113)


def count_non_human_reads(file_list, threads):
    """
    Input: list of fastq files after human read filtering
    Output: counts of non-human reads per sample (as pandas dataframe)
    
    For each sample, count the reads that did not map to the human genome:
     - forward read/mate
     - reverse read/mate
     - unpaired read
     (all summed)
    """
    print(FASTQ_COUNT_WARINING_MESSAGE)

    files_read = 0  # keep track of how many files have been read
    read_counts = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        for file, reads in zip(
            file_list, executor.map(count_sequences_in_fastq, file_list)
        ):
            if "_pR1.fq" in file:
                sample = file[: file.index("_pR1.fq")].split("/")[-1]
            elif "_pR2.fq" in file:
                sample = file[: file.index("_pR2.fq")].split("/")[-1]
            elif "_unpaired.fq" in file:
                sample = file[: file.index("_unpaired.fq")].split("/")[-1]
            else:
                print(
                    "A fastq file with has been provided for non-human read counting that does not match the expected suffix patterns (i.e. '_pR1.fq', '\_pR2.fq', or '_unpaired.fq'."
                )
                break  # stop if there is a file with an unexpected name

            if sample in read_counts:
                read_counts[sample] += reads
            else:
                read_counts[sample] = reads

            files_read += 1
            # Show the progress of counting files on the terminal:
            progress(files_read, len(file_list), status="Reading files")

    read_df = pd.DataFrame(
        {
            "Sample": [n for n in read_counts.keys()],
            "Non-human_reads": [n for n in read_counts.values()],
        }
    )

    read_df["Non-human_reads"] = read_df["Non-human_reads"].astype(int)

    print("Done counting!")

    return read_df


def sum_superkingdoms(classified_file, mapped_reads_file):
    """
    Input:
        taxonomic classifications and quantifications as
        generated by Jovian (i.e. results/all_taxClassified.tsv)
    Output:
        dataframe with the number of reads assigned to each
        superkindgom per sample, also taking into account reads
        that were not assigned a superkingdom ("not classified")
    """
    clas_df = pd.read_csv(classified_file, delimiter="\t")

    counts_df = pd.read_csv(mapped_reads_file, delimiter="\t")
    clas_df = pd.merge(
        clas_df, counts_df, how="left", on=["scaffold_name", "Sample_name"]
    )
    # If a scaffold has no reads mapping to it, set to 0:
    clas_df.fillna({"mapped_reads": 0}, inplace=True)

    # N.B. Take into account reads that were not assigned a superkingdom
    # and are therefore actually still unclassified!
    clas_df.fillna("not classified", inplace=True)

    # Count the reads assigned to Archaea, Bacteria, Eukaryota and Viruses per sample:
    superkingdom_sums = pd.DataFrame(
        clas_df.groupby(["Sample_name", "superkingdom"]).sum()[["mapped_reads"]]
    )

    # Check all samples and superkingdoms to find out for which superkingdoms there is no data
    samples_superkingdoms_dict = {}
    for i in superkingdom_sums.index:
        if i[0] not in samples_superkingdoms_dict.keys():
            # If the sample name is not yet in the dictionary, add it with the corresponding superkingdom
            samples_superkingdoms_dict[i[0]] = [i[1]]
        else:
            # If the sample is there already, add the current superkingdom
            samples_superkingdoms_dict[i[0]].append(i[1])

    # Make a dataframe for the missing data
    newdata = pd.DataFrame(index=superkingdom_sums.index, columns=[])

    for key, value in samples_superkingdoms_dict.items():
        for superkingdom in [
            "Archaea",
            "Bacteria",
            "Eukaryota",
            "Viruses",
            "not classified",
        ]:
            if superkingdom not in value:
                newdata.loc[(key, superkingdom), "mapped_reads"] = 0
            else:
                pass

    # Concatenate the missing data into the dataframe
    new_df = pd.concat([superkingdom_sums, newdata], sort=False)
    new_df.reset_index(inplace=True)
    new_df = new_df.sort_values(by=["Sample_name", "superkingdom"])
    new_df = new_df.dropna()
    new_df.reset_index(inplace=True, drop=True)

    # And make these into a proper dataframe:
    final_df = pd.DataFrame(
        new_df.pivot(index="Sample_name", columns="superkingdom", values="mapped_reads")
    )
    final_df.reset_index(inplace=True)

    return final_df


def sum_unclassified(unclassified_file, mapped_reads_file):
    """
    Input:
        table of scaffolds that were not classified by PZN
        (results/all_taxUnclassified.tsv)
    Output:
        dataframe with the number of reads assigned to 
        unclassified scaffolds per sample
    """
    unclas_df = pd.read_csv(unclassified_file, delimiter="\t")

    counts_df = pd.read_csv(mapped_reads_file, delimiter="\t")
    unclas_df = pd.merge(
        unclas_df, counts_df, how="left", on=["scaffold_name", "Sample_name"]
    )

    # Summarise the reads per sample:
    unclas_sums = pd.DataFrame(unclas_df.groupby("Sample_name").sum()[["mapped_reads"]])
    unclas_sums.reset_index(inplace=True)

    return unclas_sums


def validate_numbers(df, log=False):
    """
    Validates if the numbers add up to a number lower than the total 
    number of raw reads. (I.e. low-quality + human + other taxa +
    unclassified <= total sequences.) Reports a warning when the sum
    of these groups is too high.
    Input:
        Dataframe with all numbers of reads - profile of:
         - total/raw reads
         - low-quality reads
         - human reads
         - reads mapped to classified contigs:
             - Archaea
             - Bacteria
             - Eukaryota
             - Viruses
    Output:
        None if the numbers are within expected margins,
        Warning when the numbers are clearly wrong.
    """
    if log:
        with open(log, "a") as logfile:
            logfile.write("---\n\nNow checking if numbers add up properly...\n")
    else:
        print("Now checking if the numbers add up properly...")

    errors = []
    for i in range(len(df)):
        sample = df.loc[i, "Sample"]
        total = df.loc[i, "Total_reads"]
        reads_sum = (
            df.loc[i, "Low-quality"]
            + df.loc[i, "Human"]
            + df.loc[i, "Archaea"]
            + df.loc[i, "Bacteria"]
            + df.loc[i, "Eukaryota"]
            + df.loc[i, "Viruses"]
            + df.loc[i, "Unclassified"]
        )

        if not total >= reads_sum:
            if log:
                with open(log, "a") as logfile:
                    logfile.write(
                        "Warning! Sample %s has a bad reads sum: %i > total (%i)\n"
                        % (df.loc[i, "Sample"], reads_sum, total)
                    )
            else:
                print(
                    "Warning! Sample %s has a bad reads sum: %i > total (%i)"
                    % (df.loc[i, "Sample"], reads_sum, total)
                )
            errors.append(df.loc[i, "Sample"])
        else:
            pass

    if log:
        with open(log, "a") as logfile:
            if len(errors) == 0:
                logfile.write("All numbers okay: no sums are higher than the total!\n")
            else:
                logfile.write(
                    "%i errors have been found, these are from samples: %s\n"
                    % (len(errors), errors)
                )
    else:
        if len(errors) == 0:
            print("All numbers okay: no sums are higher than the total!\n")
        else:
            print(
                "%i errors have been found, these are from samples: %s"
                % (len(errors), errors)
            )

    return None


def draw_stacked_bars(df, perc, sample="", parts=[], outfile="", colours=COLOURS):
    """
    Takes a file of quantified read annotations by a pipeline (e.g. PZN)
    and draws a stacked bar chart using Bokeh, creating an interactive
    HTML file.
    
    Input:
     - A Pandas dataframe with all the required data
     - The column name of Sample IDs (as string)
     - The column names for the numbers ("parts", as list)
     - The title of the figure (as string)
     - The name of the output file (as string)
     - The colours to be used for the bars (as list of strings - hexcodes)
    Output:
     - Interactive (Bokeh) stacked bar chart that visualises
      the composition of each sample in your experiment
    """
    # Set a comfortable length depending on the number of samples to show:
    if len(df[sample]) > 40:
        width = len(df[sample]) * 20
    # Have the graph grow horizontally to accomodate large numbers of samples
    elif len(df[sample]) > 25:
        width = len(df[sample]) * 33
    else:
        # With a minimum of 800 to display the legend
        width = 800

    nr_fig = figure(
        x_range=df[sample],
        plot_height=800,
        plot_width=width,
        title="Composition of samples (read-based)",
        toolbar_location=None,
        tools="hover, pan",
        tooltips="@%s $name: @$name" % sample,
    )

    nr_fig.vbar_stack(
        parts,
        x=sample,
        width=0.9,
        color=colours,
        source=df,
        legend_label=str([value(x) for x in parts]),
    )

    nr_fig.y_range.start = 0
    nr_fig.x_range.range_padding = 0.1
    nr_fig.xgrid.grid_line_color = None
    nr_fig.axis.minor_tick_line_color = None
    nr_fig.outline_line_color = None
    nr_fig.legend.location = "top_left"
    nr_fig.legend.orientation = "horizontal"
    nr_fig.xaxis.major_label_orientation = 1

    nr_panel = Panel(child=nr_fig, title="Absolute number of reads")

    perc_fig = figure(
        x_range=perc[sample],
        plot_height=800,
        plot_width=width,
        title="Composition of samples (percentages)",
        toolbar_location=None,
        tools="hover, pan",
        tooltips="@%s $name: @$name" % sample,
    )

    perc_fig.vbar_stack(
        parts,
        x=sample,
        width=0.9,
        color=colours,
        source=perc,
        legend_label=str([value(x) for x in parts]),
    )

    perc_fig.y_range.start = 0
    perc_fig.x_range.range_padding = 0.1
    perc_fig.xgrid.grid_line_color = None
    perc_fig.axis.minor_tick_line_color = None
    perc_fig.outline_line_color = None
    perc_fig.legend.location = "top_left"
    perc_fig.legend.orientation = "horizontal"
    perc_fig.xaxis.major_label_orientation = 1

    perc_panel = Panel(child=perc_fig, title="Percentage of reads")

    tabs = Tabs(tabs=[nr_panel, perc_panel])

    output_file(outfile)
    save(tabs)

    return None


def main():
    """
    The main functionality of this script:
    reads files, merges them together in a table,
    save them as csv files and create graphs.
    
    Input:
      files as defined at the top (lines 26-30)
    Output:
      defined at the top (lines 34-37)
    """
    # 1. Parse and show arguments
    arguments = parse_arguments()

    message = (
        "\n"
        "These are the arguments you have provided:\n"
        "  INPUT:\n"
        "fastqc = {0},\n"
        "trimmomatic = {1}\n"
        "hugo = {2}\n"
        "classified = {3}\n"
        "unclassified: {4}\n"
        "mapped_reads: {5}\n"
        "  OUTPUT:\n"
        "counts = {6}\n"
        "percentages = {7}\n"
        "graph = {8}\n"
        "  OPTIONAL:\n"
        "log = {9}\n"
        "cores = {10}\n"
        "colours = {11}".format(
            arguments.fastqc,
            arguments.trimmomatic,
            arguments.hugo,
            arguments.classified,
            arguments.unclassified,
            arguments.mapped_reads,
            arguments.counts,
            arguments.percentages,
            arguments.graph,
            arguments.log,
            arguments.cores,
            arguments.colours,
        )
    )

    if arguments.log:
        if os.path.exists(arguments.log) and os.path.getsize(arguments.log) > 0:
            # File already exists and has been written to
            header = "\n---\n[{0}] Please find below the log of the current execution.".format(
                datetime.datetime.now()
            )
        else:
            header = "[{0}] Please find below the log of the current execution.".format(
                datetime.datetime.now()
            )
        with open(arguments.log, "a") as logfile:
            logfile.write("{0}\n---\n".format(header))
            logfile.write("{0}\n".format(message))
    else:
        print(message)

    threads = arguments.cores

    if arguments.log:
        with open(arguments.log, "a") as logfile:
            try:
                threads = snakemake.threads
                logfile.write(
                    "Variable 'threads' overwritten by snakemake = %i\n"
                    % snakemake.threads
                )
            except:
                logfile.write("Cores not provided by snakemake,\n")
                logfile.write("trying to read them from the -cpu argument...\n")

            if not isinstance(threads, int) or not threads in range(1, 32):
                threads = 4
            else:
                pass

            logfile.write(
                "Continuing with %i threads for reading fastq files.\n" % threads
            )
    else:
        try:
            threads = snakemake.threads
            print(
                "\nVariable 'threads' overwritten by snakemake = %i" % snakemake.threads
            )
        except:
            print("\nCores not provided by snakemake,")
            print("trying to read them from the -cpu argument...")

        if not isinstance(threads, int) or not threads in range(1, 32):
            threads = 4
        else:
            pass

        print("Continuing with %i threads for reading fastq files." % threads)

    # 2. Total (raw) read numbers/sample
    read_nrs = pd.read_csv(arguments.fastqc, delimiter="\t")[
        ["Sample", "Total Sequences"]
    ]
    # Note1: FastQC reports numbers of reads of pre- and post-QC fastq files,
    # I only want to have the pre-QC numbers now:
    pattern = re.compile(
        r".*(_|\.)R?1([_.].*$|$)"
    )  # A regex meaning: get anything up to "_R?1" and then either "_R?1[_.].*$" or "_R?1$".
    pre_qc = [col for col in read_nrs.Sample if bool(re.search(pattern, col))]
    read_nrs = read_nrs[read_nrs.Sample.isin(pre_qc)]
    read_nrs.Sample = read_nrs.Sample.str.replace(
        r"(_|\.)R?1([_.].*$|$)", ""
    )  # remove the "_R1" or "_1" suffix and anything that may be after it. N.B. a simple `.*`` doesn't work because otherwise it will greedily remove everything after the first occurrence of "_R1" e.g. "Sample_bla_bla_R138_ACGGT_R1" becomes "Sample_bla_bla"

    # Note2: I now only have the number of forward reads. To add reverse,
    # multiply this number by 2:
    read_nrs["Total Sequences"] = (read_nrs["Total Sequences"] * 2).astype(int)
    read_nrs.reset_index(drop=True, inplace=True)

    # Debug:
    # print(read_nrs.head())

    # 3. low-quality reads/sample (by Trimmomatic):
    lowq_nrs = pd.read_csv(arguments.trimmomatic, delimiter="\t")[
        ["Sample", "forward_only_surviving", "reverse_only_surviving", "dropped"]
    ]
    lowq_nrs.Sample = lowq_nrs.Sample.str.replace(
        r"(_|\.)R?1([_.].*$|$)", ""
    )  # remove the "_R1" or "_1" suffix and anything that may be after it. N.B. a simple `.*`` doesn't work because otherwise it will greedily remove everything after the first occurrence of "_R1" e.g. "Sample_bla_bla_R138_ACGGT_R1" becomes "Sample_bla_bla"
    # Note: trimmomatic drops read pairs ("dropped") and for some pairs
    # only one mate is of sufficiently high quality ("forward/reverse
    # only surviving"). The number of low-quality reads is calculated
    # by multiplying the dropped pairs by 2 and adding the forward
    # and reverse only mates:
    lowq_nrs["dropped_reads"] = (
        lowq_nrs["dropped"] * 2
        + lowq_nrs["forward_only_surviving"]
        + lowq_nrs["reverse_only_surviving"]
    ).astype(int)
    lowq_nrs.reset_index(drop=True, inplace=True)

    # Debug:
    # print(lowq_nrs.head())

    # 4. human reads/sample (Bowtie2 to human genome)
    human_nrs = count_non_human_reads(file_list=arguments.hugo, threads=threads)
    # Note: MultiQC does not have complete data for the Bowtie2
    # alignments. Instead I count the number of lines in each
    # of the fastq files that proceed to assembly. (Forward,
    # reverse and unpaired reads that did not map to the human
    # genome). This process may take a while!

    # Sort the dataframe by Sample ID to make it match the other dataframes
    human_nrs.sort_values(by=["Sample"], inplace=True)
    human_nrs.reset_index(inplace=True, drop=True)

    # Debug:
    # print(human_nrs.head())

    # 5. Reads classified by mapping to scaffolds/sample
    ## (Minimap reads to scaffolds > 500 nt long from SPAdes,
    ## that have been classified with Megablast to BLAST nt database
    ## and the LCA algorithm from Krona (see PZN methods).)
    classified_nrs = sum_superkingdoms(arguments.classified, arguments.mapped_reads)

    # Debug:
    # print(classified_nrs.head())

    # 6. Reads mapped to unclassified scaffolds/sample:
    ## (Minimap reads to scaffolds that Megablast could not
    ## assign.)
    unclassified_nrs = sum_unclassified(arguments.unclassified, arguments.mapped_reads)

    # Debug:
    # print(unclassified_nrs.head())

    # 7. Merge all these data into one dataframe:
    dfs = [read_nrs, lowq_nrs, human_nrs, classified_nrs, unclassified_nrs]

    nrs_df = pd.concat(dfs, axis=1, join="outer")

    # Calculate the human reads by subtracting low-quality
    # and non-human reads from the total number of reads:
    nrs_df["Human"] = nrs_df["Total Sequences"] - (
        nrs_df["dropped_reads"] + nrs_df["Non-human_reads"]
    )
    # Calculate the unclassified reads by summing the numbers
    # from the 'unclassified contigs table' and those from
    # the 'classified contigs table' that wern not assigned
    # a superkingdom:
    nrs_df["Unclassified"] = nrs_df.mapped_reads + nrs_df["not classified"]

    # And keep only the relevant columns,
    # in a logical order:
    nrs_df = nrs_df[
        [
            "Sample",
            "Total Sequences",
            "dropped_reads",
            "Human",
            "Archaea",
            "Bacteria",
            "Eukaryota",
            "Viruses",
            "Unclassified",
        ]
    ]
    nrs_df = nrs_df.iloc[
        :, ~nrs_df.columns.duplicated()
    ]  # remove duplicate columns, thanks https://stackoverflow.com/a/35798262 !
    #  and with proper names:
    nrs_df.rename(
        columns={"Total Sequences": "Total_reads", "dropped_reads": "Low-quality"},
        inplace=True,
    )
    # Replace NaN by 0:
    nrs_df.fillna(value=0, inplace=True)

    # And calculate the number of remaining reads
    # (i.e. high-quality reads that were not assembled
    # into scaffolds >= 500 nt and could thus not be
    # classified).
    nrs_df["Remaining"] = nrs_df["Total_reads"] - (
        nrs_df["Low-quality"]
        + nrs_df["Human"]
        + nrs_df["Archaea"]
        + nrs_df["Bacteria"]
        + nrs_df["Eukaryota"]
        + nrs_df["Viruses"]
        + nrs_df["Unclassified"]
    )

    # Integers with decimals are silly, remove them:
    nrs_df = nrs_df.apply(pd.to_numeric, downcast="integer", errors="ignore")

    # Debug:
    # print(nrs_df.head(20))

    validate_numbers(df=nrs_df, log=arguments.log)

    # 8. Write the numbers dataframe to a file:
    nrs_df.to_csv(arguments.counts, index=False)
    if arguments.log:
        with open(arguments.log, "a") as logfile:
            logfile.write("---\n\nWrote read-based table to:\t %s\n" % arguments.counts)
    else:
        print("Wrote read-based table to:\t %s" % arguments.counts)

    # 9. Create a dataframe with percentages:
    perc_df = pd.DataFrame()

    for header in nrs_df.columns:
        if header == "Sample":
            # Copy the Sample column
            perc_df[header] = nrs_df[header]
        elif header == "Total_reads":
            # Skip the Total_reads column (no need with percentages)
            pass
        else:
            # Create a percentage column for all other values
            perc_df[header] = nrs_df[header] / nrs_df["Total_reads"] * 100

    perc_df.to_csv(arguments.percentages, index=False)

    if arguments.log:
        with open(arguments.log, "a") as logfile:
            logfile.write("Wrote percentages table to:\t %s\n" % arguments.percentages)
    else:
        print("Wrote percentages table to:\t %s" % arguments.percentages)

    # 10. Create a stacked bar chart to visualise annotated reads per sample:
    draw_stacked_bars(
        df=nrs_df,
        perc=perc_df,
        sample="Sample",
        parts=[
            "Low-quality",
            "Human",
            "Archaea",
            "Bacteria",
            "Eukaryota",
            "Viruses",
            "Unclassified",
            "Remaining",
        ],
        outfile=arguments.graph,
    )

    if arguments.log:
        with open(arguments.log, "a") as logfile:
            logfile.write(
                "Created a read-based stacked bar chart in:\t %s\n" % arguments.graph
            )
    else:
        print("Created a read-based stacked bar chart in:\t %s" % arguments.graph)

    return None


### Run the script -----------------------------------------------

if __name__ == "__main__":
    main()
