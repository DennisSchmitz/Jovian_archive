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
# _-- 14 Nov 2018 update: changed notebook into regular Python script that is run by the Snakefile_
# 
# Required Python packages:
#  - Pandas
#  - Bokeh

#IMPORT required libraries---------------------------------
import numpy as np
import pandas as pd
from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, ColumnDataSource

#Initialise global VARIABLES-------------------------------
classifications = snakemake.input['classified']
read_pairs = snakemake.input['numbers']

COLOUR = [ "#000000" ]
PHAGE_FAMILY_LIST = [ "Myoviridae", "Siphoviridae", "Podoviridae", "Lipothrixviridae", 
              "Rudiviridae", "Ampullaviridae", "Bicaudaviridae", "Clavaviridae", 
              "Corticoviridae", "Cystoviridae", "Fuselloviridae", "Globuloviridae", 
              "Guttaviridae", "Inoviridae", "Leviviridae", "Microviridae", 
              "Plasmaviridae", "Tectiviridae" ]

#Define FUNCTIONS------------------------------------------

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
    classifications_df = pd.read_csv(classifications, delimiter='\t')

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

def merge_numbers_and_classifications(nr_df, cl_df):
    """
    Merges two dataframes to add the total number of (raw) reads per sample
      to the results table with classified scaffolds.
    Also calculates the percentage of reads mapped to each scaffold of the
      total number of reads.
    """
    merged_df = cl_df.merge(nr_df, left_on="Sample_name", right_on="Sample")
    merged_df["Percentage"] = merged_df.reads / merged_df.read_pairs * 100
    
    return(merged_df)

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
        #Count how many taxa have been reported
        for t in [ "superkingdom", "phylum", "class", "order", "family", "genus", "species" ]:
            f.write("%s\t%i\n" % (t, classifications_df[t].nunique()))
    
    return(None)

def draw_superkingdoms_heatmap(df, outfile):
    """
    Input: Dataframe with scaffold classifications and quantifications
    Output: Heatmap of superkingdoms found in the analysed samples (as html file)
    """
    #Count the percentages of Archaea, Bacteria, Eukaryota and Viruses per sample
    superkingdom_sums = pd.DataFrame(df.groupby(
    [ "Sample_name", "superkingdom" ]).sum()
                         [[ "reads", "Percentage" ]])
    
    #Create a heatmap of superkingdoms, summed per sample
    superkingdom_sums.reset_index(inplace=True) #to use MultiIndex "Sample_name" and "superkingdom" as columns
    #And save it
    superkingdom_sums.to_csv("results/Superkingdoms_quantities_per_sample.csv", index=False)
    
    def super_heatmap(df):
        """
    (Only for superkingdoms per sample.)
        """
        samples = df["Sample_name"].astype(str)
        assigned = df["superkingdom"].astype(str)
        reads = df["reads"].astype(int)
        percent_of_total = df["Percentage"].astype(float)

        colors = len(reads) * COLOUR #multiply to make an equally long list

        max_load = max(percent_of_total)
        alphas = [ min( x / float(max_load), 0.9) + 0.1 for x in percent_of_total ]

        source = ColumnDataSource(
            data = dict(samples=samples, assigned=assigned,
                        reads=reads, percent_of_total=percent_of_total, 
                        colors=colors, alphas=alphas)
        )

        TOOLS = "hover, save, pan, box_zoom, wheel_zoom, reset"

        p = figure(title = "Superkingdoms heatmap",
                  #If desired, the sample can be displayed as "Run x, sample y"
                  # uncomment the next line if desired
                  #x_range = [ "Run %s, sample %s" % (x.split('_')[0], x.split('_')[1]) for x in list(sorted(set(samples))) ],
                  x_range = list(sorted(set(samples))),
                  y_range = list(reversed(sorted(set(assigned)))), #reverse to order 'from top to bottom'
                  x_axis_location = "above",
                  toolbar_location="right",
                  tools = TOOLS)

        #Edit the size of the heatmap when there are many samples
        if len(set(samples)) > 20:
            p.plot_width = int(len(set(samples)) * 25)
        else:
            pass
        
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        if len(set(assigned)) > 15:
            p.axis.major_label_text_font_size = "10pt"
        else:
            p.axis.major_label_text_font_size = "12pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi/4
        p.title.text_color = COLOUR[0]
        p.title.text_font_size = "16pt"
        p.title.align = 'right'

        p.rect("samples", "assigned", 1, 1, source=source,
               color="colors", alpha="alphas", line_color=None)

        p.select_one(HoverTool).tooltips = [
            ('Sample', "@samples"),
            ('Taxon' , "@assigned"),
            ('Number of reads', "@reads"),
            ('Percentage of total',  "@percent_of_total %")
        ]

        output_file(snakemake.output["super"], title="Superkingdoms heatmap")
        save(p)

        print("Heatmap \"Superkingdoms heatmap\" has been written to %s" % snakemake.output["super"])
        
        return(None)

    super_heatmap(superkingdom_sums)
    
    return(None)

def filter_viruses(df, outfile):
    """
    Input: Dataframe with scaffold classification results
    Output: Dataframe with only classifications in the superkingdom "Viruses"
      (filtered rows), and a list of statistics (as tabular txt file)
    """
    virus_df = df[df.superkingdom == "Viruses"]
    
    header="taxonomic_level\tnumber_found\n"
    with open(outfile, 'w') as f:
        f.write(header)
        for t in [ "superkingdom", "phylum", "class", "order", "family", "genus", "species" ]:
            f.write("%s\t%i\n" % (t, virus_df[t].nunique()))
    
    return(virus_df)

def discard_phages(df):
    """
    Input: Dataframe with scaffold classification results of viruses
    Output: Dataframe with the same information, but phage families removed
    """
    phageless_df = df[~df['family'].isin(PHAGE_FAMILY_LIST)]
    
    return(phageless_df)

def select_phages(df):
    """
    Input: Dataframe with scaffold classification results of viruses
    Output: Dataframe with only phage families
    """
    phage_df = df[df['family'].isin(PHAGE_FAMILY_LIST)]
    
    return(phage_df)

def generate_virus_rank_list():
    """
    Input: list of html files (heatmaps) to be created (from snakemake)
    Output: list of taxonomic ranks for which to make a heatmap
    """
    files_list = snakemake.output["virus"]
    
    rank_list = []
    
    for file in files_list:
        file = file.split('/')[-1]
        rank = file.split('_')[1]
        rank_list.append(rank)

    return(rank_list)

def create_heatmaps(df, rank_list):
    """
    Input: Dataframe with classifications and quantifications for different
      taxonomic ranks, and a list of ranks for which to generate heatmaps
    Output: Heatmaps for the given taxonomic ranks (as html files)
    """
    
    def create_heatmap(df, all_samples, title, filename, taxonomic_level="species"):
        samples = df["Sample_name"].astype(str)
        scaffolds = df["#ID"].astype(str)
        assigned = df["tax_name"].astype(str)
        taxonomy = df[taxonomic_level].astype(str)
        reads = df["reads"].astype(int)
        total_reads = df["read_pairs"].astype(int)
        percent_of_total = df["Percentage"].astype(float)
        coverage = df["Avg_fold"].astype(int)
        contig_length = df["Length"].astype(int)

        colors = len(reads) * COLOUR #multiply to make an equally long list
        
        max_load = max(percent_of_total)
        alphas = [ min( x / float(max_load), 0.9) + 0.1 for x in percent_of_total ]
        
        source = ColumnDataSource(
            data = dict(samples=samples, scaffolds=scaffolds,
                        assigned=assigned, taxonomy=taxonomy,
                        reads=reads, total_reads=total_reads,
                        percent_of_total=percent_of_total, 
                        coverage=coverage,
                        contig_length=contig_length,
                        colors=colors, alphas=alphas)
        )
        
        TOOLS = "hover, save, pan, box_zoom, wheel_zoom, reset"

        p = figure(title = title,
                  #If desired, the sample can be displayed as "Run x, sample y"
                  # uncomment the next line if desired
                  #x_range = [ "Run %s, sample %s" % (x.split('_')[0], x.split('_')[1]) for x in list(sorted(set(samples))) ],
                  x_range = all_samples,
                  y_range = list(reversed(sorted(set(taxonomy)))), #reverse to order 'from top to bottom'
                  x_axis_location = "above",
                  toolbar_location="right",
                  tools = TOOLS)

        #Edit the size of the heatmap when there are many samples and/or taxa
        if len(set(samples)) > 20:
            p.plot_width = int(len(set(samples)) * 25)
        else:
            pass
        if len(set(taxonomy)) > 100:
            p.plot_height = int(p.plot_height * 3)
            p.plot_width = int(p.plot_width * 1.5)
        elif len(set(taxonomy)) > 50:
            p.plot_height = int(p.plot_height * 2)
            p.plot_width = int(p.plot_width * 1.2)
        elif len(set(taxonomy)) > 25:
            p.plot_height = int(p.plot_height * 1.2)
        else:
            pass
        
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        if len(set(assigned)) > 15:
            p.axis.major_label_text_font_size = "10pt"
        else:
            p.axis.major_label_text_font_size = "12pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi/4
        p.title.text_color = COLOUR[0]
        p.title.text_font_size = "16pt"
        p.title.align = 'right'

        p.rect("samples", "taxonomy", 1, 1, source=source,
               color="colors", alpha="alphas", line_color=None)

        p.select_one(HoverTool).tooltips = [
            ('Sample', "@samples"),
            ('Scaffold', "@scaffolds"),
            ('Taxon' , "@assigned"),
            ('Number of reads', "@reads (@percent_of_total % of sample total)"),
            ('Scaffold length', "@contig_length"),
            ('Average Depth of Coverage', "@coverage")
        ]

        output_file(filename, title=title)
        save(p)
        print("The heatmap %s has been created and written to: %s" % (title, filename))
        return(None)
    
    print("Need to make: %s" % snakemake.output["virus"])
    
    for t in rank_list:
        try:
            tmp_df = df[df[t].notnull()]
            all_samples = list(sorted(set(df["Sample_name"])))
            create_heatmap(df = tmp_df, all_samples = all_samples, 
                           title = "%s heatmap" % t,
                           filename = "results/heatmaps/Virus_%s_heatmap.html" % t,
                           taxonomic_level = t)
        except ValueError:
            print("Taxonomic level %s has no information" % t)
    
    return(None)

def filter_bacteria(df, outfile):
    """
    Input: Dataframe with scaffold classification results
    Output: Dataframe with only classifications in the superkingdom "Bacteria"
      (filtered rows), and a list of statistics (as tabular txt file)
    """
    bact_df = df[df.superkingdom == "Bacteria"]
    
    header="taxonomic_level\tnumber_found\n"
    with open(outfile, 'w') as f:
        f.write(header)
        for t in [ "superkingdom", "phylum", "class", "order", "family", "genus", "species" ]:
            f.write("%s\t%i\n" % (t, bact_df[t].nunique()))
    
    return(bact_df)

def generate_bact_rank_list():
    """
    Input: list of html files (heatmaps) to be created (from snakemake)
    Output: list of taxonomic ranks for which to make a heatmap
    """
    files_list = snakemake.output["bact"]
    
    rank_list = []
    
    for file in files_list:
        file = file.split('/')[-1]
        rank = file.split('_')[1]
        rank_list.append(rank)

    return(rank_list)

def bact_heatmaps(df, rank_list):
    """
    Input: Dataframe with classifications and quantifications for different
      taxonomic ranks, and a list of ranks for which to generate heatmaps
    Output: Heatmaps for the given taxonomic ranks (as html files)
    """
    
    def create_heatmap(df, all_samples, title, filename, taxonomic_level="species"):
        samples = df["Sample_name"].astype(str)
        scaffolds = df["#ID"].astype(str)
        assigned = df["tax_name"].astype(str)
        taxonomy = df[taxonomic_level].astype(str)
        reads = df["reads"].astype(int)
        total_reads = df["read_pairs"].astype(int)
        percent_of_total = df["Percentage"].astype(float)
        coverage = df["Avg_fold"].astype(int)
        contig_length = df["Length"].astype(int)

        colors = len(reads) * COLOUR #multiply to make an equally long list
        
        max_load = max(percent_of_total)
        alphas = [ min( x / float(max_load), 0.9) + 0.1 for x in percent_of_total ]
        
        source = ColumnDataSource(
            data = dict(samples=samples, scaffolds=scaffolds,
                        assigned=assigned, taxonomy=taxonomy,
                        reads=reads, total_reads=total_reads,
                        percent_of_total=percent_of_total, 
                        coverage=coverage,
                        contig_length=contig_length,
                        colors=colors, alphas=alphas)
        )
        
        TOOLS = "hover, save, pan, box_zoom, wheel_zoom, reset"

        p = figure(title = title,
                  #If desired, the sample can be displayed as "Run x, sample y"
                  # uncomment the next line if desired
                  #x_range = [ "Run %s, sample %s" % (x.split('_')[0], x.split('_')[1]) for x in list(sorted(set(samples))) ],
                  x_range = all_samples,
                  y_range = list(reversed(sorted(set(taxonomy)))), #reverse to order 'from top to bottom'
                  x_axis_location = "above",
                  toolbar_location="right",
                  tools = TOOLS)

        #Edit the size of the heatmap when there are many samples and/or taxa
        if len(set(samples)) > 20:
            p.plot_width = int(len(set(samples)) * 25)
        else:
            pass
        if len(set(taxonomy)) > 100:
            p.plot_height = int(p.plot_height * 3)
            p.plot_width = int(p.plot_width * 1.5)
        elif len(set(taxonomy)) > 50:
            p.plot_height = int(p.plot_height * 2)
            p.plot_width = int(p.plot_width * 1.2)
        elif len(set(taxonomy)) > 25:
            p.plot_height = int(p.plot_height * 1.2)
        else:
            pass
        
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        if len(set(assigned)) > 15:
            p.axis.major_label_text_font_size = "10pt"
        else:
            p.axis.major_label_text_font_size = "12pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi/4
        p.title.text_color = COLOUR[0]
        p.title.text_font_size = "16pt"
        p.title.align = 'right'

        p.rect("samples", "taxonomy", 1, 1, source=source,
               color="colors", alpha="alphas", line_color=None)

        p.select_one(HoverTool).tooltips = [
            ('Sample', "@samples"),
            ('Scaffold', "@scaffolds"),
            ('Taxon' , "@assigned"),
            ('Number of reads', "@reads (@percent_of_total % of sample total)"),
            ('Scaffold length', "@contig_length"),
            ('Average Depth of Coverage', "@coverage")
        ]

        output_file(filename, title=title)
        save(p)
        print("The heatmap %s has been created and written to: %s" % (title, filename))
        return(None)
    
    print("Need to make: %s" % snakemake.output["bact"])
    
    for t in rank_list:
        try:
            tmp_df = df[df[t].notnull()]
            all_samples = list(sorted(set(df["Sample_name"])))
            create_heatmap(df = tmp_df, all_samples = all_samples, 
                           title = "%s heatmap" % t,
                           filename = "results/heatmaps/Bacteria_%s_heatmap.html" % t,
                           taxonomic_level = t)
        except ValueError:
            print("Taxonomic level %s has no information" % t)
    
    return(None)

def phage_heatmaps(df, rank_list):
    """
    Input: Dataframe with classifications and quantifications for different
      taxonomic ranks, and a list of ranks for which to generate heatmaps
    Output: Heatmaps for the given taxonomic ranks (as html files)
    """
    
    def create_heatmap(df, all_samples, title, filename, taxonomic_level="species"):
        samples = df["Sample_name"].astype(str)
        scaffolds = df["#ID"].astype(str)
        assigned = df["tax_name"].astype(str)
        taxonomy = df[taxonomic_level].astype(str)
        reads = df["reads"].astype(int)
        total_reads = df["read_pairs"].astype(int)
        percent_of_total = df["Percentage"].astype(float)
        coverage = df["Avg_fold"].astype(int)
        contig_length = df["Length"].astype(int)

        colors = len(reads) * COLOUR #multiply to make an equally long list
        
        max_load = max(percent_of_total)
        alphas = [ min( x / float(max_load), 0.9) + 0.1 for x in percent_of_total ]
        
        source = ColumnDataSource(
            data = dict(samples=samples, scaffolds=scaffolds,
                        assigned=assigned, taxonomy=taxonomy,
                        reads=reads, total_reads=total_reads,
                        percent_of_total=percent_of_total, 
                        coverage=coverage,
                        contig_length=contig_length,
                        colors=colors, alphas=alphas)
        )
        
        TOOLS = "hover, save, pan, box_zoom, wheel_zoom, reset"

        p = figure(title = title,
                  #If desired, the sample can be displayed as "Run x, sample y"
                  # uncomment the next line if desired
                  #x_range = [ "Run %s, sample %s" % (x.split('_')[0], x.split('_')[1]) for x in list(sorted(set(samples))) ],
                  x_range = all_samples,
                  y_range = list(reversed(sorted(set(taxonomy)))), #reverse to order 'from top to bottom'
                  x_axis_location = "above",
                  toolbar_location="right",
                  tools = TOOLS)

        #Edit the size of the heatmap when there are many samples and/or taxa
        if len(set(samples)) > 20:
            p.plot_width = int(len(set(samples)) * 25)
        else:
            pass
        if len(set(taxonomy)) > 100:
            p.plot_height = int(p.plot_height * 3)
            p.plot_width = int(p.plot_width * 1.5)
        elif len(set(taxonomy)) > 50:
            p.plot_height = int(p.plot_height * 2)
            p.plot_width = int(p.plot_width * 1.2)
        elif len(set(taxonomy)) > 25:
            p.plot_height = int(p.plot_height * 1.2)
        else:
            pass
        
        p.grid.grid_line_color = None
        p.axis.axis_line_color = None
        p.axis.major_tick_line_color = None
        if len(set(assigned)) > 15:
            p.axis.major_label_text_font_size = "10pt"
        else:
            p.axis.major_label_text_font_size = "12pt"
        p.axis.major_label_standoff = 0
        p.xaxis.major_label_orientation = np.pi/4
        p.title.text_color = COLOUR[0]
        p.title.text_font_size = "16pt"
        p.title.align = 'right'

        p.rect("samples", "taxonomy", 1, 1, source=source,
               color="colors", alpha="alphas", line_color=None)

        p.select_one(HoverTool).tooltips = [
            ('Sample', "@samples"),
            ('Scaffold', "@scaffolds"),
            ('Taxon' , "@assigned"),
            ('Number of reads', "@reads (@percent_of_total % of sample total)"),
            ('Scaffold length', "@contig_length"),
            ('Average Depth of Coverage', "@coverage")
        ]

        output_file(filename, title=title)
        save(p)
        print("The heatmap %s has been created and written to: %s" % (title, filename))
        return(None)
    
    print("Need to make: %s" % snakemake.output["phage"])
    
    for t in rank_list:
        try:
            tmp_df = df[df[t].notnull()]
            all_samples = list(sorted(set(df["Sample_name"])))
            create_heatmap(df = tmp_df, all_samples = all_samples, 
                           title = "%s heatmap" % t,
                           filename = "results/heatmaps/Phage_%s_heatmap.html" % t,
                           taxonomic_level = t)
        except ValueError:
            print("Taxonomic level %s has no information" % t)
    
    return(None)

#Script EXECUTION--------------------------------------
if __name__ == "__main__":
    #Create the required dataframe:
    numbers_df = read_numbers(read_pairs)
    classifications_df = read_classifications(classifications)
    
    merged_df = merge_numbers_and_classifications(numbers_df, classifications_df)
    
    #Report statistics: numbers of taxonomic ranks found:
    report_taxonomic_statistics(merged_df, snakemake.output["stats"])
    
    #Create a heatmap per superkingdom found in the samples:
    draw_superkingdoms_heatmap(merged_df, snakemake.output["super"])
    
    #Create a virus-only dataframe:
    virus_df = discard_phages(
        filter_viruses(merged_df, snakemake.output["vir_stats"])
    )
    
    #And create virus heatmaps:
    create_heatmaps(virus_df,
                    generate_virus_rank_list()
                   )
    
    #Create a phage-only dataframe (based on selected families):
    phage_df = select_phages(filter_viruses(
        merged_df, snakemake.output["phage_stats"])
                            )
    
    #And create phage heatmaps:
    phage_heatmaps(phage_df,
                    generate_virus_rank_list()
                   )
    
    #Create a bacteria-only dataframe:
    bact_df = filter_bacteria(merged_df, snakemake.output["bact_stats"])
    
    #And create bacteria heatmaps:
    bact_heatmaps(bact_df, generate_bact_rank_list())
