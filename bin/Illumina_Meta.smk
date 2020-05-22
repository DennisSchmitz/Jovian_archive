"""
Authors:
    Dennis Schmitz, Sam Nooij, Florian Zwagemaker, Robert Verhagen,
    Jeroen Cremer, Thierry Janssens, Mark Kroon, Erwin van Wieringen,
    Annelies Kroneman, Harry Vennema, Marion Koopmans
Organization:
    Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
    Dutch Public Health institute (https://www.rivm.nl/en)
Department:
    Virology - Emerging and Endemic Viruses (EEV)
    Virology - Viruses of the Vaccination Program (VVP)
Date and license:
    23-08-2018, AGPL3 license
Homepage containing documentation, examples and a changelog:
    https://github.com/DennisSchmitz/Jovian
Funding:
    This project/research has received funding from the European Union’s
    Horizon 2020 research and innovation programme under grant agreement
    No. 643476. and the Dutch working group on molecular diagnostics (WMDI).
Automation:
    iRODS automatically executes this workflow for all "vir-ngs" labelled
    Illumina sequencing runs. See the internal Gitlab repo for the wrapper
    with additional information.
"""

#@################################################################################
#@#### Import config file, sample_sheet and set output folder names          #####
#@################################################################################

shell.executable("/bin/bash")

import pprint
import os
import yaml
yaml.warnings({'YAMLLoadWarning': False}) # Suppress yaml "unsafe" warnings.

# set dirs
logdir      =   "logs/"
cnf         =   "config/"
bench       =   "benchmark/"
conda_envs  =   "../envs/"
rls         =   "rules/"

configfile: f"{cnf}pipeline_parameters.yaml"
configfile: f"{cnf}variables.yaml"

datadir     =   "data/"
bindir      =   "bin/"
scrdir      =   "scripts/"
qc_pre      =   "FastQC_pretrim/"
qc_post     =   "FastQC_posttrim/"
cln         =   "cleaned_fastq/"
hugo_no_rm  =   "fastq_without_HuGo_removal/"
scf_raw     =   "scaffolds_raw/"
scf_filt    =   "scaffolds_filtered/"
taxclas     =   "taxonomic_classification/"
tbl         =   "tables/"
res         =   "results/"
hmap        =   "heatmaps/"
cnt         =   "counts/"
mqc_data    =   "multiqc_data/"
html        =   "html/"
fls         =   "files/"
fai_size    =   config["Illumina_meta"]["minlen"]   

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) # SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"

#@################################################################################
#@#### Jovian processes                                                      #####
#@################################################################################

localrules: 
    all,
    quantify_output,
    Concat_files,
    Concat_filtered_SNPs,
    HTML_IGVJs_generate_final,
    HTML_IGVJs_variable_parts


rule all:
    input:
        expand( "{p}{sample}_{read}.fq",
                p       =   f"{datadir + cln}",
                sample  =   SAMPLES,
                read    =   [   'pR1',
                                'pR2',
                                'unpaired'
                                ]
                ), # Extract unmapped & paired reads AND unpaired from HuGo alignment; i.e. cleaned fastqs
        expand( "{p}{sample}/scaffolds.fasta",
                p       =   f"{datadir + scf_filt}",
                sample  =   SAMPLES
                ), # SPAdes assembly output
        expand( "{p}{sample}_scaffolds_ge{l}nt.{extension}",
                p           =   f"{datadir + scf_filt}",
                sample      =   SAMPLES,
                l           =   config["Illumina_meta"]["minlen"],
                extension   =   [   'fasta',
                                    'fasta.fai'
                                    ]
                ), # Filtered SPAdes Scaffolds
        expand( "{p}{sample}_sorted.bam",
                p       =   f"{datadir + scf_filt}",
                sample  =   SAMPLES
                ), # BWA mem alignment for fragment length analysis
        expand( "{p}{sample}_{extension}",
                p           =   f"{datadir + scf_filt}",
                sample      =   SAMPLES,
                extension   =   [   'ORF_AA.fa',
                                    'ORF_NT.fa',
                                    'annotation.gff',
                                    'annotation.gff.gz',
                                    'annotation.gff.gz.tbi',
                                    'contig_ORF_count_list.txt'
                                    ]
                ), # Prodigal ORF prediction output
        expand( "{p}{sample}_{extension}",
                p           =   f"{datadir + scf_filt}",
                sample      =   SAMPLES,
                extension   =   [   'unfiltered.vcf',
                                    'filtered.vcf',
                                    'filtered.vcf.gz',
                                    'filtered.vcf.gz.tbi'
                                    ]
                ), # SNP calling output
        expand( "{p}{sample}_GC.bedgraph",
                p       =   f"{datadir + scf_filt}",
                sample  =   SAMPLES
                ), # Percentage GC content per specified window
        expand( "{p}{sample}.blastn",
                p       =   f"{datadir + taxclas}",
                sample  =   SAMPLES
                ), # MegablastN output
        expand( "{p}{sample}_{extension}",
                p           =   f"{datadir + tbl}",
                sample      =   SAMPLES,
                extension   =   [   'taxClassified.tsv',
                                    'taxUnclassified.tsv',
                                    'virusHost.tsv'
                                    ]
                ), # Tab seperated tables with merged data
        expand( "{p}{file}",
                p       =   f"{res}",
                file    =   [   'all_taxClassified.tsv',
                                'all_taxUnclassified.tsv',
                                'all_virusHost.tsv',
                                'all_filtered_SNPs.tsv'
                                ]
                ), # Concatenated classification, virus host and typing tool tables
        expand( "{p}{file}",
                p       =   f"{res}",
                file    =   [   'heatmaps/Superkingdoms_heatmap.html',
                                'Sample_composition_graph.html',
                                'Taxonomic_rank_statistics.tsv',
                                'Virus_rank_statistics.tsv',
                                'Phage_rank_statistics.tsv',
                                'Bacteria_rank_statistics.tsv'
                                ]
                ), # Taxonomic profile and heatmap output
        expand( "{p}{tax}_heatmap.html",
                p   =   f"{res + hmap}",
                tax =   [   'Virus',
                            'Phage',
                            'Bacteria'
                            ]
        ),
        expand( "{p}{file}.html",
                p       =   f"{res}",
                file    =   [   'multiqc',
                                'krona'
                                ]
                ), # HTML Reports
        expand( "{p}igv.html",
                p   =   f"{res}"
                ) # IGVjs index


#>############################################################################
#>#### Data quality control and cleaning                                 #####
#>############################################################################

include: f"{rls}QC_raw.smk"
include: f"{rls}CleanData.smk"
include: f"{rls}QC_clean.smk"

#>############################################################################
#>#### Removal of background host data                                   #####
#>############################################################################

include: f"{rls}BG_removal_1.smk"
include: f"{rls}BG_removal_2.smk"
include: f"{rls}BG_removal_3.smk"

#>############################################################################
#>#### De novo assembly and filtering                                    #####
#>############################################################################

include: f"{rls}assembly.smk"

#>############################################################################
#>#### Scaffold analysis and metrics                                     #####
#>############################################################################

include: f"{rls}Read2scaffold_alignment_with_rmDup_and_fraglength.smk"

include: f"{rls}SNP_calling.smk"
include: f"{rls}ORF_analysis.smk"
include: f"{rls}Contig_metrics.smk"
include: f"{rls}GC_content.smk"

#>############################################################################
#>#### Generate IGVjs HTML                                               #####
#>############################################################################

include: f"{rls}IGVjs.smk"

#>############################################################################
#>#### MultiQC report of pipeline metrics                                #####
#>############################################################################

include: f"{rls}MultiQC.smk"

#>############################################################################
#>#### Taxonomic classification & LCA                                    #####
#>############################################################################

include: f"{rls}Blast.smk"

if config["Illumina_meta"]["LCA"]["Krona"] == True:
    include: f"{rls}Krona_LCA.smk"

if config["Illumina_meta"]["LCA"]["mgkit"] == True:
    include: f"{rls}mgkit_LCA.smk"

include: f"{rls}KronaChart.smk"

#>############################################################################
#>#### Count annotated reads and visualize as stacked bar charts         #####
#>############################################################################

include: f"{rls}Count_reads.smk"
include: f"{rls}Concat_reads.smk"

#>############################################################################
#>#### Data wrangling                                                    #####
#>############################################################################

include: f"{rls}Merge_metrics.smk"
include: f"{rls}Concat_files.smk"
include: f"{rls}Concat_SNPs.smk"
include: f"{rls}Quantify.smk"

#>############################################################################
#>#### Make heatmaps for superkingdoms and viruses                       #####
#>############################################################################

include: f"{rls}Heatmaps.smk"

#@################################################################################
#@#### The `onstart` checker codeblock                                       #####
#@################################################################################

onstart:
    try:
        print("Checking if all specified files are accessible...")
        for filename in [ config["databases"]["background_ref"],
                         config["databases"]["Krona_taxonomy"],
                         config["databases"]["virusHostDB"],
                         config["databases"]["NCBI_new_taxdump_rankedlineage"],
                         config["databases"]["NCBI_new_taxdump_host"] ]:
            if not os.path.exists(filename):
                raise FileNotFoundError(filename)
    except FileNotFoundError as e:
        print("This file is not available or accessible: %s" % e)
        sys.exit(1)
    else:
        print("\tAll specified files are present!")
    shell("""
        mkdir -p results
        echo -e "\nLogging pipeline settings..."

        echo -e "\tGenerating methodological hash (fingerprint)..."
        echo -e "This is the link to the code used for this analysis:\thttps://github.com/DennisSchmitz/Jovian/tree/$(git log -n 1 --pretty=format:"%H")" > results/log_git.txt
        echo -e "This code with unique fingerprint $(git log -n1 --pretty=format:"%H") was committed by $(git log -n1 --pretty=format:"%an <%ae>") at $(git log -n1 --pretty=format:"%ad")" >> results/log_git.txt

        echo -e "\tGenerating full software list of current Conda environment (\"Jovian_master\")..."
        conda list > results/log_conda.txt

        echo -e "\tGenerating used databases log..."
        echo -e "==> User-specified background reference (default: Homo Sapiens NCBI GRch38 NO DECOY genome): <==\n$(ls -lah $(grep "    background_ref:" config/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" > results/log_db.txt
        echo -e "\n==> Virus-Host Interaction Database: <==\n$(ls -lah $(grep "    virusHostDB:" config/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" >> results/log_db.txt
        echo -e "\n==> Krona Taxonomy Database: <==\n$(ls -lah $(grep "    Krona_taxonomy:" config/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" >> results/log_db.txt
        echo -e "\n==> NCBI new_taxdump Database: <==\n$(ls -lah $(grep "    NCBI_new_taxdump_rankedlineage:" config/pipeline_parameters.yaml | cut -f 2 -d ":") $(grep "    NCBI_new_taxdump_host:" config/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" >> results/log_db.txt
        echo -e "\n==> NCBI Databases as specified in .ncbirc: <==\n$(ls -lah $(grep "BLASTDB=" .ncbirc | cut -f 2 -d "=" | tr "::" " "))\n" >> results/log_db.txt
        
        echo -e "\tGenerating config file log..."
        rm -f results/log_config.txt
        for file in config/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> results/log_config.txt
            cat ${{file}} >> results/log_config.txt
            echo -e "\n\n" >> results/log_config.txt
        done
    """)

#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################

onerror:
    shell("""
        rm -f data/scaffolds_filtered/*.html
        rm -rf data/html/
        rm -f data/taxonomic_classification/*_lca_*
        rm -f data/taxonomic_classification/*_nolca_*
    """)


onsuccess:
    shell("""
        echo -e "\nStarting virus typing, this may take a while...\n"
        bash jovian -vt all
        echo -e "Virus typing finished."

        echo -e "\nCleaning up..."
        echo -e "\tRemoving empty folders..."
        find data -depth -type d -not \( -path data/scaffolds_raw -prune \) -empty -delete

        echo -e "\tRemoving temporary files..."
        if [ "{config[remove_temp]}" != "0" ]
        then
            rm -rf data/FastQC_pretrim/
            rm -rf data/FastQC_posttrim/
            rm -rf data/cleaned_fastq/fastq_without_HuGo_removal/
            rm -f data/scaffolds_filtered/*_insert_size_histogram.pdf
            rm -f data/scaffolds_filtered/*_insert_size_metrics.txt
            rm -f data/scaffolds_filtered/*_MinLenFiltSummary.stats
            rm -f data/scaffolds_filtered/*_perMinLenFiltScaffold.stats
            rm -f data/scaffolds_filtered/*nt.fasta.sizes
            rm -f data/scaffolds_filtered/*.windows
            rm -f data/taxonomic_classification/*.taxtab
            rm -f data/taxonomic_classification/*.taxMagtab
            rm -f data/taxonomic_classification/*_lca_*
            rm -f data/taxonomic_classification/*_nolca_*
            rm -rf data/html/
        else
            echo -e "\t\tYou chose to not remove temp files: the human genome alignment files are not removed."
        fi

        echo -e "\tCreating symlinks for the interactive genome viewer..."
        bin/scripts/set_symlink.sh

        echo -e "\tGenerating HTML index of log files..."
        tree -hD --dirsfirst -H "../logs" -L 2 -T "Logs overview" --noreport --charset utf-8 -P "*" -o results/logfiles_index.html logs/

        echo -e "\tGenerating Snakemake report..."
        snakemake -s bin/Snakefile --unlock --config sample_sheet=sample_sheet.yaml
        snakemake -s bin/Snakefile --report results/snakemake_report.html --config sample_sheet=sample_sheet.yaml

        echo -e "Finished"
    """)
    
#? perORFcoverage output file van de bbtools scaffold metrics nog importeren in data wrangling part!

#@################################################################################
#@#### Specify Jovian's final output:                                        #####
#@################################################################################