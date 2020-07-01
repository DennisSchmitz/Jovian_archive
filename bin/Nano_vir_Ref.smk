"""
    introduction text
"""

#@################################################################################
#@#### Import config file, sample_sheet and set output folder names          #####
#@################################################################################

shell.executable("/bin/bash")

# Load config files
configfile: "config/pipeline_parameters.yaml"
configfile: "config/variables.yaml"

# Load libraries

import pprint
import os
import yaml
yaml.warnings({'YAMLLoadWarning': False})

from globals import *

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file)


primerfile  =   config["primer_file"]


reference           =   config["reference_file"]
reference_basename  =   os.path.splitext(os.path.basename(reference))[0]


#@################################################################################
#@#### Jovian processes                                                      #####
#@################################################################################

rule all:
    input: 
        expand( "{p}{sample}.fasta",
                p       =   f"{datadir + cons + raw}", 
                sample  =   SAMPLES
                ),
        expand("{p}{sample}.bedgraph", 
                p       =   f"{datadir + cons + filt}",
                sample  =   SAMPLES
                ),
        expand( "{p}{sample}_{filt_character}-filt_cov_ge_{thresholds}.fa",
                p               =   f"{res + seqs}",
                sample          =   SAMPLES,
                filt_character  =   [   'N',
                                        'minus'
                                        ],
                thresholds      =   [   '1',
                                        '5', 
                                        '10',
                                        '30',
                                        '100'
                                        ]),
        expand( f"{res}BoC_int.tsv"
                ),
        expand( f"{res}BoC_pct.tsv"
                ),
        expand( f"{res}igv.html"
                ),
        expand( f"{res}multiqc.html"
                )

    #>############################################################################
    #>#### Data quality control and cleaning                                 #####
    #>############################################################################

include: f"{rls}Nano_Ref_index.smk"
include: f"{rls}Nano_Ref_ORF-analysis.smk"
include: f"{rls}Nano_Ref_GC-content.smk"

include: f"{rls}Nano_Ref_adp_trim.smk"

include: f"{rls}Nano_Ref_Cut-primers.smk"

include: f"{rls}Nano_Ref_Cleanup.smk"

include: f"{rls}Nano_Ref_HuGo_removal_pt1.smk"

include: f"{rls}Nano_Ref_HuGo_removal_pt2.smk"

include: f"{rls}Nano_Ref_alignment_pt1.smk"

include: f"{rls}Nano_Ref_alignment_pt2.smk"

include: f"{rls}Nano_Ref_raw-consensus.smk"

include: f"{rls}Nano_Ref_clean-consensus.smk"

include: f"{rls}Nano_Ref_calc-boc.smk"

include: f"{rls}Nano_Ref_concat-boc.smk"

include: f"{rls}Nano_Ref_IGVjs.smk"

include: f"{rls}Nano_Ref_MultiQC.smk"


onsuccess:
    shell("""
        echo -e "\tCreating symlinks for the interactive genome viewer..."
        {bindir}{scrdir}set_symlink.sh

        echo -e "\tGenerating HTML index of log files..."
        tree -hD --dirsfirst -H "../logs" -L 2 -T "Logs overview" --noreport --charset utf-8 -P "*" -o {res}logfiles_index.html {logdir}

        echo -e "\tGenerating Snakemake report..."
        snakemake -s {bindir}Nano_vir_Ref.smk --unlock --profile config
        snakemake -s {bindir}Nano_vir_Ref.smk --report {res}snakemake_report.html --profile config

        echo -e "Finished"
    """)


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