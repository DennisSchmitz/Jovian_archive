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


primerfile  =   config["primers"]


reference           =   config["reference"]
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
        expand( f"{res}IGVjs.html"
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
        bash bin/scripts/set_symlink.sh
    """)