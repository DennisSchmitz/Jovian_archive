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
    Bacteriology - Bacterial and Parasitic Diagnostics (BPD)
Date and license:
    19-12-2019, AGPL3 license
Description:
    Originally intended for the EAV internal control of the ENNGS ringtrail,
    later adapted for the nCoV outbreak. It's a WORK IN PROGRESS. It was
    initially intended for a coarse estimation of the BoC at different
    coverage thresholds.
Homepage containing documentation, examples and a changelog:
    https://github.com/DennisSchmitz/Jovian
Funding:
    This project/research has received funding from the European Union’s
    Horizon 2020 research and innovation programme under grant agreement
    No. 643476. and the Dutch working group on molecular diagnostics (WMDI).
"""


"""
#TODO:
#! Bug? the GC-contents track on IGVjs don't load. Same in the core workflow... weird.
- realign against new consensus? or too slow?
- Maybe MM2 is a better/faster aligner?
- Make a multiple alignment of all contigs/scaffolds versus the ref-alignment method as a double check!
   - the scaffolds then also need to be put in the proper orientation for MA to work.
- Simplify code with https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-multiext-function
- Add Ti/Tv ratio to BoC rule? Transitions versus transversions?
- See tips for improving variant filtering here: https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
- See tips for improving variant filtering in mail Jeroen EMC, 20200323
- Remove unneeded temp/intermediate files
    - Via onSucces clause
- Voeg de nieuwe consensus ook toe aan de IGVjs overview! --> kan niet
#TODO Onderstaande is ook handig voor Jovian zelf
- Verwijzen naar de output van een andere rule: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-dependencies
- For JupyterNotebook integration: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration
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
yaml.warnings({'YAMLLoadWarning': False}) # Suppress yaml "unsafe" warnings.

# Import sample sheet
SAMPLES = {}
with open(config["sample_sheet_reference_alignment"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) # SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"

# The reference file is given as a snakemake CLI argument from within the wrapper, so NOT via the pipeline_parameters.yaml
REFERENCE = config["reference"]
REFERENCE_BASENAME = os.path.splitext(os.path.basename(REFERENCE))[0]    # source: https://stackoverflow.com/questions/678236/how-to-get-the-filename-without-the-extension-from-a-path-in-python

# Set input directory, this is dependent on the Jovian output dir
INPUT_DIR_FILT_READS = config["reference_alignment"]["input_dir"]

# Set dir with conda-env files
CONDA_ENVS_DIR = "envs/"

# Set output base dir and sub-folder names, useful for easily changing the output locations during development.
OUTPUT_BASE_DIR = config["reference_alignment"]["output_dir"]
OUTPUT_DIR_REFERENCE = OUTPUT_BASE_DIR + "reference/"
OUTPUT_DIR_ALIGNMENT = OUTPUT_BASE_DIR + "alignment/"
OUTPUT_DIR_CONSENSUS_RAW = OUTPUT_BASE_DIR + "consensus_seqs/raw/"
OUTPUT_DIR_CONSENSUS_FILT = OUTPUT_BASE_DIR + "consensus_seqs/"
OUTPUT_DIR_BOC_ANALYSIS = OUTPUT_BASE_DIR + "BoC_analysis/"
OUTPUT_DIR_IGVjs = OUTPUT_BASE_DIR + "html/"

# Set output dir of results
OUTPUT_DIR_RESULTS = OUTPUT_BASE_DIR + "results/"
OUTPUT_IGVjs_HTML = OUTPUT_DIR_RESULTS + "igvjs.html"

# Set output dir of logfiles
OUTPUT_DIR_LOGS = config["reference_alignment"]["log_dir"] # NB the DRMAA logs will go the same dir as Jovian-core since this is set in the config.yaml file.
OUTPUT_DIR_BENCHMARKS = OUTPUT_DIR_LOGS + "benchmark/"

#@################################################################################
#@#### Specify Jovian's final output:                                        #####
#@################################################################################


localrules: 
    all,
    RA_index_reference,
    RA_concat_BoC_metrics,
    RA_HTML_IGVJs_variable_parts,
    RA_HTML_IGVJs_generate_final

rule all:
    input:
        expand("{out}{ref_basename}{extension}", out = OUTPUT_DIR_REFERENCE, ref_basename = REFERENCE_BASENAME, extension = [ '.fasta', '.fasta.1.bt2', '.fasta.fai', '.fasta.sizes', '.windows', '_GC.bedgraph' ]), # Copy of the reference file (for standardization and easy logging), bowtie2-indices (I've only specified one, but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated) and the GC-content files.
        expand("{out}{sample}_sorted.{extension}", out = OUTPUT_DIR_ALIGNMENT, sample = SAMPLES, extension = [ 'bam', 'bam.bai' ]), # The reference alignment (bam format) files.
        expand("{out}{sample}_{extension}", out = OUTPUT_DIR_CONSENSUS_RAW, sample = SAMPLES, extension = [ 'calls.vcf.gz', 'raw_consensus.fa' ]), # A zipped vcf file contained SNPs versus the given reference and a RAW consensus sequence, see explanation below for the meaning of RAW.
        expand("{out}{sample}.bedgraph", out = OUTPUT_DIR_CONSENSUS_FILT, sample = SAMPLES), # Lists the coverage of the alignment against the reference in a bedgraph format, is used to determine the coverage mask files below.
        expand("{out}{sample}_{filt_character}-filt_cov_ge_{thresholds}.fa", out = OUTPUT_DIR_CONSENSUS_FILT, sample = SAMPLES, filt_character = [ 'N', 'minus' ], thresholds = [ '1', '5', '10', '30', '100' ]), # Consensus sequences filtered for different coverage thresholds (1, 5, 10, 30 and 100). For each threshold two files are generated, one where failed positioned are replaced with a N nucleotide and the other where its replaced with a minus character (gap).
        expand("{out}{sample}_BoC{extension}", out = OUTPUT_DIR_BOC_ANALYSIS, sample = SAMPLES, extension = [ '_int.tsv', '_pct.tsv' ] ), # Output of the BoC analysis #TODO can probably removed after the concat rule is added.
        OUTPUT_DIR_RESULTS + "BoC_integer.tsv", # Integer BoC overview in .tsv format
        OUTPUT_DIR_RESULTS + "BoC_percentage.tsv", # Percentage BoC overview in .tsv format
        expand("{out}{ref_basename}_{extension}", out = OUTPUT_DIR_REFERENCE, ref_basename = REFERENCE_BASENAME , sample = SAMPLES, extension = [ 'ORF_AA.fa', 'ORF_NT.fa', 'annotation.gff', 'annotation.gff.gz', 'annotation.gff.gz.tbi' ]), # Prodigal ORF prediction output, required for the IGVjs visualisation
        OUTPUT_IGVjs_HTML, # IGVjs output html


#@################################################################################
#@#### Reference alignment extension processes                               #####
#@################################################################################


rule RA_index_reference:
    input:
        reference= REFERENCE
    output:
        reference_copy= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + ".fasta",
        reference_index= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + ".fasta.1.bt2", # I've only specified ".fasta.1.bt2", but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated. #TODO find a way to specify all output correctly (multiext snakemake syntax?)
    conda:
        CONDA_ENVS_DIR + "RA_ref_alignment.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_index_reference.txt"
    threads: 4
    log:
        OUTPUT_DIR_LOGS + "RA_index_reference.log"
    shell: # The reference is copied to the hardcoded subdir to make it standardized and easily logged.
        """
cp {input.reference} {output.reference_copy}
bowtie2-build --threads {threads} {output.reference_copy} {output.reference_copy} >> {log} 2>&1
        """


##########################!
# Nuttig voor IGVjs vis. Gejat uit Jovian core met minor changes, kunnen we waarschijlijk efficienter doen.
rule RA_reference_ORF_analysis:
    input:
        reference= rules.RA_index_reference.output.reference_copy
    output: 
        ORF_AA_fasta= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + "_ORF_AA.fa",
        ORF_NT_fasta= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + "_ORF_NT.fa",
        ORF_annotation_gff= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + "_annotation.gff",
        zipped_gff3= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + "_annotation.gff.gz",
        index_zipped_gff3= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + "_annotation.gff.gz.tbi",
    conda:
        CONDA_ENVS_DIR + "scaffold_analyses.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_reference_ORF_analysis.txt"
    log:
        OUTPUT_DIR_LOGS + "RA_reference_ORF_analysis.log"
    threads: 1
    params: #? Currently it's using the same prodigal settings as the main workflow, I see no problems with it since it's both foremost intended for viruses.
        procedure=config["ORF_prediction"]["procedure"],
        output_format=config["ORF_prediction"]["output_format"]
    shell:
        """
prodigal -q -i {input.reference} \
-a {output.ORF_AA_fasta} \
-d {output.ORF_NT_fasta} \
-o {output.ORF_annotation_gff} \
-p {params.procedure} \
-f {params.output_format} > {log} 2>&1
bgzip -c {output.ORF_annotation_gff} 2>> {log} 1> {output.zipped_gff3}
tabix -p gff {output.zipped_gff3} >> {log} 2>&1
        """


############################!
# Nuttig voor IGVjs vis. Gejat uit Jovian core met minor changes, kunnen we waarschijnlijk efficienter doen.
rule RA_determine_GC_content:
    input:
        fasta= rules.RA_index_reference.output.reference_copy,
    output:
        fasta_fai= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + ".fasta.fai",
        fasta_sizes= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + ".fasta.sizes",
        bed_windows= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + ".windows",
        GC_bed= OUTPUT_DIR_REFERENCE + REFERENCE_BASENAME + "_GC.bedgraph",
    conda:
        CONDA_ENVS_DIR + "scaffold_analyses.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_determine_GC_content.txt"
    log:
        OUTPUT_DIR_LOGS + "RA_determine_GC_content.log"
    threads: 1
    params:
        window_size="50"
    shell:
        """
samtools faidx -o {output.fasta_fai} {input.fasta} > {log} 2>&1
cut -f 1,2 {output.fasta_fai} 2> {log} 1> {output.fasta_sizes}
bedtools makewindows \
-g {output.fasta_sizes} \
-w {params.window_size} 2>> {log} 1> {output.bed_windows}
bedtools nuc \
-fi {input.fasta} \
-bed {output.bed_windows} 2>> {log} |\
cut -f 1-3,5 2>> {log} 1> {output.GC_bed}
        """


rule RA_align_to_reference:
    input:
        pR1= INPUT_DIR_FILT_READS + "{sample}_pR1.fq",
        pR2= INPUT_DIR_FILT_READS + "{sample}_pR2.fq",
        unpaired= INPUT_DIR_FILT_READS + "{sample}_unpaired.fq",
        reference= rules.RA_index_reference.output.reference_copy
    output:
        sorted_bam= OUTPUT_DIR_ALIGNMENT + "{sample}_sorted.bam",
        sorted_bam_index= OUTPUT_DIR_ALIGNMENT + "{sample}_sorted.bam.bai",
    conda:
        CONDA_ENVS_DIR + "RA_ref_alignment.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_align_to_reference_{sample}.txt"
    threads: config["threads"]["RA_align_to_reference"]
    log:
        OUTPUT_DIR_LOGS + "RA_align_to_reference_{sample}.log"
    params:
        alignment_type="--local",
    shell:
        """
bowtie2 --time --threads {threads} {params.alignment_type} \
-x {input.reference} \
-1 {input.pR1} \
-2 {input.pR2} \
-U {input.unpaired} 2> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.sorted_bam} >> {log} 2>&1
samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        """


##################################################################################################################################
# The rule below overwrites nucleotides in the ref based on the input reads and creates a RAW consensus sequence, RAW means:
#### A consensus genome is generated even if there is 0 coverage. At positions where no reads aligned it just takes
#### the reference nucleotide and inserts that into the consensus fasta. Obviously, this is not desirable and prone
#### to misinterpretation! Therefore, in a rule below, all positions with a coverage of 0 are masked, i.e. nucleotides
#### at that position are replaced with "N" or "-". Both replacements are made, since many aligners can handle gaps ("-")
#### but they cannot handle N's. However, the file with N's can be used to check for indels since this cannot be seen
#### in the file where every missing nucleotide is replaced with a "-".
##################################################################################################################################
#? This code cannot easily be made multithreaded for monopartite sequences, see https://github.com/samtools/bcftools/issues/949
#? If we ever get additional time later we can probably look into splitting it up, but now there is no time.
##################################################################################################################################
#? You will get deprecation warning saying 'samtools mpileup option `u` is functional, but deprecated' and that you should
#? use bcftools mpileup instead. I've tried this on 20200322 with 'bcftools mpileup -d 8000 -O u -f reference.fasta input.bam',
#? the -d 8000 is to be identical to the samtools default settings (bcftools default is 250) and -O u is to give it the same
#? output as the samtools output: there was a difference in SNP calling between the two (bcftools did not call one SNP that
#? samtools did). So I'm reluctant to change it. Leaving this here for an eventual later update.
#? Versions used: samtools 1.10 and bcftools 1.10 (also see https://github.com/samtools/bcftools/issues/852 )
##################################################################################################################################
rule RA_extract_raw_consensus:
    input:
        bam= rules.RA_align_to_reference.output.sorted_bam,
        reference= rules.RA_index_reference.output.reference_copy,
    output:
        gzipped_vcf= OUTPUT_DIR_CONSENSUS_RAW + "{sample}_calls.vcf.gz",
        raw_consensus_fasta= OUTPUT_DIR_CONSENSUS_RAW + "{sample}_raw_consensus.fa",
    conda:
        CONDA_ENVS_DIR + "RA_ref_alignment.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_extract_raw_consensus_{sample}.txt"
    threads: 1
    log:
        OUTPUT_DIR_LOGS + "RA_extract_raw_consensus_{sample}.log"
    params:
    shell: # Source: https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
        """
samtools mpileup -uf {input.reference} {input.bam} 2>> {log} |\
bcftools call --ploidy 1 -mv -O z -o {output.gzipped_vcf} >> {log} 2>&1
tabix {output.gzipped_vcf} >> {log} 2>&1
cat {input.reference} 2>> {log} |\
bcftools consensus {output.gzipped_vcf} | seqtk seq - > {output.raw_consensus_fasta} 2>> {log}
        """


#TODO kijken of dit multithreaded kan worden.
rule RA_extract_clean_consensus:
    input:
        bam= rules.RA_align_to_reference.output.sorted_bam,
        reference= rules.RA_index_reference.output.reference_copy,
        raw_consensus= rules.RA_extract_raw_consensus.output.raw_consensus_fasta, # Only needed for when there are no positions in the bed with a coverage of 0; in that case the RAW fasta is actually suitable for downstream processes and it is simply copied.
    output:
        bedgraph= OUTPUT_DIR_CONSENSUS_FILT + "{sample}.bedgraph",
        filt_consensus_N_filt_ge_1= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_N-filt_cov_ge_1.fa",
        filt_consensus_N_filt_ge_5= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_N-filt_cov_ge_5.fa",
        filt_consensus_N_filt_ge_10= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_N-filt_cov_ge_10.fa",
        filt_consensus_N_filt_ge_30= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_N-filt_cov_ge_30.fa",
        filt_consensus_N_filt_ge_100= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_N-filt_cov_ge_100.fa",
        filt_consensus_minus_filt_ge_1= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_minus-filt_cov_ge_1.fa",
        filt_consensus_minus_filt_ge_5= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_minus-filt_cov_ge_5.fa",
        filt_consensus_minus_filt_ge_10= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_minus-filt_cov_ge_10.fa",
        filt_consensus_minus_filt_ge_30= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_minus-filt_cov_ge_30.fa",
        filt_consensus_minus_filt_ge_100= OUTPUT_DIR_CONSENSUS_FILT + "{sample}_minus-filt_cov_ge_100.fa",
    conda:
        CONDA_ENVS_DIR + "RA_ref_alignment.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_extract_clean_consensus_{sample}.txt"
    threads: 1
    log:
        OUTPUT_DIR_LOGS + "RA_extract_clean_consensus_{sample}.log"
    params:
        output_folder= OUTPUT_DIR_CONSENSUS_FILT,
    shell:
        """
bash bin/scripts/RA_consensus_at_diff_coverages.sh {wildcards.sample} {input.bam} {input.reference} {input.raw_consensus} \
{params.output_folder} {log} >> {log} 2>&1
        """


#TODO make a python script or bash function/include to do this more efficiently, currently it's hacky, but it works
rule RA_determine_BoC_at_diff_cov_thresholds:
    input:
        bedgraph= rules.RA_extract_clean_consensus.output.bedgraph,
        reference= rules.RA_index_reference.output.reference_copy,
    output:
        percentage_BoC_tsv= OUTPUT_DIR_BOC_ANALYSIS + "{sample}_BoC_pct.tsv",
        integer_BoC_tsv= OUTPUT_DIR_BOC_ANALYSIS + "{sample}_BoC_int.tsv",
    conda:
        CONDA_ENVS_DIR + "RA_ref_alignment.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "determine_BoC_at_diff_cov_thresholds_{sample}.txt"
    threads: 1
    log:
        OUTPUT_DIR_LOGS + "RA_determine_BoC_at_diff_cov_thresholds_{sample}.log"
    params:
    shell:
        """
bash bin/scripts/RA_BoC_analysis.sh {wildcards.sample} {input.bedgraph} {input.reference} \
{output.percentage_BoC_tsv} {output.integer_BoC_tsv} >> {log} 2>&1
        """


rule RA_concat_BoC_metrics:
    input:
        BoC_int_tsv= expand("{out}{sample}_BoC_int.tsv", out = OUTPUT_DIR_BOC_ANALYSIS, sample = SAMPLES),
        BoC_pct_tsv= expand("{out}{sample}_BoC_pct.tsv", out = OUTPUT_DIR_BOC_ANALYSIS, sample = SAMPLES),
    output:
        combined_BoC_int_tsv= OUTPUT_DIR_RESULTS + "BoC_integer.tsv",
        combined_BoC_pct_tsv= OUTPUT_DIR_RESULTS + "BoC_percentage.tsv",
    conda:
        CONDA_ENVS_DIR + "RA_ref_alignment.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "concat_BoC_metrics.txt"
    threads: 1
    log:
        OUTPUT_DIR_LOGS + "RA_concat_BoC_metrics.log"
    params:
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_int_tsv}
cat {input.BoC_int_tsv} >> {output.combined_BoC_int_tsv}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_pct_tsv}
cat {input.BoC_pct_tsv} >> {output.combined_BoC_pct_tsv}
        """


#@################################################################################
#@#### Make IGVjs html                                                       #####
#@################################################################################
#############! All code below here should be integrated with Jovian core workflow


rule RA_HTML_IGVJs_variable_parts:
    input:
        fasta= rules.RA_index_reference.output.reference_copy,
        ref_GC_bedgraph= rules.RA_determine_GC_content.output.GC_bed,
        ref_zipped_ORF_gff= rules.RA_reference_ORF_analysis.output.zipped_gff3,
        basepath_zipped_SNP_vcf= rules.RA_extract_raw_consensus.output.gzipped_vcf,
        basepath_sorted_bam= rules.RA_align_to_reference.output.sorted_bam,
    output:
        tab_output= OUTPUT_DIR_IGVjs + "2_tab_{sample}",
        div_output= OUTPUT_DIR_IGVjs + "4_html_divs_{sample}",
        js_flex_output= OUTPUT_DIR_IGVjs + "6_js_flex_{sample}",
    conda:
        CONDA_ENVS_DIR + "data_wrangling.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    log:
        OUTPUT_DIR_LOGS + "RA_HTML_IGVJs_variable_parts_{sample}.log"
    params:
    shell:
        """
bash bin/html/RA_igvjs_write_tabs.sh {wildcards.sample} {output.tab_output}

bash bin/html/RA_igvjs_write_divs.sh {wildcards.sample} {output.div_output}

bash bin/html/RA_igvjs_write_flex_js_middle.sh {wildcards.sample} {output.js_flex_output} \
{input.fasta} {input.ref_GC_bedgraph} {input.ref_zipped_ORF_gff} \
{input.basepath_zipped_SNP_vcf} {input.basepath_sorted_bam}
        """


rule RA_HTML_IGVJs_generate_final:
    input:
        expand("{out}{chunk_name}_{sample}", out = OUTPUT_DIR_IGVjs, chunk_name = [ '2_tab', '4_html_divs', '6_js_flex' ], sample = SAMPLES)
    output:
        OUTPUT_IGVjs_HTML
    conda:
        CONDA_ENVS_DIR + "data_wrangling.yaml"
    benchmark:
        OUTPUT_DIR_BENCHMARKS + "RA_HTML_IGVJs_generate_final.txt"
    threads: 1
    log:
        OUTPUT_DIR_LOGS + "RA_HTML_IGVJs_generate_final.log"
    params:
        tab_basename= OUTPUT_DIR_IGVjs + "2_tab_",
        div_basename= OUTPUT_DIR_IGVjs + "4_html_divs_",
        js_flex_output= OUTPUT_DIR_IGVjs + "6_js_flex_",
    shell:
        """
cat files/html_chunks/1_header.html > {output}
cat {params.tab_basename}* >> {output}
cat files/html_chunks/3_tab_explanation_RA.html >> {output}
cat {params.div_basename}* >> {output}
cat files/html_chunks/5_js_begin.html >> {output}
cat {params.js_flex_output}* >> {output}
cat files/html_chunks/7_js_end.html >> {output}
        """


#@################################################################################
#@#### These are the conditional cleanup rules                               #####
#@################################################################################


onsuccess:
    shell("""
        echo -e "\nCleaning up..."
        #TODO hier nog zo'n remove temp checker inbouwen.
        rm -rf {OUTPUT_DIR_IGVjs}   # Remove intermediate IGVjs html chunks.

        echo -e "\tCreating symlinks for the interactive genome viewer..."
        bin/scripts/set_symlink.sh

        echo -e "Finished"
    """)