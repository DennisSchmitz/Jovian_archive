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
- integrate with IGVjs so you can assess the read alignment
    - Add prodigal ORF prediction, or parse from reference
- Make the genomecov threshold a parameter instead of only 0 cov regions are replaced with N-nucleotides
- Maybe MM2 is a better/faster aligner?
- Make a multiple alignment of all contigs/scaffolds versus the ref-alignment method as a double check!
   - the scaffolds then also need to be put in the proper orientation for MA to work.
- Simplify code with https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#the-multiext-function
- Add Ti/Tv ratio to BoC rule? Transitions versus transversions?
- See tips for improving variant filtering here: https://github.com/samtools/bcftools/wiki/HOWTOs#consensus-calling
- Remove unneeded temp/intermediate files
    - Via onSucces clause
#TODO Onderstaande is ook handig voor Jovian zelf
- Verwijzen naar de output van een andere rule: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-dependencies
- For JupyterNotebook integration: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#jupyter-notebook-integration
"""


"""
###! First, from the `Jovian_master` env do:
###!    ./generate_sample_sheet.py private_cleaned_data > sample_sheet_ref_alignment.yaml
###! And then the actual run:
###!    snakemake -s Snakefile --latency-wait 60 --cores 300 --local-cores 4 --use-conda --cluster 'bsub -q bio -e lsf_log.stderr -o lsf_log.stdout -n {threads} -R "span[hosts=1]"'
#######! For better debugging:
###!    snakemake -s Snakefile -p --latency-wait 60 --cores 300 --local-cores 4 --use-conda --cluster 'bsub -q bio -e lsf_log.stderr -o lsf_log.stdout -n {threads} -R "span[hosts=1]"'
###! Or a shorter (using Jovian's config files):
###!     snakemake -s bin/Ref_alignment --profile config
"""


#@################################################################################
#@#### Import config file, sample_sheet and set output folder names          #####
#@################################################################################


shell.executable("/bin/bash")

configfile: "config/pipeline_parameters.yaml"
configfile: "config/variables.yaml"

import pprint
import os
import yaml
yaml.warnings({'YAMLLoadWarning': False}) # Suppress yaml "unsafe" warnings.

SAMPLES = {}
with open(config["sample_sheet_reference_alignment"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) # SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"

# This file is given as a snakemake CLI argument from within the wrapper
REFERENCE = config["reference"]

READS_INPUT_DIR = config["reference_alignment"]["input_dir"]
RESULTS_OUTPUT_DIR = config["reference_alignment"]["output_dir"]
LOGS_OUTPUT_DIR = config["reference_alignment"]["log_dir"]


#@################################################################################
#@#### Specify Jovian's final output:                                        #####
#@################################################################################


localrules: 
    all,
    RA_index_reference,
    RA_concat_BoC_metrics,


rule all:
    input:
        RESULTS_OUTPUT_DIR + "reference/" + os.path.basename(config["reference"]), # Copy of the reference fasta (for standardization and easy logging)
        RESULTS_OUTPUT_DIR + "reference/" + os.path.basename(config["reference"]) + ".1.bt2", # The bowtie2 index files (I've only specified one, but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated)
        expand("{out}alignment/{sample}_sorted.{extension}", out = RESULTS_OUTPUT_DIR, sample = SAMPLES, extension = [ 'bam', 'bam.bai' ]), # The reference alignment (bam format) files.
        expand("{out}consensus_seqs/raw/{sample}_{extension}", out = RESULTS_OUTPUT_DIR, sample = SAMPLES, extension = [ 'calls.vcf.gz', 'raw_consensus.fa' ]), # A zipped vcf file contained SNPs versus the given reference and a RAW consensus sequence, see explanation below for the meaning of RAW.
        expand("{out}consensus_seqs/{sample}.bedgraph", out = RESULTS_OUTPUT_DIR, sample = SAMPLES), # Lists the coverage of the alignment against the reference in a bedgraph format, is used to determine the coverage mask files below.
        expand("{out}consensus_seqs/{sample}_{filt_character}-filt_cov_ge_{thresholds}.fa", out = RESULTS_OUTPUT_DIR, sample = SAMPLES, filt_character = [ 'N', 'minus' ], thresholds = [ '1', '5', '10', '30', '100' ]), # Consensus sequences filtered for different coverage thresholds (1, 5, 10, 30 and 100). For each threshold two files are generated, one where failed positioned are replaced with a N nucleotide and the other where its replaced with a minus character (gap).
        expand("{out}consensus_seqs/BoC_analysis/{sample}_BoC{extension}", out = RESULTS_OUTPUT_DIR, sample = SAMPLES, extension = [ '.vcf', '_int.tsv', '_pct.tsv' ] ), # Output of the BoC analysis #TODO can probably removed after the concat rule is added.
        RESULTS_OUTPUT_DIR + "results/BoC_integer.tsv", # Integer BoC overview in .tsv format
        RESULTS_OUTPUT_DIR + "results/BoC_percentage.tsv", # Percentage BoC overview in .tsv format


#@################################################################################
#@#### Reference alignment extension processes                               #####
#@################################################################################


rule RA_index_reference:
    input:
        reference=config["reference"]
    output:
        reference_copy= RESULTS_OUTPUT_DIR + "reference/" + os.path.basename(config["reference"]),
        reference_index= RESULTS_OUTPUT_DIR + "reference/" + os.path.basename(config["reference"]) + ".1.bt2",
    conda:
        "envs/RA_ref_alignment.yaml"
    benchmark:
        LOGS_OUTPUT_DIR + "benchmark/RA_index_reference.txt"
    threads: 4
    log:
        LOGS_OUTPUT_DIR + "RA_index_reference.log"
    shell: # The reference is copied to the hardcoded subdir to make it standardized and easily logged.
        """
cp {input.reference} {output.reference_copy}
bowtie2-build --threads {threads} {output.reference_copy} {output.reference_copy} >> {log} 2>&1
        """


rule RA_align_to_reference:
    input:
        pR1= READS_INPUT_DIR + "{sample}_pR1.fq",
        pR2= READS_INPUT_DIR + "{sample}_pR2.fq",
        unpaired= READS_INPUT_DIR + "{sample}_unpaired.fq",
        reference= rules.RA_index_reference.output.reference_copy
    output:
        sorted_bam= RESULTS_OUTPUT_DIR + "alignment/{sample}_sorted.bam",
        sorted_bam_index= RESULTS_OUTPUT_DIR + "alignment/{sample}_sorted.bam.bai",
    conda:
        "envs/RA_ref_alignment.yaml"
    benchmark:
        LOGS_OUTPUT_DIR + "benchmark/RA_align_to_reference_{sample}.txt"
    threads: config["threads"]["RA_align_to_reference"]
    log:
        LOGS_OUTPUT_DIR + "RA_align_to_reference_{sample}.log"
    params:
        alignment_type="--local",
        #alignment_type=config["HuGo_removal"]["alignment_type"] # TODO voor latere versie
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


# This rule overwrites nucleotides in the ref based on the input reads and creates a RAW consensus sequence, RAW means:
#### A consensus genome is generated even if there is 0 coverage. At positions where no reads aligned it just takes
#### the reference nucleotide and inserts that into the consensus fasta. Obviously, this is not desirable and prone
#### to misinterpretation! Therefore, in a rule below, all positions with a coverage of 0 are masked, i.e. nucleotides
#### at that position are replaced with "N" or "-". Both replacements are made, since many aligners can handle gaps ("-")
#### but they cannot handle N's. However, the file with N's can be used to check for indels since this cannot be seen
#### in the file where every missing nucleotide is replaced with a "-".
#! basecalling kan niet multithreaded, hoe ze het bij HuGo doen is door het te splitten op chromosome of regio's.
#TODO bovenstaande toch even uitzoeken en bronnn aan de code toevoegen.
#TODO de output consensus fasta is nu een size-limited multi-line fasta, hier nog een two-line fasta van maken. Zie Jovian seqtk.
rule RA_extract_raw_consensus:
    input:
        bam= rules.RA_align_to_reference.output.sorted_bam,
        reference= rules.RA_index_reference.output.reference_copy,
    output:
        gzipped_vcf= RESULTS_OUTPUT_DIR + "consensus_seqs/raw/{sample}_calls.vcf.gz",
        raw_consensus_fasta= RESULTS_OUTPUT_DIR + "consensus_seqs/raw/{sample}_raw_consensus.fa",
    conda:
        "envs/RA_ref_alignment.yaml"
    benchmark:
        LOGS_OUTPUT_DIR + "benchmark/RA_extract_raw_consensus_{sample}.txt"
    threads: 1
    log:
        LOGS_OUTPUT_DIR + "RA_extract_raw_consensus_{sample}.log"
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
#TODO de fasta headers moeten nog een suffix krijgen met de cov en de filtering character, anders is het niet duidelijk bij het catten van de files.
rule RA_extract_clean_consensus:
    input:
        bam= rules.RA_align_to_reference.output.sorted_bam,
        reference= rules.RA_index_reference.output.reference_copy,
        raw_consensus= rules.RA_extract_raw_consensus.output.raw_consensus_fasta, # Only needed for when there are no positions in the bed with a coverage of 0; in that case the RAW fasta is actually suitable for downstream processes and it is simply copied.
    output:
        bedgraph= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}.bedgraph",
        filt_consensus_N_filt_ge_1= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_N-filt_cov_ge_1.fa",
        filt_consensus_N_filt_ge_5= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_N-filt_cov_ge_5.fa",
        filt_consensus_N_filt_ge_10= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_N-filt_cov_ge_10.fa",
        filt_consensus_N_filt_ge_30= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_N-filt_cov_ge_30.fa",
        filt_consensus_N_filt_ge_100= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_N-filt_cov_ge_100.fa",
        filt_consensus_minus_filt_ge_1= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_minus-filt_cov_ge_1.fa",
        filt_consensus_minus_filt_ge_5= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_minus-filt_cov_ge_5.fa",
        filt_consensus_minus_filt_ge_10= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_minus-filt_cov_ge_10.fa",
        filt_consensus_minus_filt_ge_30= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_minus-filt_cov_ge_30.fa",
        filt_consensus_minus_filt_ge_100= RESULTS_OUTPUT_DIR + "consensus_seqs/{sample}_minus-filt_cov_ge_100.fa",
    conda:
        "envs/RA_ref_alignment.yaml"
    benchmark:
        LOGS_OUTPUT_DIR + "benchmark/RA_extract_clean_consensus_{sample}.txt"
    threads: 1
    log:
        LOGS_OUTPUT_DIR + "RA_extract_clean_consensus_{sample}.log"
    params:
        output_folder= RESULTS_OUTPUT_DIR + "consensus_seqs/",
    shell:
        """
bash bin/scripts/RA_consensus_at_diff_coverages.sh {wildcards.sample} {input.bam} {input.reference} {input.raw_consensus} \
{params.output_folder} {log} >> {log} 2>&1
        """


#TODO make a python script or bash function/include to do this more efficiently, currently it's hacky, but it works
#TODO the mpileup in this rule is probably superfluous since I think that the pileup of the rule above can be reused (after minor tweaks to the BoC_script), but currently no time for it
rule RA_determine_BoC_at_diff_cov_thresholds:
    input:
        bam= rules.RA_align_to_reference.output.sorted_bam,
        reference= rules.RA_index_reference.output.reference_copy,
    output:
        BoC_vcf= RESULTS_OUTPUT_DIR + "consensus_seqs/BoC_analysis/{sample}_BoC.vcf",
        percentage_BoC_tsv= RESULTS_OUTPUT_DIR + "consensus_seqs/BoC_analysis/{sample}_BoC_pct.tsv",
        integer_BoC_tsv= RESULTS_OUTPUT_DIR + "consensus_seqs/BoC_analysis/{sample}_BoC_int.tsv",
    conda:
        "envs/RA_ref_alignment.yaml"
    benchmark:
        LOGS_OUTPUT_DIR + "RA_benchmark/determine_BoC_at_diff_cov_thresholds_{sample}.txt"
    threads: 1
    log:
        LOGS_OUTPUT_DIR + "RA_determine_BoC_at_diff_cov_thresholds_{sample}.log"
    params:
    shell: # Source: http://www.metagenomics.wiki/tools/samtools/breadth-of-coverage
        """
bash bin/scripts/RA_BoC_analysis.sh {wildcards.sample} {input.bam} {input.reference} \
{output.BoC_vcf} {output.percentage_BoC_tsv} {output.integer_BoC_tsv} >> {log} 2>&1
        """


rule RA_concat_BoC_metrics:
    input:
        BoC_int_tsv= expand("{out}consensus_seqs/BoC_analysis/{sample}_BoC_int.tsv", out = RESULTS_OUTPUT_DIR, sample = SAMPLES),
        BoC_pct_tsv= expand("{out}consensus_seqs/BoC_analysis/{sample}_BoC_pct.tsv", out = RESULTS_OUTPUT_DIR, sample = SAMPLES),
    output:
        combined_BoC_int_tsv= RESULTS_OUTPUT_DIR + "results/BoC_integer.tsv",
        combined_BoC_pct_tsv= RESULTS_OUTPUT_DIR + "results/BoC_percentage.tsv",
    conda:
        "envs/RA_ref_alignment.yaml"
    benchmark:
        LOGS_OUTPUT_DIR + "RA_benchmark/concat_BoC_metrics.txt"
    threads: 1
    log:
        LOGS_OUTPUT_DIR + "RA_concat_BoC_metrics.log"
    params:
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_int_tsv}
cat {input.BoC_int_tsv} >> {output.combined_BoC_int_tsv}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_pct_tsv}
cat {input.BoC_pct_tsv} >> {output.combined_BoC_pct_tsv}
        """
