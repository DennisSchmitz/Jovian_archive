"""
Authors:
    Dennis Schmitz, Sam Nooij, Florian Zwagemaker, Robert Verhagen,
    Jeroen Cremer, Thierry Janssens, Mark Kroon, Erwin van Wieringen,
    Annelies Kroneman, Harry Vennema, Marion Koopmans
Acknowledgement:
    Jeroen van Rooij, André Uitterlinden
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
Automation:
    iRODS automatically executes this workflow for all "vir-ngs" labelled
    Illumina sequencing runs. See the internal Gitlab repo for the wrapper
    with additional information.
    #TODO see below
    N.B. this is currently a hacky implementation, it automatically assumes
    samples are nCoV. Which is true for now but must be solved more elegantly
    at a later time.
    #TODO see above
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
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) # SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"

# The reference file is given as a snakemake CLI argument from within the wrapper, so NOT via the pipeline_parameters.yaml
reference = config["reference"]
reference_basename = os.path.splitext(os.path.basename(reference))[0]    # source: https://stackoverflow.com/questions/678236/how-to-get-the-filename-without-the-extension-from-a-path-in-python

# Set base directories
indir = config["reference_alignment"]["input_dir"]
outdir = config["reference_alignment"]["output_dir"]
logdir = config["reference_alignment"]["log_dir"]
bench = "benchmark/"

# Set dir with conda-env files
conda_envs = "envs/"
glob = "global/"
illuminaref = "Illumina_Ref/"
illuminameta = "Illlumina_meta/"
nanoref = "Nano_Ref/"

refdir = "reference/"
aln = "alignment/"
cons = "consensus/"
raw = "raw/"
boc = "BoC_analysis/"
igv = "igv/"
res = "results/"
mqc = "multiqc_data/"

igv_rep = outdir + res + "igvjs.html"
mqc_rep = outdir + res + "multiqc.html"




#@################################################################################
#@#### Specify Jovian's final output:                                        #####
#@################################################################################


localrules: 
    all,
    RA_index_reference,
    RA_determine_BoC_at_diff_cov_thresholds,
    RA_concat_BoC_metrics,
    RA_HTML_IGVJs_variable_parts,
    RA_HTML_IGVJs_generate_final

rule all:
    input:
        expand("data/cleaned_fastq/{sample}_{read}.fq",
                sample = SAMPLES,
                read = [ 'pR1', 'pR2', 'unpaired' ]
                ), # Extract unmapped & paired reads AND unpaired from HuGo alignment; i.e. cleaned fastqs #TODO omschrijven naar betere smk syntax
        expand("{out}{ref}{extension}",
                out = outdir + refdir,
                ref = reference_basename,
                extension = [ '.fasta', '.fasta.1.bt2', '.fasta.fai', '.fasta.sizes', '.windows', '_GC.bedgraph' ]
                ), # Copy of the reference file (for standardization and easy logging), bowtie2-indices (I've only specified one, but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated) and the GC-content files.
        expand("{out}{sample}_sorted.{extension}",
                out = outdir + aln,
                sample = SAMPLES,
                extension = [ 'bam', 'bam.bai', 'MarkDup_metrics' ]
                ), # The reference alignment (bam format) files.
        expand("{out}{sample}_{extension}",
                out = outdir + cons + raw,
                sample = SAMPLES,
                extension = [ 'calls.vcf.gz', 'raw_consensus.fa' ]
                ), # A zipped vcf file contained SNPs versus the given reference and a RAW consensus sequence, see explanation below for the meaning of RAW.
        expand("{out}{sample}.bedgraph",
                out = outdir + cons,
                sample = SAMPLES), # Lists the coverage of the alignment against the reference in a bedgraph format, is used to determine the coverage mask files below.
        expand("{out}{sample}_{filt_character}-filt_cov_ge_{thresholds}.fa",
                out = outdir + res + cons,
                sample = SAMPLES,
                filt_character = [ 'N', 'minus' ],
                thresholds = [ '1', '5', '10', '30', '100' ]
                ), # Consensus sequences filtered for different coverage thresholds (1, 5, 10, 30 and 100). For each threshold two files are generated, one where failed positioned are replaced with a N nucleotide and the other where its replaced with a minus character (gap).
        expand("{out}{sample}_BoC{extension}",
                out = outdir + boc,
                sample = SAMPLES,
                extension = [ '_int.tsv', '_pct.tsv' ]
                ), # Output of the BoC analysis #TODO can probably removed after the concat rule is added.
        outdir + res + "BoC_integer.tsv", # Integer BoC overview in .tsv format
        outdir + res + "BoC_percentage.tsv", # Percentage BoC overview in .tsv format
        expand("{out}{ref}_{extension}",
                out = outdir + refdir,
                ref = reference_basename,
                sample = SAMPLES,
                extension = [ 'ORF_AA.fa', 'ORF_NT.fa', 'annotation.gff', 'annotation.gff.gz', 'annotation.gff.gz.tbi' ]
                ), # Prodigal ORF prediction output, required for the IGVjs visualisation
        igv_rep, # IGVjs output html
        mqc_rep, # MultiQC report


#@################################################################################
#@#### The `onstart` checker codeblock                                       #####
#@################################################################################

onstart:
    shell("""
        mkdir -p {outdir}{res} 
        echo -e "\nLogging pipeline settings..."

        echo -e "\tGenerating methodological hash (fingerprint)..."
        echo -e "This is the link to the code used for this analysis:\thttps://github.com/DennisSchmitz/Jovian/tree/$(git log -n 1 --pretty=format:"%H")" > {outdir}{res}/log_git.txt
        echo -e "This code with unique fingerprint $(git log -n1 --pretty=format:"%H") was committed by $(git log -n1 --pretty=format:"%an <%ae>") at $(git log -n1 --pretty=format:"%ad")" >> {outdir}{res}/log_git.txt

        echo -e "\tGenerating full software list of current Conda environment (\"Jovian_master\")..."
        conda list > {outdir}{res}/log_conda.txt
        
        echo -e "\tGenerating config file log..."
        rm -f {outdir}{res}/log_config.txt
        for file in config/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> {outdir}{res}/log_config.txt
            cat ${{file}} >> {outdir}{res}/log_config.txt
            echo -e "\n\n" >> {outdir}{res}/log_config.txt
        done
    """)


#@################################################################################
#@#### Reference alignment extension processes                               #####
#@################################################################################

#! rules via include statements are shared between core workflow and RA workflow
#>############################################################################
#>#### Data quality control and cleaning                                 #####
#>############################################################################

include: "rules/QC_raw.smk"
include: "rules/CleanData.smk"
include: 'rules/QC_clean.smk'

#>############################################################################
#>#### Removal of background host data                                   #####
#>############################################################################

include: "rules/BG_removal_1.smk"
include: "rules/BG_removal_2.smk"
include: "rules/BG_removal_3.smk"


###########! nuttig om contig metrics rule ook toe te voegen?


#>############################################################################
#>#### Process the reference                                             #####
#>############################################################################
rule RA_index_reference:
    input:
        reference = reference
    output:
        reference_copy = outdir + refdir + reference_basename + ".fasta",
        reference_index = outdir + refdir + reference_basename + ".fasta.1.bt2", # I've only specified ".fasta.1.bt2", but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated. #TODO find a way to specify all output correctly (multiext snakemake syntax?)
    conda:
        conda_envs + "RA_ref_alignment.yaml"
    benchmark:
        logdir + bench + "RA_index_reference.txt"
    threads: 4
    log:
        logdir + "RA_index_reference.log"
    shell: # The reference is copied to the hardcoded subdir to make it standardized and easily logged. Convert it to a two-line fasta for easier downstream processing.
        """
cat {input.reference} | seqtk seq - > {output.reference_copy}
bowtie2-build --threads {threads} {output.reference_copy} {output.reference_copy} >> {log} 2>&1
        """


##########################!
# Nuttig voor IGVjs vis. Gejat uit Jovian core met minor changes, kunnen we waarschijlijk efficienter doen. Bijvoorbeeld door gewoon een goed gecureerde ORF annotatie toe te voegen bij starten van analyse.
rule RA_reference_ORF_analysis:
    input:
        reference = rules.RA_index_reference.output.reference_copy
    output: 
        ORF_AA_fasta = outdir + refdir + reference_basename + "_ORF_AA.fa",
        ORF_NT_fasta = outdir + refdir + reference_basename + "_ORF_NT.fa",
        ORF_annotation_gff = outdir + refdir + reference_basename + "_annotation.gff",
        zipped_gff3 = outdir + refdir + reference_basename + "_annotation.gff.gz",
        index_zipped_gff3 = outdir + refdir + reference_basename + "_annotation.gff.gz.tbi",
    conda:
        conda_envs + "scaffold_analyses.yaml"
    benchmark:
        logdir + bench + "RA_reference_ORF_analysis.txt"
    log:
        logdir + "RA_reference_ORF_analysis.log"
    threads: 1
    params: #? Currently it's using the same prodigal settings as the main workflow, I see no problems with it since it's both foremost intended for viruses.
        procedure = config["ORF_prediction"]["procedure"],
        output_format = config["ORF_prediction"]["output_format"]
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


rule RA_determine_GC_content:
    input:
        fasta = rules.RA_index_reference.output.reference_copy,
    output:
        fasta_fai = outdir + refdir + reference_basename + ".fasta.fai",
        fasta_sizes = outdir + refdir + reference_basename + ".fasta.sizes",
        bed_windows = outdir + refdir + reference_basename + ".windows",
        GC_bed = outdir + refdir + reference_basename + "_GC.bedgraph",
    conda:
        conda_envs + "scaffold_analyses.yaml"
    benchmark:
        logdir + bench + "RA_determine_GC_content.txt"
    log:
        logdir + "RA_determine_GC_content.log"
    threads: 1
    params:
        window_size = "50"
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


#>##################################################################################################
#>#### Align to ref, mark and optionally remove duplicates, call SNPs, generate new consensus  #####
#>##################################################################################################
rule RA_align_to_reference:
    input:
        pR1 = rules.HuGo_removal_pt2_extract_paired_unmapped_reads.output.fastq_R1,
        pR2 = rules.HuGo_removal_pt2_extract_paired_unmapped_reads.output.fastq_R2,
        unpaired = rules.HuGo_removal_pt3_extract_unpaired_unmapped_reads.output,
        reference = rules.RA_index_reference.output.reference_copy
    output:
        sorted_bam = outdir + aln + "{sample}_sorted.bam",
        sorted_bam_index = outdir + aln + "{sample}_sorted.bam.bai",
        dup_metrics = outdir + aln + "{sample}_sorted.MarkDup_metrics" #TODO deze toevoegen aan MultiQC?
    conda:
        conda_envs + "RA_ref_alignment.yaml"
    benchmark:
        logdir + bench + "RA_align_to_reference_{sample}.txt"
    threads: config["threads"]["RA_align_to_reference"]
    log:
        logdir + "RA_align_to_reference_{sample}.log"
    params:
        alignment_type = "--local",
        remove_dups = "-r", #! Don't change this, see this gotcha with duplicate marked reads in bedtools genomecov (which is used downstream): https://groups.google.com/forum/#!msg/bedtools-discuss/wJNC2-icIb4/wflT6PnEHQAJ . bedtools genomecov is not able to filter them out and includes those dup-reads in it's coverage metrics. So the downstream BoC analysis and consensus at diff cov processes require dups to be HARD removed.
        markdup_mode = "t",
        max_read_length = "300", # This is the default value and also the max read length of Illumina in-house sequencing.
    shell:
        """
bowtie2 --time --threads {threads} {params.alignment_type} \
-x {input.reference} \
-1 {input.pR1} \
-2 {input.pR2} \
-U {input.unpaired} 2> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools collate -@ {threads} -O - 2>> {log} |\
samtools fixmate -@ {threads} -m - - 2>> {log} |\
samtools sort -@ {threads} - -o - 2>> {log} |\
samtools markdup -@ {threads} -l {params.max_read_length} -m {params.markdup_mode} {params.remove_dups} -f {output.dup_metrics} - {output.sorted_bam} >> {log} 2>&1
samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        """


rule RA_extract_raw_consensus:
    input:
        bam = rules.RA_align_to_reference.output.sorted_bam,
        reference = rules.RA_index_reference.output.reference_copy,
    output: #TODO check if it can use bcf output instead of vcf for downstream processing, saves diskspace, but not a huge differences for small viral genomes. Would require changes in the igvjs index
        gzipped_vcf = outdir + cons + raw + "{sample}_calls.vcf.gz",
        raw_consensus_fasta = outdir + cons + raw + "{sample}_raw_consensus.fa",
    conda:
        conda_envs + "RA_ref_alignment.yaml"
    benchmark:
        logdir + bench + "RA_extract_raw_consensus_{sample}.txt"
    threads: 1 # Increasing this makes no differences for monopartite references/viruses. I think it splits different chromosomes to different threads #TODO check this in future version that is compatible with segmented viruses.
    log:
        logdir + "RA_extract_raw_consensus_{sample}.log"
    params: #TODO move this param to pipeline_variables.yaml when we assess and optimize this value for different viral families.
        calling_prior = "1.1e-3" # From manual: mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]. Also see https://samtools.github.io/bcftools/howtos/variant-calling.html --> higher value is less strict and vice versa #TODO this can be (has to be?) adapted to the different virus mutation rate. Assess later and optimize for different viral families
    shell:
        """
bcftools mpileup --threads {threads} --ignore-RG -O u -d 10000 -f {input.reference} {input.bam} 2>> {log} |\
bcftools call --threads {threads} -m --prior {params.calling_prior} --ploidy 1 -mv -O z 2>> {log} |\
bcftools norm -m -both -O z -f {input.reference} -o {output.gzipped_vcf} >> {log} 2>&1
tabix {output.gzipped_vcf} >> {log} 2>&1
cat {input.reference} 2>> {log} |\
bcftools consensus {output.gzipped_vcf} | seqtk seq - > {output.raw_consensus_fasta} 2>> {log}
        """


#TODO kijken of dit multithreaded kan worden.
rule RA_extract_clean_consensus:
    input:
        bam = rules.RA_align_to_reference.output.sorted_bam,
        raw_consensus = rules.RA_extract_raw_consensus.output.raw_consensus_fasta, # Only needed for when there are no positions in the bed with a coverage of 0; in that case the RAW fasta is actually suitable for downstream processes and it is simply copied.
    output:
        bedgraph = outdir + cons + "{sample}.bedgraph",
        filt_consensus_N_filt_ge_1 = outdir + res + cons + "{sample}_N-filt_cov_ge_1.fa",
        filt_consensus_N_filt_ge_5 = outdir + res + cons + "{sample}_N-filt_cov_ge_5.fa",
        filt_consensus_N_filt_ge_10 = outdir + res + cons + "{sample}_N-filt_cov_ge_10.fa",
        filt_consensus_N_filt_ge_30 = outdir + res + cons + "{sample}_N-filt_cov_ge_30.fa",
        filt_consensus_N_filt_ge_100 = outdir + res + cons + "{sample}_N-filt_cov_ge_100.fa",
        filt_consensus_minus_filt_ge_1 = outdir + res + cons + "{sample}_minus-filt_cov_ge_1.fa",
        filt_consensus_minus_filt_ge_5 = outdir + res + cons + "{sample}_minus-filt_cov_ge_5.fa",
        filt_consensus_minus_filt_ge_10 = outdir + res + cons + "{sample}_minus-filt_cov_ge_10.fa",
        filt_consensus_minus_filt_ge_30 = outdir + res + cons + "{sample}_minus-filt_cov_ge_30.fa",
        filt_consensus_minus_filt_ge_100 = outdir + res + cons + "{sample}_minus-filt_cov_ge_100.fa",
    conda:
        conda_envs + "RA_ref_alignment.yaml"
    benchmark:
        logdir + bench + "RA_extract_clean_consensus_{sample}.txt"
    threads: 1
    log:
        logdir + "RA_extract_clean_consensus_{sample}.log"
    params:
        output_data_folder = outdir + cons,
        output_results_folder = outdir + res + cons
    shell:
        """
bash bin/scripts/RA_consensus_at_diff_coverages.sh {wildcards.sample} {input.bam} {input.raw_consensus} \
{params.output_data_folder} {params.output_results_folder} {log} >> {log} 2>&1
        """


#TODO make a python script or bash function/include to do this more efficiently, currently it's hacky, but it works
rule RA_determine_BoC_at_diff_cov_thresholds:
    input:
        bedgraph = rules.RA_extract_clean_consensus.output.bedgraph,
        reference = rules.RA_index_reference.output.reference_copy,
    output:
        percentage_BoC_tsv = outdir + boc + "{sample}_BoC_pct.tsv",
        integer_BoC_tsv = outdir + boc + "{sample}_BoC_int.tsv",
    conda:
        conda_envs + "RA_ref_alignment.yaml"
    benchmark:
        logdir + bench + "RA_determine_BoC_at_diff_cov_thresholds_{sample}.txt"
    threads: 1
    log:
        logdir + "RA_determine_BoC_at_diff_cov_thresholds_{sample}.log"
    params:
    shell:
        """
bash bin/scripts/RA_BoC_analysis.sh {wildcards.sample} {input.bedgraph} {input.reference} \
{output.percentage_BoC_tsv} {output.integer_BoC_tsv} >> {log} 2>&1
        """


rule RA_concat_BoC_metrics:
    input:
        BoC_int_tsv = expand("{out}{sample}_BoC_int.tsv", 
                            out = outdir + boc,
                            sample = SAMPLES
                            ),
        BoC_pct_tsv = expand("{out}{sample}_BoC_pct.tsv",
                            out = outdir + boc,
                            sample = SAMPLES
                            ),
    output:
        combined_BoC_int_tsv = OUTPUT_DIR_RESULTS + "BoC_integer.tsv",
        combined_BoC_pct_tsv = OUTPUT_DIR_RESULTS + "BoC_percentage.tsv",
    conda:
        conda_envs + "RA_ref_alignment.yaml"
    benchmark:
        logdir + bench + "RA_concat_BoC_metrics.txt"
    threads: 1
    log:
        logdir + "RA_concat_BoC_metrics.log"
    params:
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_int_tsv}
cat {input.BoC_int_tsv} >> {output.combined_BoC_int_tsv}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_pct_tsv}
cat {input.BoC_pct_tsv} >> {output.combined_BoC_pct_tsv}
        """


##########################!
# Gejat uit Jovian core met minor changes, kunnen we waarschijlijk efficienter doen.
# TODO the report is still a bit dirty since we include two bowtie2 metric files:
#### TODO one for the hugo removal
#### TODO another for the ref alignment
#### TODOD hence the '-d' flag in the multiqc command based on https://multiqc.info/docs/#directory-names
rule RA_MultiQC_report:
    input:
        expand("data/FastQC_pretrim/{sample}_{read}_fastqc.zip",
                sample = SAMPLES,
                read = "R1 R2".split()
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand("data/FastQC_posttrim/{sample}_{read}_fastqc.zip",
                sample = SAMPLES,
                read = "pR1 pR2 uR1 uR2".split()
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand("logs/Clean_the_data_{sample}.log",
                sample = SAMPLES
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand("logs/HuGo_removal_pt1_alignment_{sample}.log",
                sample = SAMPLES
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand("{out}RA_align_to_reference_{sample}.log",
                out = logdir,
                sample = SAMPLES
                ), # TODO dit moet nog verbetert worden qua smk syntax
    output:
        mqc_rep,
        expand("{out}multiqc_{program}.txt",
                out = outdir + res + mqc,
                program = ['trimmomatic','bowtie2','fastqc']
                ),
    conda:
        conda_envs + glob + "MultiQC_report.yaml"
    benchmark:
        logdir + bench + "RA_MultiQC_report.txt"
    threads: 1
    params:
        config_file = "files/multiqc_config.yaml",
        output_dir = OUTPUT_DIR_RESULTS,
    log:
        logdir + "RA_MultiQC_report.log"
    shell:
        """
multiqc -d --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
        """


#@################################################################################
#@#### Make IGVjs html                                                       #####
#@################################################################################


rule RA_HTML_IGVJs_variable_parts:
    input:
        fasta = rules.RA_index_reference.output.reference_copy,
        ref_GC_bedgraph = rules.RA_determine_GC_content.output.GC_bed,
        ref_zipped_ORF_gff = rules.RA_reference_ORF_analysis.output.zipped_gff3,
        basepath_zipped_SNP_vcf = rules.RA_extract_raw_consensus.output.gzipped_vcf,
        basepath_sorted_bam = rules.RA_align_to_reference.output.sorted_bam,
    output:
        tab_output = outdir + igv + "2_tab_{sample}",
        div_output = outdir + igv + "4_html_divs_{sample}",
        js_flex_output = outdir + igv + "6_js_flex_{sample}",
    conda:
        conda_envs + glob + "data_wrangling.yaml"
    benchmark:
        logdir + bench + "RA_HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    log:
        logdir + "RA_HTML_IGVJs_variable_parts_{sample}.log"
    params:
    shell:
        """
bash bin/html/igvjs_write_tabs.sh {wildcards.sample} {output.tab_output}

bash bin/html/igvjs_write_divs.sh {wildcards.sample} {output.div_output}

bash bin/html/RA_igvjs_write_flex_js_middle.sh {wildcards.sample} {output.js_flex_output} \
{input.fasta} {input.ref_GC_bedgraph} {input.ref_zipped_ORF_gff} \
{input.basepath_zipped_SNP_vcf} {input.basepath_sorted_bam}
        """


rule RA_HTML_IGVJs_generate_final:
    input:
        expand("{out}{chunk_name}_{sample}",
                out = outdir + igv,
                chunk_name = [ '2_tab', '4_html_divs', '6_js_flex' ],
                sample = SAMPLES
                )
    output:
        igv_rep
    conda:
        conda_envs + glob + "data_wrangling.yaml"
    benchmark:
        logdir + bench + "RA_HTML_IGVJs_generate_final.txt"
    threads: 1
    log:
        logdir + "RA_HTML_IGVJs_generate_final.log"
    params:
        tab_basename = outdir + igv + "2_tab_",
        div_basename = outdir + igv + "4_html_divs_",
        js_flex_output = outdir + igv + "6_js_flex_",
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
        
        echo -e "\tRemoving temporary files..."
        if [ "{config[remove_temp]}" != "0" ]; then
            rm -rf {outdir + igv}   # Remove intermediate IGVjs html chunks.
        else
            echo -e "\t\tYou chose not to remove temp files, skipping..."
        fi

        echo -e "\tCreating symlinks for the interactive genome viewer..."
        bash bin/scripts/set_symlink.sh

        echo -e "\tGenerating Snakemake report..."
        snakemake -s bin/Ref_alignment.smk --unlock --config config --config reference={reference}
        snakemake -s bin/Ref_alignment.smk --report {outdir}{res}snakemake_report.html --config config --config reference={reference}

        echo -e "Finished"
    """)