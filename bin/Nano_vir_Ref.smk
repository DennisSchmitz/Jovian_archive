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

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file)

primerfile = config["primers"]

reference = config["reference"]
reference_basename = os.path.splitext(os.path.basename(reference))[0]


# set dirs
conda_envs = "envs/"
logdir = "logs/"
benchmarkdir = "benchmark/"

datadir = "data/"
refdir = "reference/"
trims = "trimmed/"
prdir = "primers/"
cln = "cleaned_data/"
withouthugo = "without_HuGo_removal/"
QCdata = "data/"
QChtml = "html/"
QCjson = "json/"
mapping = "mapped/"
aln = "alignment/"
bf = "bam-files/"
vf = "vcf-files/"
cons = "consensus/"
raw = "raw/"
filt = "filt/"
seqs = "sequences/"
boc = "BoC/"
res = "results/"
igv = "igv/"
mqc = "multiqc_data/"
MULTIQC_OUTPUT = "results/multiqc.html"


#@################################################################################
#@#### Jovian processes                                                      #####
#@################################################################################

rule all:
    input: 
        expand("{path}{sample}.fasta",
                path = datadir + cons + raw, 
                sample = SAMPLES
                ),
        expand("{path}{sample}.bedgraph", 
                path = datadir + cons + filt,
                sample = SAMPLES
                ),
        expand("{path}BoC_int.tsv",
                path = res
                ),
        expand("{path}BoC_pct.tsv",
                path = res
                ),
        expand("{path}IGVjs.html",
                path = res
                ),
        expand("{path}multiqc.html",
                path = res
                )

    #>############################################################################
    #>#### Data quality control and cleaning                                 #####
    #>############################################################################

rule Index_ref:
    input:
        ref = reference
    output:
        refcopy =  datadir + refdir + reference_basename + ".fasta",
        refcopy_index = datadir + refdir + reference_basename + ".fasta.1.bt2",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    log:
        logdir + "Index_ref.log"
    benchmark:
        logdir + benchmarkdir + "Index_ref.txt"
    threads: 4
    shell:
        """
cat {input.ref} | seqtk seq - > {output.refcopy}
bowtie2-build --threads {threads} {output.refcopy} {output.refcopy} >> {log} 2>&1
        """


rule Adapter_trimming:
    input: 
        lambda wildcards: SAMPLES[wildcards.sample]
    output: 
        trimmeddata = datadir + trims + "{sample}.fastq"
    conda:
        conda_envs + "QC_and_clean.yaml"
    benchmark:
        logdir + benchmarkdir + "Adapter_trimming_{sample}.txt"
    log:
        logdir + "Adapter_trimming_{sample}.log"
    threads: 26
    shell:
        """
porechop \
-i {input} \
-o {output.trimmeddata} \
--threads {threads} > {log} 2>&1
        """

rule Cut_primers:
    input:
        fastq = rules.Adapter_trimming.output.trimmeddata,
        primers = primerfile,
    output:
        cleaneddata_pt1 = datadir + cln + prdir + "{sample}.fastq"
    conda:
        conda_envs + "QC_and_clean.yaml"
    benchmark:
        logdir + benchmarkdir + "Primer_removal_{sample}.txt"
    log:
        logdir + "Primer_removal_{sample}.log"
    threads: 26
    params:
        primer_cutoff_plus = config["Nanopore_ref"]["Primer_cutoff_plus"],
        primer_cutoff_minus = config["Nanopore_ref"]["Primer_cutoff_minus"]
    shell:
        """
cutadapt \
--cores={threads} \
--cut {params.primer_cutoff_plus} \
--cut {params.primer_cutoff_minus} \
--revcomp -b file:{input.primers} \
-o {output.cleaneddata_pt1} {input.fastq} > {log} 2>&1
        """

rule Cleanup:
    input:
        fastq = rules.Cut_primers.output.cleaneddata_pt1
    output:
        fastq = datadir + cln + QCdata + "{sample}.fastq",
        html = datadir + cln + QChtml + "{sample}.html",
        json = datadir + cln + QCjson + "{sample}.fastp.json",
    conda:
        conda_envs + "QC_and_clean.yaml"
    benchmark:
        logdir + benchmarkdir + "Data_Cleanup_{sample}.txt"
    log:
        logdir + "Data_Cleanup_{sample}.log"
    threads: 26
    params:
        QualityFilter = config["Nanopore_ref"]["Quality_score"]
    shell:
        """
fastp \
-i {input.fastq} \
-q {params.QualityFilter} \
-o {output.fastq} \
-h {output.html} \
-j {output.json} > {log} 2>&1
        """


rule Hugo_removal_pt1:
    input:
        bg = config["databases"]["background_ref"],
        unmapped_fastq = rules.Cleanup.output.fastq,
    output: 
        sorted_bam = datadir + cln + withouthugo + mapping + "{sample}.bam",
        sorted_bam_index = datadir + cln + withouthugo + mapping + "{sample}.bam.bai",
    conda:
        conda_envs + "HuGo_removal.yaml"
    benchmark:
        logdir + benchmarkdir + "Hugo_removal_pt1_{sample}.txt"
    log:
        logdir + "Hugo_removal_pt1_{sample}.log"
    threads: config["threads"]["HuGo_removal"]
    shell: 
        """
minimap2 -ax map-ont {input.bg} {input.unmapped_fastq} 2>> {log} |\
samtools view -uS 2>> {log} |\
samtools sort -o {output.sorted_bam} 2>> {log}
samtools index {output.sorted_bam} >> {log} 2>&1
        """

rule Hugo_removal_pt2:
    input:
        bam = rules.Hugo_removal_pt1.output.sorted_bam,
        bam_index = rules.Hugo_removal_pt1.output.sorted_bam_index,
    output:
        cleanedfastq = datadir + cln + "{sample}.fq"
    conda:
        conda_envs + "HuGo_removal.yaml"
    benchmark:
        logdir + benchmarkdir + "Hugo_removal_pt2_{sample}.txt"
    log:
        logdir + "Hugo_removal_pt2_{sample}.log"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
samtools view -b -F 1 -f 4 {input.bam} 2>> {log} |\
samtools sort -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output.cleanedfastq} >> {log} 2>&1
        """


rule Align_to_reference_pt1:
    input:
        ref = reference,
        fastq = rules.Hugo_removal_pt2.output.cleanedfastq,
    output:
        bam = datadir + aln + bf + "{sample}.bam",
        indexed_bam = datadir + aln + bf + "{sample}.bam.bai",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + benchmarkdir + "Align_to_reference_pt1_{sample}.txt"
    log:
        logdir + "Align_to_reference_pt1_{sample}.log"
    threads: 26
    shell:
        """
minimap2 -ax map-ont {input.ref} {input.fastq} 2>> {log} |\
samtools view -@ {threads} -F 256 -F 512 -F 4 -F 2048 -uS 2>> {log} |\
samtools sort -o {output.bam} >> {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
        """ 

rule Align_to_reference_pt2:
    input:
        ref = rules.Index_ref.output.refcopy,
        bam = rules.Align_to_reference_pt1.output.bam,
    output:
        vcf = datadir + aln + vf + "{sample}.vcf.gz",
        vcf_index = datadir + aln + vf + "{sample}.vcf.gz.csi",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + benchmarkdir + "Align_to_reference_pt2_{sample}.txt"
    log:
        logdir + "Align_to_reference_pt2_{sample}.log"
    threads: 26
    shell:
        """
bcftools mpileup --ignore-RG -Ou -d 10000 -f {input.ref} {input.bam} 2>> {log} |\
bcftools call --ploidy 1 -mv -Oz 2>> {log} |\
bcftools norm -m -both -O z -f {input.ref} -o {output.vcf} >> {log} 2>&1
tabix {output.vcf} >> {log} 2>&1
bcftools index {output.vcf}
        """ 

rule Create_raw_consensus:
    input:
        ref = reference,
        vcf = rules.Align_to_reference_pt2.output.vcf,
    output: 
        consensus = datadir + cons + raw + "{sample}.fasta"
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + benchmarkdir + "Create_raw_consensus_{sample}.txt"
    log:
        logdir + "Create_raw_consensus_{sample}.log"
    threads: 26
    shell:
        """
cat {input.ref} | bcftools consensus {input.vcf} 1> {output.consensus} 2> {log}
        """

rule extract_cleaned_consensus:
    input:
        raw_consensus = rules.Create_raw_consensus.output.consensus,
        bam = rules.Align_to_reference_pt1.output.bam,
    output:
        bedgraph = datadir + cons + filt + "{sample}.bedgraph",
        filt_consensus_N_filt_ge_1 = datadir + cons + seqs + "{sample}_N-filt_cov_ge_1.fa",
        filt_consensus_N_filt_ge_5 = datadir + cons + seqs + "{sample}_N-filt_cov_ge_5.fa",
        filt_consensus_N_filt_ge_10 = datadir + cons + seqs + "{sample}_N-filt_cov_ge_10.fa",
        filt_consensus_N_filt_ge_30 = datadir + cons + seqs + "{sample}_N-filt_cov_ge_30.fa",
        filt_consensus_N_filt_ge_100 = datadir + cons + seqs + "{sample}_N-filt_cov_ge_100.fa",
        filt_consensus_minus_filt_ge_1 = datadir + cons + seqs + "{sample}_minus-filt_cov_ge_1.fa",
        filt_consensus_minus_filt_ge_5 = datadir + cons + seqs + "{sample}_minus-filt_cov_ge_5.fa",
        filt_consensus_minus_filt_ge_10 = datadir + cons + seqs + "{sample}_minus-filt_cov_ge_10.fa",
        filt_consensus_minus_filt_ge_30 = datadir + cons + seqs + "{sample}_minus-filt_cov_ge_30.fa",
        filt_consensus_minus_filt_ge_100 = datadir + cons + seqs + "{sample}_minus-filt_cov_ge_100.fa",
    params:
        output_data_folder = datadir + cons + filt,
        output_results_folder = datadir + cons + seqs
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    threads: 26
    log:
        logdir + "Extract_cleaned_consensus_{sample}.log"
    benchmark:
        logdir + benchmarkdir + "Extract_cleaned_consensus_{sample}.txt"
    shell:
        """
bash bin/scripts/RA_consensus_at_diff_coverages.sh {wildcards.sample} {input.bam} {input.raw_consensus} \
{params.output_data_folder} {params.output_results_folder} {log} >> {log} 2>&1
        """

rule calculate_BoC:
    input:
        bedgraph = rules.extract_cleaned_consensus.output.bedgraph,
        reference = rules.Index_ref.output.refcopy,
    output:
        pct_boc_tsv = datadir + cons + boc + "{sample}_BoC_pct.tsv",
        int_boc_tsv = datadir + cons + boc + "{sample}_BoC_int.tsv",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + benchmarkdir + "Calculate_boc_{sample}.txt"
    log:
        logdir + "Calculate_boc_{sample}.log"
    threads: 1
    shell:
        """
bash bin/scripts/RA_BoC_analysis.sh {wildcards.sample} {input.bedgraph} {input.reference} \
{output.pct_boc_tsv} {output.int_boc_tsv} >> {log} 2>&1
        """

rule concat_boc:
    input:
        boc_int = expand("{path}{sample}_BoC_int.tsv",
                            path = datadir + cons + boc,
                            sample = SAMPLES
                            ),
        boc_pct = expand("{path}{sample}_BoC_pct.tsv",
                            path = datadir + cons + boc,
                            sample = SAMPLES
                            ),
    output:
        conc_boc_int = res + "BoC_int.tsv",
        conc_boc_pct = res + "BoC_pct.tsv",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + benchmarkdir + "Concat_boc.txt"
    log:
        logdir + "Concat_boc.log"
    threads: 1
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.conc_boc_int}
cat {input.boc_int} >> {output.conc_boc_int}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.conc_boc_pct}
cat {input.boc_pct} >> {output.conc_boc_pct}
        """ 

rule ORF_Analysis:
    input: 
        ref = rules.Index_ref.output.refcopy
    output:
        ORF_AA = datadir + refdir + reference_basename + "_ORF_AA.fa",
        ORF_NT = datadir + refdir + reference_basename + "_ORF_NT.fa",
        ORF_gff = datadir + refdir + reference_basename + "_annotation.gff",
        gff_zip = datadir + refdir + reference_basename + "_annotation.gff.gz",
        gff_ind = datadir + refdir + reference_basename + "_annotation.gff.gz.tbi",
    conda:
        conda_envs + "Sequence_analysis.yaml"
    log:
        logdir + "ORF_Analysis.log"
    benchmark:
        logdir + benchmarkdir + "ORF_Analysis.txt"
    threads: 1
    params:
        procedure = config["Global"]["ORF_procedure"],
        output_format = config["Global"]["ORF_output_format"]
    shell:
        """
prodigal -q -i {input.ref} \
-a {output.ORF_AA} \
-d {output.ORF_NT} \
-o {output.ORF_gff} \
-p {params.procedure} \
-f {params.output_format} > {log} 2>&1
bgzip -c {output.ORF_gff} 2>> {log} 1> {output.gff_zip}
tabix -p gff {output.gff_zip} >> {log} 2>&1
        """ 

rule determine_GC_content:
    input:
        ref = rules.Index_ref.output.refcopy,
    output:
        fai = datadir + refdir + reference_basename + ".fasta.fai",
        sizes = datadir + refdir + reference_basename + ".fasta.sizes",
        bed_windows = datadir + refdir + reference_basename + ".windows",
        GC_bed = datadir + refdir + reference_basename + "_GC.bedgraph",
    conda:
        conda_envs + "Sequence_analysis.yaml"
    log:
        logdir + "Determine_GC_content.log"
    benchmark:
        logdir + benchmarkdir + "Determine_GC_content.txt"
    threads: 1
    params:
        window_size = config["Global"]["GC_window_size"]
    shell:
        """
samtools faidx -o {output.fai} {input.ref} > {log} 2>&1
cut -f 1,2 {output.fai} 2> {log} 1> {output.sizes}
bedtools makewindows \
-g {output.sizes} \
-w {params.window_size} 2>> {log} 1> {output.bed_windows}
bedtools nuc \
-fi {input.ref} \
-bed {output.bed_windows} 2>> {log} |\
cut -f 1-3,5 2>> {log} 1> {output.GC_bed}
        """ 


rule HTML_IGVJs_variable_parts:
    input:
        ref = rules.Index_ref.output.refcopy,
        GC_bed = rules.determine_GC_content.output.GC_bed,
        ORF_gff = rules.ORF_Analysis.output.gff_zip,
        vcf = rules.Align_to_reference_pt2.output.vcf,
        bam = rules.Align_to_reference_pt1.output.bam,
    output:
        tab = datadir + igv + "2_tab_{sample}",
        div = datadir + igv + "4_html_divs_{sample}",
        js = datadir + igv + "6_js_flex_{sample}",
    conda:
        conda_envs + "data_wrangling.yaml"
    log:
        logdir + "HTML_IGVJs_variable_parts_{sample}.log"
    benchmark:
        logdir + benchmarkdir + "HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_tabs.sh {wildcards.sample} {output.tab}

bash bin/html/igvjs_write_divs.sh {wildcards.sample} {output.div}

bash bin/html/igvjs_write_flex_js_middle.sh {wildcards.sample} {output.js} \
{input.ref} {input.GC_bed} {input.ORF_gff} \
{input.vcf} {input.bam}
        """

rule HTML_IGVJs_generate_file:
    input:
        expand("{path}{b}_{sample}",
                path = datadir + igv,
                b = [ '2_tab', '4_html_divs', '6_js_flex' ],
                sample = SAMPLES
                )
    output:
        html = res + "IGVjs.html"
    conda:
        conda_envs + "data_wrangling.yaml"
    log:
        logdir + "HTML_IGVJs_generate_file.log"
    benchmark:
        logdir + benchmarkdir + "HTML_IGVJs_generate_file.txt"
    threads: 1
    params:
        tab_basename = datadir + igv + "2_tab_",
        div_basename = datadir + igv + "4_html_divs_",
        js_flex_output = datadir + igv + "6_js_flex_",
    shell:
        """
cat files/html_chunks/1_header.html > {output.html}
cat {params.tab_basename}* >> {output.html}
cat files/html_chunks/3_tab_explanation_RA.html >> {output.html}
cat {params.div_basename}* >> {output.html}
cat files/html_chunks/5_js_begin.html >> {output.html}
cat {params.js_flex_output}* >> {output.html}
cat files/html_chunks/7_js_end.html >> {output.html}
        """ 

rule MultiQC_report:
    input: 
        #expand("{l}Cleanup_{sample}.log", l = logdir, sample = SAMPLES),
        expand("{path}{sample}.fastp.json", path = datadir + cln + QCjson, sample = SAMPLES),
        expand("{l}Primer_removal_{sample}.log", l = logdir, sample = SAMPLES),
    output: 
        MULTIQC_OUTPUT,
        #expand("{out}{program}_multiqc.txt", out = res + mqc, program = ['fastp','cutadapt']),
    conda:
        conda_envs + "MultiQC_report.yaml"
    log:
        logdir + "MultiQC_report.log"
    benchmark:
        logdir + benchmarkdir + "MultiQC_report.txt"
    threads: 1
    params:
        config_file = "files/multiqc_config.yaml",
        output_dir = res,
    shell:
        """
multiqc -d --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
        """


onsuccess:
    shell("""
        echo -e "\tCreating symlinks for the interactive genome viewer..."
        bash bin/scripts/set_symlink.sh
    """)