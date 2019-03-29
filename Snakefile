"""
Authors: Dennis Schmitz, Robert Verhagen, Sam Nooij, Thierry Janssens, Jeroen Cremer, Mark Kroon
Organisation: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Virology - Emerging and Endemic Viruses (EEV)
Date: 23-08-2018
Changelog, examples, installation guide and explanation on:
   https://github.com/DennisSchmitz/Jovian
"""

#################################################################################
##### The `onstart` checker codeblock                                       #####
#################################################################################

shell.executable("/bin/bash")

onstart:
    shell("echo -e '\n\tSnakemake is starting...\n\t\tPlaceholder for required configfile checker code.\n\t\tPlacholder for .ncbirc BLAST db aliases.\n\t\tPlaceholder for Jupyter config checker.\n\t\tPlaceholder Jupyter Notebook theme settings\n'")

#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

configfile: "profile/pipeline_parameters.yaml"

import pprint
import yaml

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) # SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"

# Pretty-print sample sheet for debugging / logging.
#print("Samples:")
#pprint.pprint(SAMPLES)

#################################################################################
##### Specify Jovian's final output:                                        #####
#################################################################################

localrules: 
    all,
    Generate_index_html,
    Generate_IGVjs_html_file,
    quantify_output,
    Viral_typing,
    Concat_files,
    Concat_TT_files,
    Concat_filtered_SNPs,

rule all:
    input:
        expand("data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_{read}.fastq", sample = SAMPLES, read = [ 'pR1', 'pR2', 'uR1', 'uR2' ]), # Trimmomatic output, N.B. this still contains HuGo data!
        expand("data/cleaned_fastq/{sample}_{read}.fq", sample = SAMPLES, read = [ 'pR1', 'pR2', 'unpaired' ]), # Extract unmapped & paired reads AND unpaired from HuGo alignment; i.e. cleaned fastqs
        expand("data/scaffolds_raw/{sample}/scaffolds.fasta", sample = SAMPLES), # SPAdes assembly output
        expand("data/scaffolds_filtered/{sample}_scaffolds_ge{len}nt.{extension}", sample = SAMPLES, len = config["scaffold_minLen_filter"]["minlen"], extension = [ 'fasta', 'fasta.fai', 'fasta.sizes' ]), # Filtered SPAdes Scaffolds
        expand("data/scaffolds_filtered/{sample}_sorted.bam", sample = SAMPLES), # BWA mem alignment for fragment length analysis
        expand("data/scaffolds_filtered/{sample}_insert_size_{extension}", sample = SAMPLES, extension = [ 'metrics.txt', 'histogram.pdf' ]), # Picard's CollectInsertSizeMetrics output txt and pdf
        expand("data/scaffolds_filtered/{sample}_{type}.stats", sample = SAMPLES, type = [ 'MinLenFiltSummary', 'perMinLenFiltScaffold', 'perORFcoverage' ]), ### TEMP
        expand("data/scaffolds_filtered/{sample}_{extension}", sample = SAMPLES, extension = [ 'ORF_AA.fa', 'ORF_NT.fa', 'annotation.gff', 'annotation.gff.gz', 'annotation.gff.gz.tbi' ]), # Prodigal ORF prediction output
        expand("data/scaffolds_filtered/{sample}_{extension}", sample = SAMPLES, extension = [ 'unfiltered.vcf', 'filtered.vcf', 'filtered.vcf.gz', 'filtered.vcf.gz.tbi' ]), # SNP calling output
        expand("data/scaffolds_filtered/{sample}{extension}", sample = SAMPLES, extension = [ '.windows', '_GC.bedgraph' ]), # Percentage GC content per specified window
        expand("data/scaffolds_filtered/{sample}_IGVjs.html", sample = SAMPLES), # IGVjs html's
        expand("data/taxonomic_classification/{sample}.blastn", sample = SAMPLES), # MegablastN output
        expand("data/taxonomic_classification/{sample}.{extension}", sample = SAMPLES, extension = [ 'taxtab', 'taxMagtab' ]), # Krona output
        expand("data/tables/{sample}_{extension}", sample = SAMPLES, extension = [ 'taxClassified.tsv', 'taxUnclassified.tsv', 'virusHost.tsv' ]), # Tab seperated tables with merged data
        expand("data/virus_typing_tables/{sample}_{virus}.{extension}", sample = SAMPLES, virus = [ 'NoV', 'EV' ], extension = [ 'fa', 'xml', 'csv' ]), # Virus typingtool output tables
        expand("results/{file}", file = [ 'all_taxClassified.tsv', 'all_taxUnclassified.tsv', 'all_virusHost.tsv', 'all_NoV-TT.tsv', 'all_EV-TT.tsv', 'all_filtered_SNPs.tsv' ]), # Concatenated classification, virus host and typing tool tables
        expand("results/{file}", file = [ 'heatmaps/Superkingdoms_heatmap.html', 'Sample_composition_graph.html', 'Taxonomic_rank_statistics.tsv', 'Virus_rank_statistics.tsv', 'Phage_rank_statistics.tsv', 'Bacteria_rank_statistics.tsv' ]), # Taxonomic profile and heatmap output
        expand("results/heatmaps/Virus_{rank}_heatmap.html", rank=[ "order", "family", "genus", "species" ]), # Virus (excl. phages) order|family|genus|species level heatmap for the entire run
        expand("results/heatmaps/Phage_{rank}_heatmap.html", rank=[ "order", "family", "genus", "species" ]), # Phage order|family|genus|species heatmaps for the entire run (based on a selection of phage families)
        expand("results/heatmaps/Bacteria_{rank}_heatmap.html", rank=[ "phylum", "class", "order", "family", "genus", "species" ]), # Bacteria phylum|class|order|family|genus|species level heatmap for the entire run
        expand("results/{file}.html", file = [ 'multiqc', 'krona', 'Heatmap_index', 'IGVjs_index' ]), # Reports and heatmap and IGVjs index.html

#################################################################################
##### Jovian sub-processes                                                  #####
#################################################################################

    #############################################################################
    ##### Data quality control and cleaning                                 #####
    #############################################################################

rule QC_raw_data:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][wildcards.read]
    output:
        html=temp("data/FastQC_pretrim/{sample}_{read}_fastqc.html"),
        zip=temp("data/FastQC_pretrim/{sample}_{read}_fastqc.zip")
    conda:
        "envs/QC_and_clean.yaml"
    benchmark:
        "logs/benchmark/QC_raw_data_{sample}_{read}.txt"
    threads: 1
    log:
        "logs/QC_raw_data_{sample}_{read}.log"
    shell:
        """
fastqc --quiet --outdir data/FastQC_pretrim/ {input} > {log} 2>&1
        """

rule Clean_the_data:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][i] for i in ("R1", "R2"))
    output:
        r1="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR1.fastq",
        r2="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR2.fastq",
        r1_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR1.fastq",
        r2_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR2.fastq",
    conda:
        "envs/QC_and_clean.yaml"
    benchmark:
        "logs/benchmark/Clean_the_data_{sample}.txt"
    threads: config["threads"]["Clean_the_data"]
    log:
        "logs/Clean_the_data_{sample}.log"
    params:
        adapter_removal_config=config["Trimmomatic"]["adapter_removal_config"],
        quality_trimming_config=config["Trimmomatic"]["quality_trimming_config"],
        minimum_length_config=config["Trimmomatic"]["minimum_length_config"],
    shell:
        """
trimmomatic PE -threads {threads} \
{input[0]:q} {input[1]:q} \
{output.r1} {output.r1_unpaired} \
{output.r2} {output.r2_unpaired} \
{params.adapter_removal_config} \
{params.quality_trimming_config} \
{params.minimum_length_config} > {log} 2>&1
touch -r {output.r1} {output.r1_unpaired}
touch -r {output.r2} {output.r2_unpaired}
        """

rule QC_clean_data:
    input:
        "data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_{read}.fastq"
    output:
        html=temp("data/FastQC_posttrim/{sample}_{read}_fastqc.html"),
        zip=temp("data/FastQC_posttrim/{sample}_{read}_fastqc.zip")
    conda:
        "envs/QC_and_clean.yaml"
    benchmark:
        "logs/benchmark/QC_clean_data_{sample}_{read}.txt"
    threads: 1
    log:
        "logs/QC_clean_data_{sample}_{read}.log"
    shell:
        """
if [ -s "{input}" ] # If file exists and is NOT empty (i.e. filesize > 0) do...
then
    fastqc --quiet --outdir data/FastQC_posttrim/ {input} > {log} 2>&1
else
    touch {output.html}
    touch {output.zip}
fi
        """

    #############################################################################
    ##### Removal of human (privacy sensitive) host data                    #####
    #############################################################################
    
rule HuGo_removal_pt1_alignment:
    input:
        HuGo_ref=config["databases"]["HuGo_ref"],
        r1="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR1.fastq",
        r2="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR2.fastq",
        r1_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR1.fastq",
        r2_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR2.fastq",
    output:
        sorted_bam=temp("data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam"),
        sorted_bam_index=temp("data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai"),
    conda:
        "envs/HuGo_removal.yaml"
    benchmark:
        "logs/benchmark/HuGo_removal_pt1_alignment_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    params:
        alignment_type=config["HuGo_removal"]["alignment_type"]
    log:
        "logs/HuGo_removal_pt1_alignment_{sample}.log"
    shell:
        """
bowtie2 --time --threads {threads} {params.alignment_type} \
-x {input.HuGo_ref} \
-1 {input.r1} \
-2 {input.r2} \
-U {input.r1_unpaired} \
-U {input.r2_unpaired} 2> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.sorted_bam} >> {log} 2>&1
samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        """

rule HuGo_removal_pt2_extract_paired_unmapped_reads:
    input:
         bam="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam",
         bam_index="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai",
    output:
         fastq_R1="data/cleaned_fastq/{sample}_pR1.fq",
         fastq_R2="data/cleaned_fastq/{sample}_pR2.fq",
    conda:
        "envs/HuGo_removal.yaml"
    log:
        "logs/HuGo_removal_pt2_extract_paired_unmapped_reads_{sample}.log"
    benchmark:
        "logs/benchmark/HuGo_removal_pt2_extract_paired_unmapped_reads_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
samtools view -@ {threads} -b -f 1 -f 4 -f 8 {input.bam} 2> {log} |\
samtools sort -@ {threads} -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output.fastq_R1} -fq2 {output.fastq_R2} >> {log} 2>&1
        """

rule HuGo_removal_pt3_extract_unpaired_unmapped_reads:
    input:
        bam="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam",
        bam_index="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai"
    output:
        "data/cleaned_fastq/{sample}_unpaired.fq"
    conda:
        "envs/HuGo_removal.yaml"
    log:
        "logs/HuGo_removal_pt3_extract_unpaired_unmapped_reads_{sample}.log"
    benchmark:
        "logs/benchmark/HuGo_removal_pt3_extract_unpaired_unmapped_reads_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
samtools view -@ {threads} -b -F 1 -f 4 {input.bam} 2> {log} |\
samtools sort -@ {threads} -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output} >> {log} 2>&1
        """

    #############################################################################
    ##### De novo assembly and filtering                                    #####
    #############################################################################

rule De_novo_assembly:
    input:
        fastq_pR1="data/cleaned_fastq/{sample}_pR1.fq",
        fastq_pR2="data/cleaned_fastq/{sample}_pR2.fq",
        fastq_unpaired="data/cleaned_fastq/{sample}_unpaired.fq"
    output:
        all_scaffolds="data/scaffolds_raw/{sample}/scaffolds.fasta",
        filt_scaffolds="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
    conda:
        "envs/de_novo_assembly.yaml"
    benchmark:
        "logs/benchmark/De_novo_assembly_{sample}.txt"
    threads: config["threads"]["De_novo_assembly"]
    log:
        "logs/De_novo_assembly_{sample}.log"
    params:
        max_GB_RAM="100",
        kmersizes=config["SPAdes"]["kmersizes"],
        outputfoldername="data/scaffolds_raw/{sample}/",
        minlength=config["scaffold_minLen_filter"]["minlen"],
    shell:
        """
spades.py --meta \
-1 {input.fastq_pR1} \
-2 {input.fastq_pR2} \
-s {input.fastq_unpaired} \
-t {threads} \
-m {params.max_GB_RAM} \
-k {params.kmersizes} \
-o {params.outputfoldername} > {log} 2>&1
seqtk seq {output.all_scaffolds} 2>> {log} |\
gawk -F "_" '/^>/ {{if ($4 >= {params.minlength}) {{print $0; getline; print $0}};}}' 2>> {log} 1> {output.filt_scaffolds} 
        """

    #############################################################################
    ##### Scaffold analyses: fragment length analysis, SNP-calling,        #####
    #####                    ORF prediction and QC-metrics                 #####
    #############################################################################

rule Fragment_length_analysis:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        pR1="data/cleaned_fastq/{sample}_pR1.fq",
        pR2="data/cleaned_fastq/{sample}_pR2.fq",
    output:
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        bam_bai="data/scaffolds_filtered/{sample}_sorted.bam.bai",
        txt="data/scaffolds_filtered/{sample}_insert_size_metrics.txt",
        pdf="data/scaffolds_filtered/{sample}_insert_size_histogram.pdf"
    conda:
        "envs/scaffold_analyses.yaml"
    log:
        "logs/Fragment_length_analysis_{sample}.log"
    benchmark:
        "logs/benchmark/Fragment_length_analysis_{sample}.txt"
    threads: config["threads"]["Fragment_length_analysis"]
    shell:
        """
bwa index {input.fasta} > {log} 2>&1
bwa mem -t {threads} {input.fasta} \
{input.pR1} \
{input.pR2} 2>> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.bam} >> {log} 2>&1
samtools index -@ {threads} {output.bam} >> {log} 2>&1
picard -Dpicard.useLegacyParser=false CollectInsertSizeMetrics \
-I {output.bam} \
-O {output.txt} \
-H {output.pdf} >> {log} 2>&1
        """

rule SNP_calling:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        bam_bai="data/scaffolds_filtered/{sample}_sorted.bam.bai"
    output:
        fasta_fai="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.fai" % config["scaffold_minLen_filter"]["minlen"],
        unfilt_vcf="data/scaffolds_filtered/{sample}_unfiltered.vcf",
        filt_vcf="data/scaffolds_filtered/{sample}_filtered.vcf",
        zipped_vcf="data/scaffolds_filtered/{sample}_filtered.vcf.gz",
        zipped_vcf_index="data/scaffolds_filtered/{sample}_filtered.vcf.gz.tbi"
    conda:
        "envs/scaffold_analyses.yaml"
    log:
        "logs/SNP_calling_{sample}.log"
    benchmark:
        "logs/benchmark/SNP_calling_{sample}.txt"
    threads: config["threads"]["SNP_calling"]
    params:
        max_cov=config["SNP_calling"]["max_cov"],
        minimum_AF=config["SNP_calling"]["minimum_AF"]
    shell:
        """
samtools faidx -o {output.fasta_fai} {input.fasta} > {log} 2>&1
lofreq call-parallel -d {params.max_cov} \
--no-default-filter \
--pp-threads {threads} \
-f {input.fasta} \
-o {output.unfilt_vcf} \
{input.bam} >> {log} 2>&1
lofreq filter -a {params.minimum_AF} \
-i {output.unfilt_vcf} \
-o {output.filt_vcf} >> {log} 2>&1
bgzip -c {output.filt_vcf} 2>> {log} 1> {output.zipped_vcf}
tabix -p vcf {output.zipped_vcf} >> {log} 2>&1
        """

rule ORF_analysis:
    input:
        "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"]
    output:
        ORF_AA_fasta="data/scaffolds_filtered/{sample}_ORF_AA.fa",
        ORF_NT_fasta="data/scaffolds_filtered/{sample}_ORF_NT.fa",
        ORF_annotation_gff="data/scaffolds_filtered/{sample}_annotation.gff",
        zipped_gff3="data/scaffolds_filtered/{sample}_annotation.gff.gz",
        index_zipped_gff3="data/scaffolds_filtered/{sample}_annotation.gff.gz.tbi",
    conda:
        "envs/scaffold_analyses.yaml"
    log:
        "logs/ORF_prediction_{sample}.log"
    benchmark:
        "logs/benchmark/ORF_prediction_{sample}.txt"
    threads: 1
    params:
        procedure=config["ORF_prediction"]["procedure"],
        output_format=config["ORF_prediction"]["output_format"]
    shell:
        """
prodigal -q -i {input} \
-a {output.ORF_AA_fasta} \
-d {output.ORF_NT_fasta} \
-o {output.ORF_annotation_gff} \
-p {params.procedure} \
-f {params.output_format} > {log} 2>&1
bgzip -c {output.ORF_annotation_gff} 2>> {log} 1> {output.zipped_gff3}
tabix -p gff {output.zipped_gff3} >> {log} 2>&1
        """

rule Generate_contigs_metrics:
    input:
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        ORF_NT_fasta="data/scaffolds_filtered/{sample}_ORF_NT.fa",
    output:
        summary=temp("data/scaffolds_filtered/{sample}_MinLenFiltSummary.stats"),
        perScaffold=temp("data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats"),
        perORFcoverage=temp("data/scaffolds_filtered/{sample}_perORFcoverage.stats"),
    conda:
        "envs/scaffold_analyses.yaml"
    log:
        "logs/Generate_contigs_metrics_{sample}.log"
    benchmark:
        "logs/benchmark/Generate_contigs_metrics_{sample}.txt"
    params:
        ""
    threads: 1
    shell:
        """
pileup.sh in={input.bam} \
ref={input.fasta} \
fastaorf={input.ORF_NT_fasta} \
outorf={output.perORFcoverage} \
out={output.perScaffold} 2> {output.summary} 1> {log}
        """

    #############################################################################
    ##### Determine the GC content and DoC per scaffolds                    #####
    #############################################################################

# Hier nog iets voor maken met IGVjs ipv plotly

rule Determine_GC_content:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        fasta_fai="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.fai" % config["scaffold_minLen_filter"]["minlen"],
    output:
        fasta_sizes=temp("data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.sizes" % config["scaffold_minLen_filter"]["minlen"]),
        bed_windows=temp("data/scaffolds_filtered/{sample}.windows"),
        GC_bed="data/scaffolds_filtered/{sample}_GC.bedgraph"
    conda:
        "envs/scaffold_analyses.yaml"
    log:
        "logs/Determine_GC_content_{sample}.log"
    benchmark:
        "logs/benchmark/Determine_GC_content_{sample}.txt"
    threads: 1
    params:
        window_size="50"
    shell:
        """
cut -f 1,2 {input.fasta_fai} 2> {log} 1> {output.fasta_sizes}
bedtools makewindows \
-g {output.fasta_sizes} \
-w {params.window_size} 2>> {log} 1> {output.bed_windows}
bedtools nuc \
-fi {input.fasta} \
-bed {output.bed_windows} 2>> {log} |\
cut -f 1-3,5 2>> {log} 1> {output.GC_bed}
        """

    #############################################################################
    ##### Generate IGVjs index HTML                                         #####
    #############################################################################

# Dit is nu broken door de 403 errors, ik kan het niet debuggen. Wachten tot Robert het heeft opgelost.
### Daarna de DoC en GC% toevoegen aan de html
rule Generate_IGVjs_html_file:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        fastaFai="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.fai" % config["scaffold_minLen_filter"]["minlen"],
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        bam_bai="data/scaffolds_filtered/{sample}_sorted.bam.bai",
        zipped_vcf="data/scaffolds_filtered/{sample}_filtered.vcf.gz",
        zipped_vcf_index="data/scaffolds_filtered/{sample}_filtered.vcf.gz.tbi",
        zipped_gff3="data/scaffolds_filtered/{sample}_annotation.gff.gz",
        zipped_gff3_index="data/scaffolds_filtered/{sample}_annotation.gff.gz.tbi",
        GC_bed="data/scaffolds_filtered/{sample}_GC.bedgraph",
    output:
        "data/scaffolds_filtered/{sample}_IGVjs.html"
    conda:
        "envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/Generate_IGVjs_html_file_{sample}.txt"
    threads: 1
    log:
        "logs/Generate_IGVjs_html_file_{sample}.log"
    shell:
        """
bash bin/generate_html_template_igv.sh {wildcards.sample} \
{input.fasta} \
{input.bam} \
{input.bam_bai} \
{output} \
Jovian \
{input.zipped_vcf} \
{input.zipped_vcf_index} \
{input.zipped_gff3} \
{input.zipped_gff3_index} \
{input.GC_bed} > {log} 2>&1
        """

    #############################################################################
    ##### MultiQC report of pipeline metrics                                #####
    #############################################################################

rule MultiQC_report:
    input:
        expand("data/FastQC_pretrim/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = "R1 R2".split()),
        expand("data/FastQC_posttrim/{sample}_{read}_fastqc.zip", sample = SAMPLES, read = "pR1 pR2 uR1 uR2".split()),
        expand("data/scaffolds_filtered/{sample}_insert_size_metrics.txt", sample = SAMPLES),
        expand("logs/Clean_the_data_{sample}.log", sample = SAMPLES),
        expand("logs/HuGo_removal_pt1_alignment_{sample}.log", sample = SAMPLES),
    output:
        "results/multiqc.html",
        expand("results/multiqc_data/multiqc_{program}.txt", program = ['trimmomatic','bowtie2','fastqc']),
    conda:
        "envs/MultiQC_report.yaml"
    benchmark:
        "logs/benchmark/MultiQC_report.txt"
    threads: 1
    params:
        config_file="files/multiqc_config.yaml"
    log:
        "logs/MultiQC_report.log"
    shell:
        """
multiqc --force --config {params.config_file} \
-o results/ -n multiqc.html {input} > {log} 2>&1
        """

    #############################################################################
    ##### Taxonomic classification                                          #####
    #############################################################################

rule Scaffold_classification:
    input:
        "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"]
    output:
        "data/taxonomic_classification/{sample}.blastn"
    conda:
        "envs/scaffold_classification.yaml"
    benchmark:
        "logs/benchmark/Scaffold_classification_{sample}.txt"
    threads: config["threads"]["Taxonomic_classification_of_scaffolds"]
    log:
        "logs/Scaffold_classification_{sample}.log"
    params:
        outfmt="6 std qseqid sseqid staxids sscinames stitle",
        evalue=config["taxonomic_classification"]["evalue"],
        max_target_seqs=config["taxonomic_classification"]["max_target_seqs"]
    shell:
        """
blastn -task megablast \
-outfmt "{params.outfmt}" \
-query {input} \
-evalue {params.evalue} \
-max_target_seqs {params.max_target_seqs} \
-db nt \
-num_threads {threads} \
-out {output} > {log} 2>&1
        """

    #############################################################################
    ##### Send scaffolds to their respective typingtools                    #####
    #############################################################################

rule Viral_typing:
    input:
        "data/tables/{sample}_taxClassified.tsv",
    output:
        query_scaffolds_for_NoV_TT="data/virus_typing_tables/{sample}_NoV.fa",
        query_scaffolds_for_EV_TT="data/virus_typing_tables/{sample}_EV.fa",
        XML_result_NoV_TT="data/virus_typing_tables/{sample}_NoV.xml",
        XML_result_EV_TT="data/virus_typing_tables/{sample}_EV.xml",
        parsed_CSV_result_NoV_TT="data/virus_typing_tables/{sample}_NoV.csv",
        parsed_CSV_result_EV_TT="data/virus_typing_tables/{sample}_EV.csv"
    resources:
        TT_connections=1
    conda:
        "envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/Viral_typing_{sample}.txt"
    threads: 1
    params:
        NoV_TT_http=config["typingtool_http"]["NoV_TT_http"],
        EV_TT_http=config["typingtool_http"]["EV_TT_http"],
    log:
        "logs/Viral_typing_{sample}.log"
    shell: 
        """
awk -F "\\t" '$6 == "Norwalk virus" {{print ">" $2 "\\n" $24}}' < {input} 2> {log} 1> {output.query_scaffolds_for_NoV_TT}
if [ -s "{output.query_scaffolds_for_NoV_TT}" ]
then
    curl -s --data-urlencode fasta-sequence@{output.query_scaffolds_for_NoV_TT} {params.NoV_TT_http} 2>> {log} 1> {output.XML_result_NoV_TT}
    python bin/typingtool_NoV_XML_to_csv_parser.py {wildcards.sample} {output.XML_result_NoV_TT} {output.parsed_CSV_result_NoV_TT} 2>> {log}
else
    echo -e "No scaffolds with species == Norwalk Virus in sample:\t{wildcards.sample}." >> {log}
    touch {output.XML_result_NoV_TT}
    touch {output.parsed_CSV_result_NoV_TT}
fi

awk -F "\\t" '$8 == "Picornaviridae" {{print ">" $2 "\\n" $24}}' < {input} 2>> {log} 1> {output.query_scaffolds_for_EV_TT}
if [ -s "{output.query_scaffolds_for_EV_TT}" ]
then
    curl -s --data-urlencode fasta-sequence@{output.query_scaffolds_for_EV_TT} {params.EV_TT_http} 2>> {log} 1> {output.XML_result_EV_TT}
    python bin/typingtool_EV_XML_to_csv_parser.py {wildcards.sample} {output.XML_result_EV_TT} {output.parsed_CSV_result_EV_TT} 2>> {log}
else
    echo -e "No scaffolds with family == Picornaviridae in sample:\t {wildcards.sample}." >> {log}
    touch {output.XML_result_EV_TT}
    touch {output.parsed_CSV_result_EV_TT}
fi
        """

    #############################################################################
    ##### Generate interactive taxonomy plot. LCA, magnitudes and plot      #####
    #############################################################################

rule Krona_chart_and_LCA:
    input:
        classification="data/taxonomic_classification/{sample}.blastn",
        stats="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats"
    output:
        taxtab=temp("data/taxonomic_classification/{sample}.taxtab"),
        taxMagtab=temp("data/taxonomic_classification/{sample}.taxMagtab"),
    conda:
        "envs/Krona_plot.yaml"
    benchmark:
        "logs/benchmark/Krona_chart_and_LCA_{sample}.txt"
    threads: 1
    log:
        "logs/Krona_chart_and_LCA_{sample}.log"
    params:
        bitscoreDeltaLCA=config["taxonomic_classification_LCA"]["bitscoreDeltaLCA"],
        krona_tax_db=config["databases"]["Krona_taxonomy"]
    shell:  # We rm the [conda_path]/opt/krona/taxonomy folder and replace that to our specified krona_taxonomy path, updated weekly via crontab, see scripts Robert
        """
if [ ! -L $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g') ] # If symlink to Krona db does not exist...
then # Clean and make symlink to Krona db from the current Conda env (which has a unique and unpredictable hash, therefore, the which command)
    rm -rf $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
    ln -s {params.krona_tax_db} $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
fi
ktClassifyBLAST -o {output.taxtab} -t {params.bitscoreDeltaLCA} {input.classification} > {log} 2>&1
python bin/krona_magnitudes.py {output.taxtab} {input.stats} {output.taxMagtab} >> {log} 2>&1
        """

rule Krona_chart_combine:
    input:
        sorted(expand("data/taxonomic_classification/{sample}.taxMagtab", sample = set(SAMPLES))),
    output:
        "results/krona.html"
    conda:
        "envs/Krona_plot.yaml"
    benchmark:
        "logs/benchmark/Krona_chart_combine.txt"
    threads: 1
    log:
        "logs/Krona_chart_combine.log"
    params:
        krona_tax_db=config["databases"]["Krona_taxonomy"]
    shell:
        """
if [ ! -L $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g') ] # If symlink to Krona db does not exist...
then # Clean and make symlink to Krona db from the current Conda env (which has a unique and unpredictable hash, therefore, the which command)
    rm -rf $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
    ln -s {params.krona_tax_db} $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
fi
ktImportTaxonomy {input} -i -k -m 4 -o {output} > {log} 2>&1
        """

    #############################################################################
    ##### Count annotated reads and visualise as stacked bar charts         #####
    #############################################################################
        
rule quantify_output:
    input:
        fastqc = "results/multiqc_data/multiqc_fastqc.txt",
        trimmomatic = "results/multiqc_data/multiqc_trimmomatic.txt",
        hugo = expand("data/cleaned_fastq/{sample}_{suffix}.fq",
                      sample = set(SAMPLES),
                      suffix = [ "pR1", "pR2", "unpaired" ]),
        classified = "results/all_taxClassified.tsv",
        unclassified = "results/all_taxUnclassified.tsv"
    output:
        counts = "results/profile_read_counts.csv",
        percentages = "results/profile_percentages.csv",
        graph = "results/Sample_composition_graph.html"
    conda:
        "envs/heatmaps.yaml"
    benchmark:
        "logs/benchmark/quantify_output.txt"
    threads: config["threads"]["quantify_output"]
    log:
        "logs/quantify_output.log"
    script:
        "bin/quantify_profiles.py"

    #############################################################################
    ##### Make heatmaps for superkingdoms and viruses                       #####
    #############################################################################

rule draw_heatmaps:
    input:
        classified = "results/all_taxClassified.tsv",
        numbers = "results/multiqc_data/multiqc_trimmomatic.txt"
    output:
        super="results/heatmaps/Superkingdoms_heatmap.html",
        virus=expand("results/heatmaps/Virus_{rank}_heatmap.html", rank=[ "order", "family", "genus", "species" ]),
        phage=expand("results/heatmaps/Phage_{rank}_heatmap.html", rank=[ "order", "family", "genus", "species" ]),
        bact=expand("results/heatmaps/Bacteria_{rank}_heatmap.html", rank=[ "phylum", "class", "order", "family", "genus", "species" ]),
        stats="results/Taxonomic_rank_statistics.tsv",
        vir_stats="results/Virus_rank_statistics.tsv",
        phage_stats="results/Phage_rank_statistics.tsv",
        bact_stats="results/Bacteria_rank_statistics.tsv"
    conda:
        "envs/heatmaps.yaml"
    benchmark:
        "logs/benchmark/draw_heatmaps.txt"
    threads: 1
    log:
        "logs/draw_heatmaps.log"
    script:
        "bin/draw_heatmaps.py"

    #############################################################################
    ##### Data wrangling                                                    #####
    #############################################################################

rule Merge_all_metrics_into_single_tsv:
    input:
        bbtoolsFile="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats",
        kronaFile="data/taxonomic_classification/{sample}.taxtab",
        minLenFiltScaffolds="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        virusHostDB=config["databases"]["virusHostDB"],
        NCBI_new_taxdump_rankedlineage=config["databases"]["NCBI_new_taxdump_rankedlineage"],
        NCBI_new_taxdump_host=config["databases"]["NCBI_new_taxdump_host"],
    output:
        taxClassifiedTable="data/tables/{sample}_taxClassified.tsv",
        taxUnclassifiedTable="data/tables/{sample}_taxUnclassified.tsv",
        virusHostTable="data/tables/{sample}_virusHost.tsv",
    conda:
        "envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/Merge_all_metrics_into_single_tsv_{sample}.txt"
    threads: 1
    log:
        "logs/Merge_all_metrics_into_single_tsv_{sample}.log"
    shell:
        """
python bin/merge_data.py {wildcards.sample} \
{input.bbtoolsFile} \
{input.kronaFile} \
{input.minLenFiltScaffolds} \
{input.virusHostDB} \
{input.NCBI_new_taxdump_rankedlineage} \
{input.NCBI_new_taxdump_host} \
{output.taxClassifiedTable} \
{output.taxUnclassifiedTable} \
{output.virusHostTable} > {log} 2>&1
        """

rule Concat_files:
    input:
        expand("data/tables/{sample}_{extension}", sample = SAMPLES, extension = ['taxClassified.tsv','taxUnclassified.tsv','virusHost.tsv']),
    output:
        taxClassified="results/all_taxClassified.tsv",
        taxUnclassified="results/all_taxUnclassified.tsv",
        virusHost="results/all_virusHost.tsv",
    benchmark:
        "logs/benchmark/Concat_files.txt"
    threads: 1
    log:
        "logs/Concat_files.log"
    params:
        search_folder="data/tables/",
        classified_glob="*_taxClassified.tsv",
        unclassified_glob="*_taxUnclassified.tsv",
        virusHost_glob="*_virusHost.tsv",
    shell:
        """
find {params.search_folder} -type f -name "{params.classified_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + 2> {log} 1> {output.taxClassified}
find {params.search_folder} -type f -name "{params.unclassified_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + 2>> {log} 1> {output.taxUnclassified}
find {params.search_folder} -type f -name "{params.virusHost_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + 2>> {log} 1> {output.virusHost}
        """

rule Generate_index_html:
    input:
        "results/heatmaps/Superkingdoms_heatmap.html",
        expand("results/heatmaps/Virus_{rank}_heatmap.html", rank=[ "order", "family", "genus", "species" ]),
        expand("data/scaffolds_filtered/{sample}_IGVjs.html", sample = SAMPLES),
    output:
        heatmap_index="results/Heatmap_index.html",
        IGVjs_index="results/IGVjs_index.html",
    benchmark:
        "logs/benchmark/Generate_index_html.txt"
    threads: 1
    log:
        "logs/Generate_index_html.log"
    params:
        heatmap_title=config["HTML_index_titles"]["heatmap_title"],
        igvjs_title=config["HTML_index_titles"]["IGVjs_title"],
        http_adress=config["server_info"]["http_adress"],
        port=config["server_info"]["port"],
    shell:
        """
tree -H "heatmaps" -L 1 -T "{params.heatmap_title}" --noreport --charset utf-8 -P "*.html" -o {output.heatmap_index} results/heatmaps/ > {log} 2>&1 
tree -H "{params.http_adress}:{params.port}/igv/Jovian/data/scaffolds_filtered" -L 1 -T "{params.igvjs_title}" --noreport --charset utf-8 -P "*.html" -o {output.IGVjs_index} data/scaffolds_filtered/ >> {log} 2>&1
        """

rule Concat_TT_files:
    input:
        expand("data/virus_typing_tables/{sample}_{virus}.csv", sample = SAMPLES, virus = ['NoV','EV'])
    output:
        NoV="results/all_NoV-TT.tsv",
        EV="results/all_EV-TT.tsv",
    benchmark:
        "logs/benchmark/Concat_TT_files.txt"
    threads: 1
    log:
        "logs/Concat_TT_files.log"
    params:
        search_folder="data/virus_typing_tables/",
        NoV_glob="*_NoV.csv",
        EV_glob="*_EV.csv",
    shell:
        """
find {params.search_folder} -type f -name "{params.NoV_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + 2> {log} 1> {output.NoV}
find {params.search_folder} -type f -name "{params.EV_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + 2>> {log} 1> {output.EV} 
        """

rule Concat_filtered_SNPs:
    input:
        expand("data/scaffolds_filtered/{sample}_filtered.vcf", sample = SAMPLES)
    output:
        "results/all_filtered_SNPs.tsv"
    conda:
        "envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/Concat_filtered_SNPs.txt"
    threads: 1
    params:
        vcf_folder_glob="data/scaffolds_filtered/\*_filtered.vcf"
    log:
        "logs/Concat_filtered_SNPs.log"
    shell:
        """
python bin/concat_filtered_vcf.py {params.vcf_folder_glob} {output} > {log} 2>&1
        """

#################################################################################
##### These are the conditional cleanup rules                               #####
#################################################################################

onsuccess:
    shell("""
echo -e "\nCleaning up..."
echo -e "\tRemoving empty folders..."
find data -depth -type d -not \( -path data/scaffolds_raw -prune \) -empty -delete

echo -e "\tRemoving empty typing tool files..."
find data/virus_typing_tables/ -type f -empty -delete

echo -e "\tGenerating rulegraph plot of the pipeline..."
snakemake --unlock --config sample_sheet=sample_sheet.yaml
snakemake --rulegraph --config sample_sheet=sample_sheet.yaml | dot -Tpng > results/rulegraph_Jovian.png

echo -e "\tGenerating Snakemake report..."
snakemake --report results/snakemake_report.html --config sample_sheet=sample_sheet.yaml
         """)
    
# perORFcoverage output file van de bbtools scaffold metrics nog importeren in data wrangling part!