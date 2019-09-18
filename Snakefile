"""
Authors: Dennis Schmitz, Sam Nooij, Robert Verhagen, Thierry Janssens, Jeroen Cremer, Florian Zwagemaker, Mark Kroon, Erwin van Wieringen, Annelies Kroneman, Harry Vennema, Marion Koopmans
Organisation: Rijksinstituut voor Volksgezondheid en Milieu (RIVM)
Department: Virology - Emerging and Endemic Viruses (EEV)
Date: 23-08-2018
Changelog, examples, installation guide and explanation on:
   https://github.com/DennisSchmitz/Jovian
"""

#################################################################################
##### Import config file, sample_sheet and set output folder names          #####
#################################################################################

shell.executable("/bin/bash")

configfile: "profile/pipeline_parameters.yaml"
configfile: "profile/variables.yaml"

import pprint
import os
import yaml
yaml.warnings({'YAMLLoadWarning': False}) # Suppress yaml "unsafe" warnings.

SAMPLES = {}
with open(config["sample_sheet"]) as sample_sheet_file:
    SAMPLES = yaml.load(sample_sheet_file) # SAMPLES is a dict with sample in the form sample > read number > file. E.g.: SAMPLES["sample_1"]["R1"] = "x_R1.gz"

#################################################################################
##### The `onstart` checker codeblock                                       #####
#################################################################################

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
        echo -e "==> User-specified background reference (default: Homo Sapiens NCBI GRch38 NO DECOY genome): <==\n$(ls -lah $(grep "    background_ref:" profile/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" > results/log_db.txt
        echo -e "\n==> Virus-Host Interaction Database: <==\n$(ls -lah $(grep "    virusHostDB:" profile/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" >> results/log_db.txt
        echo -e "\n==> Krona Taxonomy Database: <==\n$(ls -lah $(grep "    Krona_taxonomy:" profile/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" >> results/log_db.txt
        echo -e "\n==> NCBI new_taxdump Database: <==\n$(ls -lah $(grep "    NCBI_new_taxdump_rankedlineage:" profile/pipeline_parameters.yaml | cut -f 2 -d ":") $(grep "    NCBI_new_taxdump_host:" profile/pipeline_parameters.yaml | cut -f 2 -d ":"))\n" >> results/log_db.txt
        echo -e "\n==> NCBI Databases as specified in ~/.ncbirc: <==\n$(ls -lah $(grep "BLASTDB=" ~/.ncbirc | cut -f 2 -d "=" | tr "::" " "))\n" >> results/log_db.txt
        
        echo -e "\tGenerating config file log..."
        rm -f results/log_config.txt
        for file in profile/*.yaml
        do
            echo -e "\n==> Contents of file \"${{file}}\": <==" >> results/log_config.txt
            cat ${{file}} >> results/log_config.txt
            echo -e "\n\n" >> results/log_config.txt
        done
    """)

#    shell("\n\t\tPlaceholder for Jupyter config checker.\n\t\tPlaceholder Jupyter Notebook theme settings\n'")

#################################################################################
##### Specify Jovian's final output:                                        #####
#################################################################################

localrules: 
    all,
    Generate_index_html,
    Generate_IGVjs_html_file,
    quantify_output,
    Concat_files,
    Concat_filtered_SNPs,

rule all:
    input:
        expand("data/cleaned_fastq/{sample}_{read}.fq", sample = SAMPLES, read = [ 'pR1', 'pR2', 'unpaired' ]), # Extract unmapped & paired reads AND unpaired from HuGo alignment; i.e. cleaned fastqs
        expand("data/scaffolds_raw/{sample}/scaffolds.fasta", sample = SAMPLES), # SPAdes assembly output
        expand("data/scaffolds_filtered/{sample}_scaffolds_ge{len}nt.{extension}", sample = SAMPLES, len = config["scaffold_minLen_filter"]["minlen"], extension = [ 'fasta', 'fasta.fai' ]), # Filtered SPAdes Scaffolds
        expand("data/scaffolds_filtered/{sample}_sorted.bam", sample = SAMPLES), # BWA mem alignment for fragment length analysis
        expand("data/scaffolds_filtered/{sample}_{extension}", sample = SAMPLES, extension = [ 'ORF_AA.fa', 'ORF_NT.fa', 'annotation.gff', 'annotation.gff.gz', 'annotation.gff.gz.tbi', 'contig_ORF_count_list.txt' ]), # Prodigal ORF prediction output
        expand("data/scaffolds_filtered/{sample}_{extension}", sample = SAMPLES, extension = [ 'unfiltered.vcf', 'filtered.vcf', 'filtered.vcf.gz', 'filtered.vcf.gz.tbi' ]), # SNP calling output
        expand("data/scaffolds_filtered/{sample}_GC.bedgraph", sample = SAMPLES), # Percentage GC content per specified window
        expand("data/scaffolds_filtered/{sample}_IGVjs.html", sample = SAMPLES), # IGVjs html's
        expand("data/taxonomic_classification/{sample}.blastn", sample = SAMPLES), # MegablastN output
        expand("data/tables/{sample}_{extension}", sample = SAMPLES, extension = [ 'taxClassified.tsv', 'taxUnclassified.tsv', 'virusHost.tsv' ]), # Tab seperated tables with merged data
        expand("results/{file}", file = [ 'all_taxClassified.tsv', 'all_taxUnclassified.tsv', 'all_virusHost.tsv', 'all_filtered_SNPs.tsv' ]), # Concatenated classification, virus host and typing tool tables
        expand("results/{file}", file = [ 'heatmaps/Superkingdoms_heatmap.html', 'Sample_composition_graph.html', 'Taxonomic_rank_statistics.tsv', 'Virus_rank_statistics.tsv', 'Phage_rank_statistics.tsv', 'Bacteria_rank_statistics.tsv' ]), # Taxonomic profile and heatmap output
        "results/heatmaps/Virus_heatmap.html", # Virus (excl. phages) order|family|genus|species level heatmap for the entire run
        "results/heatmaps/Phage_heatmap.html", # Phage order|family|genus|species heatmaps for the entire run (based on a selection of phage families)
        "results/heatmaps/Bacteria_heatmap.html", # Bacteria phylum|class|order|family|genus|species level heatmap for the entire run
        expand("results/{file}.html", file = [ 'multiqc', 'krona', 'IGVjs_index' ]), # Reports and IGVjs index.html

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
        html="data/FastQC_pretrim/{sample}_{read}_fastqc.html",
        zip="data/FastQC_pretrim/{sample}_{read}_fastqc.zip"
    conda:
        "envs/QC_and_clean.yaml"
    benchmark:
        "logs/benchmark/QC_raw_data_{sample}_{read}.txt"
    threads: 1
    log:
        "logs/QC_raw_data_{sample}_{read}.log"
    params:
        output_dir="data/FastQC_pretrim/"
    shell:
        """
bash bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log}
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
        html="data/FastQC_posttrim/{sample}_{read}_fastqc.html",
        zip="data/FastQC_posttrim/{sample}_{read}_fastqc.zip"
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
        background_ref=config["databases"]["background_ref"],
        r1="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR1.fastq",
        r2="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR2.fastq",
        r1_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR1.fastq",
        r2_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR2.fastq",
    output:
        sorted_bam="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam",
        sorted_bam_index="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai",
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
-x {input.background_ref} \
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
        contig_ORF_count_list="data/scaffolds_filtered/{sample}_contig_ORF_count_list.txt"
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

egrep "^>" {output.ORF_NT_fasta} | sed 's/_/ /6' | tr -d ">" | cut -f 1 -d " " | uniq -c > {output.contig_ORF_count_list}
        """

rule Generate_contigs_metrics:
    input:
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        ORF_NT_fasta="data/scaffolds_filtered/{sample}_ORF_NT.fa",
    output:
        summary="data/scaffolds_filtered/{sample}_MinLenFiltSummary.stats",
        perScaffold="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats",
        perORFcoverage="data/scaffolds_filtered/{sample}_perORFcoverage.stats",
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
out={output.perScaffold} \
secondary=f \
samstreamer=t \
2> {output.summary} 1> {log}
        """

rule Determine_GC_content:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        fasta_fai="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.fai" % config["scaffold_minLen_filter"]["minlen"],
    output:
        fasta_sizes="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.sizes" % config["scaffold_minLen_filter"]["minlen"],
        bed_windows="data/scaffolds_filtered/{sample}.windows",
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
        max_target_seqs=config["taxonomic_classification"]["max_target_seqs"],
        max_hsps=config["taxonomic_classification"]["max_hsps"]
    shell:
        """
blastn -task megablast \
-outfmt "{params.outfmt}" \
-query {input} \
-evalue {params.evalue} \
-max_target_seqs {params.max_target_seqs} \
-max_hsps {params.max_hsps} \
-db nt \
-num_threads {threads} \
-out {output} > {log} 2>&1
        """

if config["taxonomic_classification_LCA"]["Krona"] == True:
    include: "rules/Krona_LCA.smk"

if config["taxonomic_classification_LCA"]["mgkit"] == True:
    include: "rules/mgkit_LCA.smk"

    #############################################################################
    ##### Send scaffolds to their respective typingtools                    #####
    #############################################################################

# See issue #29 (https://github.com/DennisSchmitz/Jovian/issues/29)
#    Will be added to the pipeline again at a later date.
#    For now, these analysis can be done via the scripts attached along with Jovian when needed. For instructions, see the link above.

    #############################################################################
    ##### Generate interactive taxonomy plot. LCA, magnitudes and plot      #####
    #############################################################################

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
        
rule count_mapped_reads:
    input:
        expand("data/scaffolds_filtered/{sample}_sorted.bam", sample = SAMPLES)
    output:
        "results/Mapped_read_counts.tsv"
    conda:
        "envs/scaffold_analyses.yaml"
    benchmark:
        "logs/benchmark/count_mapped_reads.txt"
    threads: 1
    log:
        "logs/count_mapped_reads.txt"
    params:
        input_dir="data/scaffolds_filtered/"
    shell:
        """
bash bin/count_mapped_reads.sh {params.input_dir} > {output} 2> {log}
        """

rule quantify_output:
    input:
        fastqc = "results/multiqc_data/multiqc_fastqc.txt",
        trimmomatic = "results/multiqc_data/multiqc_trimmomatic.txt",
        hugo = expand("data/cleaned_fastq/{sample}_{suffix}.fq",
                      sample = set(SAMPLES),
                      suffix = [ "pR1", "pR2", "unpaired" ]),
        classified = "results/all_taxClassified.tsv",
        unclassified = "results/all_taxUnclassified.tsv",
        mapped_reads = "results/Mapped_read_counts.tsv"
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
    shell:
        """
python bin/quantify_profiles.py \
-f {input.fastqc} \
-t {input.trimmomatic} \
-hg {input.hugo} \
-c {input.classified} \
-u {input.unclassified} \
-m {input.mapped_reads} \
-co {output.counts} \
-p {output.percentages} \
-g {output.graph} \
-cpu {threads} \
-l {log}
        """

    #############################################################################
    ##### Make heatmaps for superkingdoms and viruses                       #####
    #############################################################################

rule draw_heatmaps:
    input:
        classified = "results/all_taxClassified.tsv",
        numbers = "results/multiqc_data/multiqc_trimmomatic.txt"
    output:
        super_quantities="results/Superkingdoms_quantities_per_sample.csv",
        super="results/heatmaps/Superkingdoms_heatmap.html",
        virus="results/heatmaps/Virus_heatmap.html",
        phage="results/heatmaps/Phage_heatmap.html",
        bact="results/heatmaps/Bacteria_heatmap.html",
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
    shell:
        """
python bin/draw_heatmaps.py -c {input.classified} -n {input.numbers} -sq {output.super_quantities} -st {output.stats} -vs {output.vir_stats} -ps {output.phage_stats} -bs {output.bact_stats} -s {output.super} -v {output.virus} -p {output.phage} -b {output.bact} > {log} 2>&1
        """

    #############################################################################
    ##### Data wrangling                                                    #####
    #############################################################################

rule Merge_all_metrics_into_single_tsv:
    input:
        bbtoolsFile="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats",
        kronaFile="data/taxonomic_classification/{sample}.taxtab",
        minLenFiltScaffolds="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        scaffoldORFcounts="data/scaffolds_filtered/{sample}_contig_ORF_count_list.txt",
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
{input.scaffoldORFcounts} \
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
        expand("data/scaffolds_filtered/{sample}_IGVjs.html", sample = SAMPLES),
    output:
        IGVjs_index="results/IGVjs_index.html",
    benchmark:
        "logs/benchmark/Generate_index_html.txt"
    threads: 1
    log:
        "logs/Generate_index_html.log"
    params:
        igvjs_title=config["HTML_index_titles"]["IGVjs_title"],
        http_adress=config["Server_host"]["hostname"],
        port=config["server_info"]["port"],
    shell:
        """
bin/create_igv_index.sh
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

onerror:
    shell("""
        rm -f data/scaffolds_filtered/*.html
        rm -f results/IGVjs_index.html
    """)


onsuccess:
    shell("""
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
        else
            echo -e "\t\tYou chose to not remove temp files: the human genome alignment files are not removed."
        fi

        echo -e "\tCreating symlinks for the interactive genome viewer..."
        bin/set_symlink.sh

        echo -e "\tGenerating HTML index of log files..."
        tree -hD --dirsfirst -H "../logs" -L 2 -T "Logs overview" --noreport --charset utf-8 -P "*" -o results/logfiles_index.html logs/

        echo -e "\tGenerating Snakemake report..."
        snakemake --unlock --config sample_sheet=sample_sheet.yaml
        snakemake --report results/snakemake_report.html --config sample_sheet=sample_sheet.yaml

        echo -e "Finished"
    """)
    
# perORFcoverage output file van de bbtools scaffold metrics nog importeren in data wrangling part!
