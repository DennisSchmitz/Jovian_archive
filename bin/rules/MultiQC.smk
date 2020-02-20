
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
        "../envs/MultiQC_report.yaml"
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