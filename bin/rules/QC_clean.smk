
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule QC_clean_data:
    input:
        "data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_{read}.fastq"
    output:
        html="data/FastQC_posttrim/{sample}_{read}_fastqc.html",
        zip="data/FastQC_posttrim/{sample}_{read}_fastqc.zip"
    conda:
        "../envs/QC_and_clean.yaml"
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