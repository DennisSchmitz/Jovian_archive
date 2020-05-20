
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule QC_raw_data:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][wildcards.read]
    output:
        html    =   "data/FastQC_pretrim/{sample}_{read}_fastqc.html",
        zip     =   "data/FastQC_pretrim/{sample}_{read}_fastqc.zip"
    conda:
        conda_envs + "QC_and_clean.yaml"
    benchmark:
        "logs/benchmark/QC_raw_data_{sample}_{read}.txt"
    threads: 1
    log:
        "logs/QC_raw_data_{sample}_{read}.log"
    params:
        output_dir  =   "data/FastQC_pretrim/"
    shell:
        """
bash bin/scripts/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log}
        """