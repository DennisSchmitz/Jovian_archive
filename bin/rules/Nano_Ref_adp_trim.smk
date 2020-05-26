

rule Adapter_trimming:
    input: 
        lambda wildcards: SAMPLES[wildcards.sample]
    output: 
        trimmeddata = datadir + cln + trims + "{sample}.fastq"
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    benchmark:
        logdir + bench + "Adapter_trimming_{sample}.txt"
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