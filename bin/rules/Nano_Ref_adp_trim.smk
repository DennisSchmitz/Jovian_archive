

rule Adapter_trimming:
    input: 
        lambda wildcards: SAMPLES[wildcards.sample]
    output: 
        trimmeddata =   f"{datadir + cln + trims}" + "{sample}.fastq"
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "Adapter_trimming_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Adapter_trimming_{sample}.txt"
    threads: config["threads"]["Adapter_trimming_cpu"]
    shell:
        """
porechop \
-i {input} \
-o {output.trimmeddata} \
--threads {threads} > {log} 2>&1
        """