rule raw_quality_control:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        html    =   f"{datadir + qc_pre}" + "{sample}_fastqc.html",
        zip     =   f"{datadir + qc_pre}" + "{sample}_fastqc.zip",
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "FASTQC_PRECLEAN_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "FASTQC_PRECLEAN_{sample}.txt"
    threads: 1
    params:
        output_dir  =   f"{datadir + qc_pre}"
    shell:
        """
bash bin/scripts/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log}
        """