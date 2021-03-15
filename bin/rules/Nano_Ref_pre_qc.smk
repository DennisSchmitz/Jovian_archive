rule raw_quality_control:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]
    output:
        html    =   f"{datadir + qc_pre}" + "{sample}_fastqc.html",
        zip     =   f"{datadir + qc_pre}" + "{sample}_fastqc.zip",
    conda:
        f"{conda_envs}Nano_clean.yaml"
    log:
        f"{logdir}" + "raw_quality_control_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "raw_quality_control_{sample}.txt"
    threads: config["threads"]["Nanopore_QC"]
    resources:
        memory = config["threads"]["Nanopore_QC"] * 4
    params:
        output_dir  =   f"{datadir + qc_pre}"
    shell:
        """
bash bin/scripts/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log}
        """