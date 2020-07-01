


rule MultiQC_report:
    input: 
        #expand("{l}Cleanup_{sample}.log", l = logdir, sample = SAMPLES),
        expand( "{p}{sample}.fastp.json",
                p       =   f"{datadir + cln + json}",
                sample  =   SAMPLES
                ),
        expand( "{l}Primer_removal_{sample}.log",
                l       =   logdir,
                sample  =   SAMPLES
                )
    output: 
        f"{res}multiqc.html"
    conda:
        f"{conda_envs}MultiQC_report.yaml"
    log:
        f"{logdir}MultiQC_report.log"
    benchmark:
        f"{logdir + bench}MultiQC_report.txt"
    threads: 1
    params:
        config_file = "files/multiqc_config.yaml",
        output_dir = res,
    shell:
        """
multiqc -d --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
        """