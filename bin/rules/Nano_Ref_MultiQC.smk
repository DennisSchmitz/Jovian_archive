


rule MultiQC_report:
    input: 
        expand( rules.raw_quality_control.output.zip,
                sample  =   SAMPLES
                ),
        expand( rules.cleaned_quality_control.output.zip,
                sample  =   SAMPLES
                ),
        expand( rules.Cut_primers.log,
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
    resources:
        memory = 12
    params:
        config_file = "files/multiqc_config.yaml",
        output_dir = res,
    shell:
        """
multiqc -d --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
        """