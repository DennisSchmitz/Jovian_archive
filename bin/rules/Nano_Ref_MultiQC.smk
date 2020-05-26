


rule MultiQC_report:
    input: 
        #expand("{l}Cleanup_{sample}.log", l = logdir, sample = SAMPLES),
        expand("{path}{sample}.fastp.json", path = datadir + cln + json, sample = SAMPLES),
        expand("{l}Primer_removal_{sample}.log", l = logdir, sample = SAMPLES),
    output: 
        f"{res}multiqc.html"
        #expand("{out}{program}_multiqc.txt", out = res + mqc, program = ['fastp','cutadapt']),
    conda:
        conda_envs + "MultiQC_report.yaml"
    log:
        logdir + "MultiQC_report.log"
    benchmark:
        logdir + bench + "MultiQC_report.txt"
    threads: 1
    params:
        config_file = "files/multiqc_config.yaml",
        output_dir = res,
    shell:
        """
multiqc -d --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
        """