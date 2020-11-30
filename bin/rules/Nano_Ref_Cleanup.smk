rule Cleanup:
    input:
        fastq   =   rules.Remove_Adapters_pt2.output
    output:
        qc_fastq    =   f"{datadir + cln + datadir}" + "{sample}.fastq",
        qc_html     =   f"{datadir + cln + html}" + "{sample}.fastp.html",
        qc_json     =   f"{datadir + cln + json}" + "{sample}.fastp.json"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    log:
        f"{logdir}" + "Data_Cleanup_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Data_Cleanup_{sample}.txt"
    threads: config["threads"]["Nanopore_cleanup"]
    params:
        QualityFilter   =   config["Nanopore_ref"]["Quality_score"]
    shell:
        """
fastp \
-i {input.fastq} \
-q {params.QualityFilter} \
-o {output.qc_fastq} \
-h {output.qc_html} \
-j {output.qc_json} > {log} 2>&1
        """