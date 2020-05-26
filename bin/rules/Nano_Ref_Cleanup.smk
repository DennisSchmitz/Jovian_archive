



rule Cleanup:
    input:
        fastq   =   rules.Cut_primers.output.cleaneddata_pt1
    output:
        qc_fastq    =   f"{datadir + cln + datadir}" + "{sample}.fastq",
        qc_html     =   f"{datadir + cln + html}" + "{sample}.html",
        qc_json     =   f"{datadir + cln + json}" + "{sample}.fastp.json"
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "Data_Cleanup_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Data_Cleanup_{sample}.txt"
    threads: 26
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