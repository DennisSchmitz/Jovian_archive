rule cleaned_quality_control:
    input: 
        rules.Cleanup.output.qc_fastq
    output:
        html    =   f"{datadir + qc_post}" + "{sample}_fastqc.html",
        zip     =   f"{datadir + qc_post}" + "{sample}_fastqc.zip"
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "cleaned_quality_control_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "cleaned_quality_control_{sample}.txt"
    threads: 1
    params:
        outdir  =   f"{datadir + qc_post}"
    shell:
        """
if [ -s "{input}" ] # If file exists and is NOT empty (i.e. filesize > 0) do...
then
    fastqc --quiet --outdir {params.outdir} {input} > {log} 2>&1
else
    touch {output.html}
    touch {output.zip}
fi
        """

