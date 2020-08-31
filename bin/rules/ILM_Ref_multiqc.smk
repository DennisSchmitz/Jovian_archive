##########################!
# Gejat uit Jovian core met minor changes, kunnen we waarschijlijk efficienter doen.
# TODO the report is still a bit dirty since we include two bowtie2 metric files:
#### TODO one for the hugo removal
#### TODO another for the ref alignment
#### TODOD hence the '-d' flag in the multiqc command based on https://multiqc.info/docs/#directory-names
rule Illumina_MultiQC_report:
    input:
        expand( "{p}{sample}_{read}_fastqc.zip",
                p       =   f"{datadir + qc_pre}",
                sample  =   SAMPLES,
                read    =   "R1 R2".split()
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand( "{p}{sample}_{read}_fastqc.zip",
                p       =   f"{datadir + qc_post}",
                sample  =   SAMPLES,
                read    =   "pR1 pR2 uR1 uR2".split()
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand( "{p}Clean_the_data_{sample}.log",
                p       = f"{logdir}",
                sample  = SAMPLES
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand( "{p}HuGo_removal_pt1_alignment_{sample}.log",
                p       = f"{logdir}",
                sample  = SAMPLES
                ), # TODO dit moet nog verbetert worden qua smk syntax
        expand( "{p}Illumina_align_to_reference_it2_{sample}.log",
                p       = f"{logdir}",
                sample  = SAMPLES
                ) # TODO dit moet nog verbetert worden qua smk syntax
    output:
        f"{res}multiqc.html",
        expand( "{p}multiqc_{program}.txt",
                p       =   f"{res + mqc_data}",
                program =   [   'trimmomatic',
                                'bowtie2',
                                'fastqc'
                                ]
                )
    conda:
        f"{conda_envs}MultiQC_report.yaml"
    log:
        f"{logdir}Illumina_MultiQC_report.log"
    benchmark:
        f"{logdir + bench}Illumina_MultiQC_report.txt"
    threads: 1
    params:
        config_file =   f"{fls}multiqc_config.yaml",
        output_dir  =   f"{res}",
    shell:
        """
multiqc -d --force --config {params.config_file} \
-o {params.output_dir} -n multiqc.html {input} > {log} 2>&1
        """