
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule MultiQC_report:
    input:
        expand( rules.QC_raw_data.output.zip,
                sample  =   SAMPLES,
                read    =   "R1 R2".split()
                ),
        expand( rules.QC_clean_data.output.zip,
                sample  =   SAMPLES,
                read    =   "pR1 pR2 uR1 uR2".split()
                ),
        expand( rules.Read2scaffold_alignment_with_rmDup_and_fraglength.output.txt,
                sample  =   SAMPLES
                ),
        expand( "{p}Clean_the_data_{sample}.log",
                p       =   f"{logdir}",
                sample  =   SAMPLES
                ),
        expand( "{p}HuGo_removal_pt1_alignment_{sample}.log",
                p       =   f"{logdir}",
                sample  =   SAMPLES
                ),
    output:
        f"{res}multiqc.html",
        expand( "{p}multiqc_{program}.txt",
                p       =   f"{res + mqc_data}",
                program =   [   'trimmomatic',
                                'bowtie2',
                                'fastqc'
                            ]
                ),
    conda:
        f"{conda_envs}MultiQC_report.yaml"
    log:
        f"{logdir}MultiQC_report.log"
    benchmark:
        f"{logdir + bench}MultiQC_report.txt"
    threads: 1
    params:
        config_file =   f"{fls}multiqc_config.yaml"
    shell:
        """
multiqc --force --config {params.config_file} \
-o results/ -n multiqc.html {input} > {log} 2>&1
        """