
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
        expand( "logs/Clean_the_data_{sample}.log",
                sample  =   SAMPLES
                ),
        expand( "logs/HuGo_removal_pt1_alignment_{sample}.log",
                sample  =   SAMPLES
                ),
    output:
        "results/multiqc.html",
        expand( "results/multiqc_data/multiqc_{program}.txt",
                program =   [   'trimmomatic',
                                'bowtie2',
                                'fastqc'
                            ]
                ),
    conda:
        conda_envs + "MultiQC_report.yaml"
    log:
        "logs/MultiQC_report.log"
    benchmark:
        "logs/benchmark/MultiQC_report.txt"
    threads: 1
    params:
        config_file =   "files/multiqc_config.yaml"
    shell:
        """
multiqc --force --config {params.config_file} \
-o results/ -n multiqc.html {input} > {log} 2>&1
        """