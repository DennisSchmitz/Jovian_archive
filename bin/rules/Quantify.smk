
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule quantify_output:
    input:
        classified      =   rules.Concat_files.output.taxClassified,
        unclassified    =   rules.Concat_files.output.taxUnclassified,
        mapped_reads    =   rules.concatenate_read_counts.output,
        fastqc          =   f"{res + mqc_data}multiqc_fastqc.txt",
        trimmomatic     =   f"{res + mqc_data}multiqc_trimmomatic.txt",
        hugo            =   expand( "data/cleaned_fastq/{sample}_{suffix}.fq",
                                    p       =   f"{datadir + cln}",
                                    sample  =   set(SAMPLES),
                                    suffix  =   [   "pR1",
                                                    "pR2",
                                                    "unpaired"
                                                    ]
                                    )
    output:
        counts      =   f"{res}profile_read_counts.csv",
        percentages =   f"{res}profile_percentages.csv",
        graph       =   f"{res}Sample_composition_graph.html"
    conda:
        f"{conda_envs}heatmaps.yaml"
    log:
        f"{logdir}quantify_output.log"
    benchmark:
        f"{logdir + bench}quantify_output.txt"
    threads: config["threads"]["quantify_output"]
    resources:
        memory = config["threads"]["quantify_output"] * 4
    shell:
        """
python bin/scripts/quantify_profiles.py \
-f {input.fastqc} \
-t {input.trimmomatic} \
-hg {input.hugo} \
-c {input.classified} \
-u {input.unclassified} \
-m {input.mapped_reads} \
-co {output.counts} \
-p {output.percentages} \
-g {output.graph} \
-cpu {threads} \
-l {log}
        """