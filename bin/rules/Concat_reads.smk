
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule concatenate_read_counts:
    input:
        expand( rules.count_mapped_reads.output,
                sample = SAMPLES
                )
    output:
        f"{res + cnt}Mapped_read_counts.tsv"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}concatenate_read_counts.txt"
    benchmark:
        f"{logdir + bench}concatenate_read_counts.txt"
    threads: 1
    resources:
        memory = 4
    shell:
        """
        bin/scripts/concatenate_mapped_read_counts.py \
        -i {input} \
        -o {output} \
        > {log} 2>&1
        """