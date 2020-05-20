
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
        "results/counts/Mapped_read_counts.tsv"
    conda:
        conda_envs + "data_wrangling.yaml"
    benchmark:
        "logs/benchmark/concatenate_read_counts.txt"
    threads: 1
    log:
        "logs/concatenate_read_counts.txt"
    shell:
        """
        bin/scripts/concatenate_mapped_read_counts.py \
        -i {input} \
        -o {output} \
        > {log} 2>&1
        """