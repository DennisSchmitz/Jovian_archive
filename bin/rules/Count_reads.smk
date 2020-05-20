
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule count_mapped_reads:
    input:
        rules.Read2scaffold_alignment_with_rmDup_and_fraglength.output.bam
    output:
        "results/counts/Mapped_read_counts-{sample}.tsv"
    conda:
        conda_envs + "Sequence_analysis.yaml"
    benchmark:
        "logs/benchmark/count_mapped_reads-{sample}.txt"
    threads: 1
    log:
        "logs/count_mapped_reads-{sample}.txt"
    shell:
        """
bash bin/scripts/count_mapped_reads.sh {input} > {output} 2> {log}
        """