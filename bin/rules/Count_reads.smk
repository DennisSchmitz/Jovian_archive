
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule count_mapped_reads:
    input:
        rules.Read2scaffold_alignment_with_rmDup_and_fraglength.output.bam
    output:
        f"{res + cnt}" + "Mapped_read_counts-{sample}.tsv"
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}" + "count_mapped_reads-{sample}.txt"
    benchmark:
        f"{logdir + bench}" + "count_mapped_reads-{sample}.txt"
    threads: 1
    resources:
        memory = 12
    shell:
        """
bash bin/scripts/count_mapped_reads.sh {input} > {output} 2> {log}
        """