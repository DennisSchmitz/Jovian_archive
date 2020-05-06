
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule count_mapped_reads:
    input:
        "data/scaffolds_filtered/{sample}_sorted.bam"
        #expand("data/scaffolds_filtered/{sample}_sorted.bam", sample = SAMPLES)
    output:
        "results/counts/Mapped_read_counts-{sample}.tsv"
    conda:
        "../envs/scaffold_analyses.yaml"
    benchmark:
        "logs/benchmark/count_mapped_reads-{sample}.txt"
    threads: 1
    log:
        "logs/count_mapped_reads-{sample}.txt"
    shell:
        """
bash bin/scripts/count_mapped_reads.sh {input} > {output} 2> {log}
        """