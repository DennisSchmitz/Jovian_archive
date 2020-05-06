
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Scaffold_classification:
    input:
        "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"]
    output:
        "data/taxonomic_classification/{sample}.blastn"
    conda:
        "../envs/scaffold_classification.yaml"
    benchmark:
        "logs/benchmark/Scaffold_classification_{sample}.txt"
    threads: config["threads"]["Taxonomic_classification_of_scaffolds"]
    log:
        "logs/Scaffold_classification_{sample}.log"
    params:
        outfmt="6 std qseqid sseqid staxids sscinames stitle",
        evalue=config["taxonomic_classification"]["evalue"],
        max_target_seqs=config["taxonomic_classification"]["max_target_seqs"],
        max_hsps=config["taxonomic_classification"]["max_hsps"]
    shell:
        """
blastn -task megablast \
-outfmt "{params.outfmt}" \
-query {input} \
-evalue {params.evalue} \
-max_target_seqs {params.max_target_seqs} \
-max_hsps {params.max_hsps} \
-db nt \
-num_threads {threads} \
-out {output} > {log} 2>&1
        """