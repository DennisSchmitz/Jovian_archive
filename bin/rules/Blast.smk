
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Scaffold_classification:
    input:
        "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["Illumina_meta"]["minlen"]
    output:
        "data/taxonomic_classification/{sample}.blastn"
    conda:
        "../envs/scaffold_classification.yaml"
    benchmark:
        "logs/benchmark/Scaffold_classification_{sample}.txt"
    threads: config["threads"]["Classification_of_scaffolds"]
    log:
        "logs/Scaffold_classification_{sample}.log"
    params:
        outfmt = "6 std qseqid sseqid staxids sscinames stitle",
        evalue = config["Illumina_meta"]["Classification"]["e_value"],
        max_target_seqs = config["Illumina_meta"]["Classification"]["max_target_seqs"],
        max_hsps = config["Illumina_meta"]["Classification"]["max_hsps"]
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