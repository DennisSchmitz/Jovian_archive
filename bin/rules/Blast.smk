
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Scaffold_classification:
    input:
        rules.De_novo_assembly.output.filt_scaffolds
    output:
        f"{datadir + taxclas}" + "{sample}.blastn"
    conda:
        f"{conda_envs}scaffold_classification.yaml"
    log:
        f"{logdir}" + "Scaffold_classification_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "logs/benchmark/Scaffold_classification_{sample}.txt"
    threads: config["threads"]["Classification_of_scaffolds"]
    resources: 
        memory = config["threads"]["Classification_of_scaffolds"] * 12
    params:
        outfmt          =   "6 std qseqid sseqid staxids sscinames stitle",
        evalue          =   config["Illumina_meta"]["Classification"]["e_value"],
        max_target_seqs =   config["Illumina_meta"]["Classification"]["max_target_seqs"],
        max_hsps        =   config["Illumina_meta"]["Classification"]["max_hsps"]
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