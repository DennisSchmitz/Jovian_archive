
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule quantify_output:
    input:
        fastqc = "results/multiqc_data/multiqc_fastqc.txt",
        trimmomatic = "results/multiqc_data/multiqc_trimmomatic.txt",
        hugo = expand("data/cleaned_fastq/{sample}_{suffix}.fq",
                      sample = set(SAMPLES),
                      suffix = [ "pR1", "pR2", "unpaired" ]),
        classified = "results/all_taxClassified.tsv",
        unclassified = "results/all_taxUnclassified.tsv",
        mapped_reads = "results/counts/Mapped_read_counts.tsv"
    output:
        counts = "results/profile_read_counts.csv",
        percentages = "results/profile_percentages.csv",
        graph = "results/Sample_composition_graph.html"
    conda:
        "../envs/heatmaps.yaml"
    benchmark:
        "logs/benchmark/quantify_output.txt"
    threads: config["threads"]["quantify_output"]
    log:
        "logs/quantify_output.log"
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