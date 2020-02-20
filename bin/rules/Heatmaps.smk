
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule draw_heatmaps:
    input:
        classified = "results/all_taxClassified.tsv",
        numbers = "results/multiqc_data/multiqc_trimmomatic.txt"
    output:
        super_quantities="results/Superkingdoms_quantities_per_sample.csv",
        super="results/heatmaps/Superkingdoms_heatmap.html",
        virus="results/heatmaps/Virus_heatmap.html",
        phage="results/heatmaps/Phage_heatmap.html",
        bact="results/heatmaps/Bacteria_heatmap.html",
        stats="results/Taxonomic_rank_statistics.tsv",
        vir_stats="results/Virus_rank_statistics.tsv",
        phage_stats="results/Phage_rank_statistics.tsv",
        bact_stats="results/Bacteria_rank_statistics.tsv"
    conda:
        "../envs/heatmaps.yaml"
    benchmark:
        "logs/benchmark/draw_heatmaps.txt"
    threads: 1
    log:
        "logs/draw_heatmaps.log"
    shell:
        """
python bin/scripts/draw_heatmaps.py -c {input.classified} -n {input.numbers} -sq {output.super_quantities} -st {output.stats} -vs {output.vir_stats} -ps {output.phage_stats} -bs {output.bact_stats} -s {output.super} -v {output.virus} -p {output.phage} -b {output.bact} > {log} 2>&1
        """