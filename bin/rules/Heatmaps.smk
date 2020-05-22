
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule draw_heatmaps:
    input:
        classified  =   rules.Concat_files.output.taxClassified,
        numbers     =   f"{res + mqc_data}multiqc_trimmomatic.txt"
    output:
        super_quantities    =   f"{res}Superkingdoms_quantities_per_sample.csv",
        super               =   f"{res + hmap}Superkingdoms_heatmap.html",
        virus               =   f"{res + hmap}Virus_heatmap.html",
        phage               =   f"{res + hmap}Phage_heatmap.html",
        bact                =   f"{res + hmap}Bacteria_heatmap.html",
        stats               =   f"{res}Taxonomic_rank_statistics.tsv",
        vir_stats           =   f"{res}Virus_rank_statistics.tsv",
        phage_stats         =   f"{res}Phage_rank_statistics.tsv",
        bact_stats          =   f"{res}Bacteria_rank_statistics.tsv"
    conda:
        f"{conda_envs}heatmaps.yaml"
    log:
        f"{logdir}draw_heatmaps.log"
    benchmark:
        f"{logdir + bench}draw_heatmaps.txt"
    threads: 1
    shell:
        """
python bin/scripts/draw_heatmaps.py -c {input.classified} -n {input.numbers} -sq {output.super_quantities} -st {output.stats} -vs {output.vir_stats} -ps {output.phage_stats} -bs {output.bact_stats} -s {output.super} -v {output.virus} -p {output.phage} -b {output.bact} > {log} 2>&1
        """