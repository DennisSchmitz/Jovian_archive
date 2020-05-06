
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Concat_filtered_SNPs:
    input:
        expand("data/scaffolds_filtered/{sample}_filtered.vcf", sample = SAMPLES)
    output:
        "results/all_filtered_SNPs.tsv"
    conda:
        "../envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/Concat_filtered_SNPs.txt"
    threads: 1
    params:
        vcf_folder_glob="data/scaffolds_filtered/\*_filtered.vcf"
    log:
        "logs/Concat_filtered_SNPs.log"
    shell:
        """
python bin/scripts/concat_filtered_vcf.py {params.vcf_folder_glob} {output} > {log} 2>&1
        """