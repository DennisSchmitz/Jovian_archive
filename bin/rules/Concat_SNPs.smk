
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
    shell: #TODO the sorting line below is dirty, but it works. Should be integrated into the python script on next code review.
        """
python bin/scripts/concat_filtered_vcf.py {params.vcf_folder_glob} results/all_filtered_SNPs.temp > {log} 2>&1
cat results/all_filtered_SNPs.temp | (read header; echo "$header"; sort -t$'\t' -k1,1 -k2,2 -k3,3n) 1> {output} 2>> {log}
rm results/all_filtered_SNPs.temp
        """