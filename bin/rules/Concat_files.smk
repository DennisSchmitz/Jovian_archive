
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Concat_files:
    input:
        expand("data/tables/{sample}_{extension}", sample = SAMPLES, extension = ['taxClassified.tsv','taxUnclassified.tsv','virusHost.tsv']),
    output:
        taxClassified="results/all_taxClassified.tsv",
        taxUnclassified="results/all_taxUnclassified.tsv",
        virusHost="results/all_virusHost.tsv",
    benchmark:
        "logs/benchmark/Concat_files.txt"
    threads: 1
    log:
        "logs/Concat_files.log"
    params:
        search_folder="data/tables/",
        classified_glob="*_taxClassified.tsv",
        unclassified_glob="*_taxUnclassified.tsv",
        virusHost_glob="*_virusHost.tsv",
    shell:
        """
find {params.search_folder} -type f -name "{params.classified_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + |\
(read header; echo "$header"; sort -t$'\t' -k 1,1 -k 15,15nr) 2> {log} 1> {output.taxClassified}

find {params.search_folder} -type f -name "{params.unclassified_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + |\
(read header; echo "$header"; sort -t$'\t' -k 1,1 -k 15,15nr) 2>> {log} 1> {output.taxUnclassified}

find {params.search_folder} -type f -name "{params.virusHost_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + |\
(read header; echo "$header"; sort -t$'\t' -k 1,1 -k 15,15nr) 2>> {log} 1> {output.virusHost}
        """