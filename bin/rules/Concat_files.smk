
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Concat_files:
    input:
        expand( rules.Merge_all_metrics_into_single_tsv.output,
                sample = SAMPLES,
                extension = [   'taxClassified.tsv',
                                'taxUnclassified.tsv',
                                'virusHost.tsv'
                                ]
                ),
    output:
        taxClassified   =   f"{res}all_taxClassified.tsv",
        taxUnclassified =   f"{res}all_taxUnclassified.tsv",
        virusHost       =   f"{res}all_virusHost.tsv"
    log:
        f"{logdir}Concat_files.log"
    benchmark:
        f"{logdir + bench}Concat_files.txt"
    threads: 1
    resources:
        memory = 4
    params:
        search_folder       =   f"{datadir + tbl}",
        classified_glob     =   "*_taxClassified.tsv",
        unclassified_glob   =   "*_taxUnclassified.tsv",
        virusHost_glob      =   "*_virusHost.tsv"
    shell:
        """
find {params.search_folder} -type f -name "{params.classified_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + |\
(read header; echo "$header"; sort -t$'\t' -k 1,1 -k 15,15nr) 2> {log} 1> {output.taxClassified}

find {params.search_folder} -type f -name "{params.unclassified_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + |\
(read header; echo "$header"; sort -t$'\t' -k 1,1 -k 4,4nr) 2>> {log} 1> {output.taxUnclassified}

find {params.search_folder} -type f -name "{params.virusHost_glob}" -exec awk 'NR==1 || FNR!=1' {{}} + |\
(read header; echo "$header"; sort -t$'\t' -k 1,1 -k 2,2) 2>> {log} 1> {output.virusHost}
        """