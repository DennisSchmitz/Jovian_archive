
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Concat_filtered_SNPs:
    input:
        expand( rules.SNP_calling.output.filt_vcf, 
                sample = SAMPLES
                )
    output:
        f"{res}all_filtered_SNPs.tsv"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}Concat_filtered_SNPs.log"
    benchmark:
        f"{logdir + bench}Concat_filtered_SNPs.txt"
    threads: 1
    params:
        vcf_folder_glob =   f"{datadir + scf_filt}/\*_filtered.vcf"
    shell:
        """
python bin/scripts/concat_filtered_vcf.py {params.vcf_folder_glob} {output} > {log} 2>&1
        """