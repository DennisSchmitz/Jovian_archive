rule Illumina_concat_BoC_metrics:
    input:
        BoC_int_tsv =   expand( "{p}{sample}_BoC_int.tsv", 
                                p       =   f"{datadir + boc}",
                                sample  =   SAMPLES
                                ),
        BoC_pct_tsv =   expand( "{p}{sample}_BoC_pct.tsv",
                                p       =   f"{datadir + boc}",
                                sample  =   SAMPLES
                                )
    output:
        combined_BoC_int_tsv    =   f"{datadir + res}BoC_integer.tsv",
        combined_BoC_pct_tsv    =   f"{datadir + res}BoC_percentage.tsv",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}Illumina_concat_BoC_metrics.log"
    benchmark:
        f"{logdir + bench}Illumina_concat_BoC_metrics.txt"
    threads: 1
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_int_tsv}
cat {input.BoC_int_tsv} >> {output.combined_BoC_int_tsv}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.combined_BoC_pct_tsv}
cat {input.BoC_pct_tsv} >> {output.combined_BoC_pct_tsv}
        """