


rule concat_boc:
    input:
        boc_int =   expand( "{p}{sample}_BoC_int.tsv",
                            p = f"{datadir + cons + boc}",
                            sample = SAMPLES
                            ),
        boc_pct =   expand( "{p}{sample}_BoC_pct.tsv",
                            p = f"{datadir + cons + boc}",
                            sample = SAMPLES
                            )
    output:
        conc_boc_int    =   f"{res}BoC_int.tsv",
        conc_boc_pct    =   f"{res}BoC_pct.tsv",
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}Concat_boc.log"
    benchmark:
        f"{logdir + bench}Concat_boc.txt"
    threads: 1
    resources:
        memory = 6
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.conc_boc_int}
cat {input.boc_int} >> {output.conc_boc_int}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.conc_boc_pct}
cat {input.boc_pct} >> {output.conc_boc_pct}
        """ 