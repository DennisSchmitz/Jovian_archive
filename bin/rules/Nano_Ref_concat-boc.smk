


rule concat_boc:
    input:
        boc_int = expand("{path}{sample}_BoC_int.tsv",
                            path = datadir + cons + boc,
                            sample = SAMPLES
                            ),
        boc_pct = expand("{path}{sample}_BoC_pct.tsv",
                            path = datadir + cons + boc,
                            sample = SAMPLES
                            ),
    output:
        conc_boc_int = res + "BoC_int.tsv",
        conc_boc_pct = res + "BoC_pct.tsv",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + bench + "Concat_boc.txt"
    log:
        logdir + "Concat_boc.log"
    threads: 1
    shell:
        """
echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.conc_boc_int}
cat {input.boc_int} >> {output.conc_boc_int}

echo -e "Sample_name\tTotal_ref_size\tBoC_at_coverage_threshold_1\tBoC_at_coverage_threshold_5\tBoC_at_coverage_threshold_10\tBoC_at_coverage_threshold_30\tBoC_at_coverage_threshold_100" > {output.conc_boc_pct}
cat {input.boc_pct} >> {output.conc_boc_pct}
        """ 