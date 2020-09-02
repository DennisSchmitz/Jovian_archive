#TODO make a python script or bash function/include to do this more efficiently, currently it's hacky, but it works
rule Illumina_determine_BoC_at_diff_cov_thresholds:
    input:
        bedgraph    =   rules.Illumina_extract_clean_consensus.output.bedgraph,
        reference   =   rules.Illumina_extract_raw_consensus_it1.output.reference_copy_it2
    output:
        percentage_BoC_tsv  =   f"{datadir + boc}" + "{sample}_BoC_pct.tsv",
        integer_BoC_tsv     =   f"{datadir + boc}" + "{sample}_BoC_int.tsv",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_determine_BoC_at_diff_cov_thresholds_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_determine_BoC_at_diff_cov_thresholds_{sample}.txt"
    threads: 1
    shell:
        """
bash bin/scripts/BoC_analysis.sh {wildcards.sample} {input.bedgraph} {input.reference} \
{output.percentage_BoC_tsv} {output.integer_BoC_tsv} >> {log} 2>&1
        """