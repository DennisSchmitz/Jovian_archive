


rule calculate_BoC:
    input:
        bedgraph    =   rules.genomecoverage.output.bedgraph,
        reference   =   rules.Index_ref.output.refcopy
    output:
        pct_boc_tsv =   f"{datadir + cons + boc}" + "{sample}_BoC_pct.tsv",
        int_boc_tsv =   f"{datadir + cons + boc}" + "{sample}_BoC_int.tsv"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}" + "Calculate_boc_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Calculate_boc_{sample}.txt"
    threads: 1
    resources:
        memory = 8
    shell:
        """
bash bin/scripts/BoC_analysis.sh {wildcards.sample} {input.bedgraph} {input.reference} \
{output.pct_boc_tsv} {output.int_boc_tsv} >> {log} 2>&1
        """