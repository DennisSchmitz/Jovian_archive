


rule calculate_BoC:
    input:
        bedgraph = rules.extract_cleaned_consensus.output.bedgraph,
        reference = rules.Index_ref.output.refcopy,
    output:
        pct_boc_tsv = datadir + cons + boc + "{sample}_BoC_pct.tsv",
        int_boc_tsv = datadir + cons + boc + "{sample}_BoC_int.tsv",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + bench + "Calculate_boc_{sample}.txt"
    log:
        logdir + "Calculate_boc_{sample}.log"
    threads: 1
    shell:
        """
bash bin/scripts/RA_BoC_analysis.sh {wildcards.sample} {input.bedgraph} {input.reference} \
{output.pct_boc_tsv} {output.int_boc_tsv} >> {log} 2>&1
        """