
rule SNP_table:
    input:
        expand( rules.Illumina_extract_raw_consensus_it1.output.majorSNP_vcf_table,
                sample  =   SAMPLES
                )
    output: f"{res}SNPs.tsv",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}Concatenated_SNPs.log"
    benchmark:
        f"{logdir + bench}Concatenated_SNPs.txt"
    threads: 1
    shell:
        """
echo -e "# The events below are incorporated into the consensus genome\nSample\tReference AccessionID\tPosition\tType\tReference\tAlternative\tQuality" > {output}

cat {input} >> {output}
        """