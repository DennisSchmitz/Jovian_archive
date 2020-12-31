rule SNP_table:
    input:
        expand( rules.Align_to_reference_pt2.output.vcf_table,
                sample  =   SAMPLES
                )
    output: f"{res}SNPs.tsv"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}Concatenated_SNPs.log"
    benchmark:
        f"{logdir + bench}Concatenated_SNPs.txt"
    threads: 1
    shell:
        """
echo -e "Sample\tReference AccessionID\tPosition\tType\tReference\tAlternative\tQuality" > {output}

cat {input} >> {output}
        """

rule concat_ins:
    input:
        expand( rules.consensus_cov_1.output.ins,
                sample = SAMPLES
            )
    output: f"{res}INS.tsv"
    log:
        f"{logdir}concat_ins.log"
    benchmark:
        f"{logdir + bench}concat_ins.txt"
    threads: 1
    shell:
        """
echo -e "Sample\tHas Insertions?\tPositions\tPercentage of reads" > {output}

cat {input} >> {output}
        """