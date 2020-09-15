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