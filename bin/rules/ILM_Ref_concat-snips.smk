
rule SNP_table:
    input:
        expand( "{path}{sample}_calls.vcf.gz",
                path = f"{datadir + cons + raw}",
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
echo -e "Sample\tReference AccessionID\tPosition\tType\tReference\tAlternative\tQuality" > {output}

for x in {input}; do
    bcftools query -f '[%SAMPLE]\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT{{0}}\t%QUAL\n' $x >> {output}
done
        """