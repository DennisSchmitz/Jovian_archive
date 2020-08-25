



rule Align_to_reference_pt2:
    input:
        ref =   rules.Index_ref.output.refcopy,
        bam =   rules.Align_to_reference_pt1.output.bam
    output:
        vcf         =   f"{datadir + aln + vf}" + "{sample}.vcf.gz",
        vcf_index   =   f"{datadir + aln + vf}" + "{sample}.vcf.gz.csi",
        vcf_table   =   f"{datadir + aln + vf}" + "{sample}.vcf.tsv"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}" + "Align_to_reference_pt2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Align_to_reference_pt2_{sample}.txt"
    threads: 1
    shell:
        """
bcftools mpileup --ignore-RG -Ou -d 10000 -f {input.ref} {input.bam} 2>> {log} |\
bcftools call --ploidy 1 -mv -Oz 2>> {log} |\
bcftools norm -m -both -O z -f {input.ref} -o {output.vcf} >> {log} 2>&1
tabix {output.vcf} >> {log} 2>&1
bcftools index {output.vcf}
bcftools query -f '{wildcards.sample}\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT{{0}}\t%QUAL\n' {output.vcf} > {output.vcf_table}
        """ 