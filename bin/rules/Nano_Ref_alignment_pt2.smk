



rule Align_to_reference_pt2:
    input:
        ref = rules.Index_ref.output.refcopy,
        bam = rules.Align_to_reference_pt1.output.bam,
    output:
        vcf = datadir + aln + vf + "{sample}.vcf.gz",
        vcf_index = datadir + aln + vf + "{sample}.vcf.gz.csi",
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + bench + "Align_to_reference_pt2_{sample}.txt"
    log:
        logdir + "Align_to_reference_pt2_{sample}.log"
    threads: 26
    shell:
        """
bcftools mpileup --ignore-RG -Ou -d 10000 -f {input.ref} {input.bam} 2>> {log} |\
bcftools call --ploidy 1 -mv -Oz 2>> {log} |\
bcftools norm -m -both -O z -f {input.ref} -o {output.vcf} >> {log} 2>&1
tabix {output.vcf} >> {log} 2>&1
bcftools index {output.vcf}
        """ 