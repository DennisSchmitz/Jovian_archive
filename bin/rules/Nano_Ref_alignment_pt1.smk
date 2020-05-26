


rule Align_to_reference_pt1:
    input:
        ref     =   reference,
        fastq   =   rules.Hugo_removal_pt2.output.cleanedfastq
    output:
        bam         =   f"{datadir + aln + bf}" + "{sample}.bam",
        indexed_bam =   f"{datadir + aln + bf}" + "{sample}.bam.bai"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}" + "Align_to_reference_pt1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Align_to_reference_pt1_{sample}.txt"
    threads: 26
    shell:
        """
minimap2 -ax map-ont {input.ref} {input.fastq} 2>> {log} |\
samtools view -@ {threads} -F 256 -F 512 -F 4 -F 2048 -uS 2>> {log} |\
samtools sort -o {output.bam} >> {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
        """ 