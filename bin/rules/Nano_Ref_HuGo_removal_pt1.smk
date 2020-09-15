



rule Hugo_removal_pt1:
    input:
        bg              =   config["databases"]["background_ref"],
        unmapped_fastq  =   rules.Cleanup.output.qc_fastq
    output: 
        sorted_bam          =   f"{datadir + cln + hugo_no_rm + raw}" + "{sample}.bam",
        sorted_bam_index    =   f"{datadir + cln + hugo_no_rm + raw}" + "{sample}.bam.bai"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    log:
        f"{logdir}" + "Hugo_removal_pt1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Hugo_removal_pt1_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    shell: 
        """
minimap2 -ax map-ont {input.bg} {input.unmapped_fastq} 2>> {log} |\
samtools view -uS 2>> {log} |\
samtools sort -o {output.sorted_bam} 2>> {log}
samtools index {output.sorted_bam} >> {log} 2>&1
        """