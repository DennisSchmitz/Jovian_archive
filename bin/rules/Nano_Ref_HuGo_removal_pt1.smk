



rule Hugo_removal_pt1:
    input:
        bg              =   config["databases"]["background_ref"],
        unmapped_fastq  =   rules.Cleanup.output.qc_fastq
    output: 
        sorted_bam          =   datadir + cln + hugo_no_rm + raw + "{sample}.bam",
        sorted_bam_index    =   datadir + cln + hugo_no_rm + raw + "{sample}.bam.bai"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    benchmark:
        logdir + bench + "Hugo_removal_pt1_{sample}.txt"
    log:
        logdir + "Hugo_removal_pt1_{sample}.log"
    threads: config["threads"]["HuGo_removal"]
    shell: 
        """
minimap2 -ax map-ont {input.bg} {input.unmapped_fastq} 2>> {log} |\
samtools view -uS 2>> {log} |\
samtools sort -o {output.sorted_bam} 2>> {log}
samtools index {output.sorted_bam} >> {log} 2>&1
        """