
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HuGo_removal_pt2_extract_paired_unmapped_reads:
    input:
        bam        =   rules.HuGo_removal_pt1_alignment.output.sorted_bam,
        bam_index  =   rules.HuGo_removal_pt1_alignment.output.sorted_bam_index
    output:
        fastq_R1   =   f"{datadir + cln}" + "{sample}_pR1.fq",
        fastq_R2   =   f"{datadir + cln}" + "{sample}_pR2.fq"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    log:
        f"{logdir}" + "HuGo_removal_pt2_extract_paired_unmapped_reads_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "HuGo_removal_pt2_extract_paired_unmapped_reads_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    resources: 
        memory = config["threads"]["HuGo_removal"] * 12
    shell:
        """
samtools view -@ {threads} -b -f 1 -f 4 -f 8 {input.bam} 2> {log} |\
samtools sort -@ {threads} -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output.fastq_R1} -fq2 {output.fastq_R2} >> {log} 2>&1
        """