
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HuGo_removal_pt2_extract_paired_unmapped_reads:
    input:
         bam="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam",
         bam_index="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai",
    output:
         fastq_R1="data/cleaned_fastq/{sample}_pR1.fq.gz",
         fastq_R2="data/cleaned_fastq/{sample}_pR2.fq.gz",
    conda:
        "../envs/HuGo_removal.yaml"
    log:
        "logs/HuGo_removal_pt2_extract_paired_unmapped_reads_{sample}.log"
    benchmark:
        "logs/benchmark/HuGo_removal_pt2_extract_paired_unmapped_reads_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
reformat.sh in={input.bam} out1={output.fastq_R1} out2={output.fastq_R2} requiredbits=13 ow=t -da >> {log} 2>&1
        """
