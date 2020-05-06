
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HuGo_removal_pt3_extract_unpaired_unmapped_reads:
    input:
        bam="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam",
        bam_index="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai"
    output:
        "data/cleaned_fastq/{sample}_unpaired.fq"
    conda:
        "../envs/HuGo_removal.yaml"
    log:
        "logs/HuGo_removal_pt3_extract_unpaired_unmapped_reads_{sample}.log"
    benchmark:
        "logs/benchmark/HuGo_removal_pt3_extract_unpaired_unmapped_reads_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
samtools view -@ {threads} -b -F 1 -f 4 {input.bam} 2> {log} |\
samtools sort -@ {threads} -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output} >> {log} 2>&1
        """