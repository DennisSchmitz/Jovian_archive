
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HuGo_removal_pt1_alignment:
    input:
        background_ref=config["databases"]["background_ref"],
        r1="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR1.fastq",
        r2="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_pR2.fastq",
        r1_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR1.fastq",
        r2_unpaired="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_uR2.fastq",
    output:
        sorted_bam="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam",
        sorted_bam_index="data/cleaned_fastq/fastq_without_HuGo_removal/{sample}_sorted.bam.bai",
    conda:
        "../envs/HuGo_removal.yaml"
    benchmark:
        "logs/benchmark/HuGo_removal_pt1_alignment_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    params:
        alignment_type=config["HuGo_removal"]["alignment_type"]
    log:
        "logs/HuGo_removal_pt1_alignment_{sample}.log"
    shell:
        """
bowtie2 --time --threads {threads} {params.alignment_type} \
-x {input.background_ref} \
-1 {input.r1} \
-2 {input.r2} \
-U {input.r1_unpaired} \
-U {input.r2_unpaired} 2> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.sorted_bam} >> {log} 2>&1
samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        """