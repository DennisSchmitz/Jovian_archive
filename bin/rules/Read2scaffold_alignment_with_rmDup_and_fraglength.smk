
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Read2scaffold_alignment_with_rmDup_and_fraglength:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        pR1="data/cleaned_fastq/{sample}_pR1.fq",
        pR2="data/cleaned_fastq/{sample}_pR2.fq",
    output:
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        bam_bai="data/scaffolds_filtered/{sample}_sorted.bam.bai",
        dup_metrics= "data/scaffolds_filtered/{sample}_sorted.MarkDup_metrics", #TODO deze toevoegen aan MultiQC?
        txt="data/scaffolds_filtered/{sample}_insert_size_metrics.txt",
        pdf="data/scaffolds_filtered/{sample}_insert_size_histogram.pdf"
    conda:
        "../envs/scaffold_analyses.yaml"
    log:
        "logs/Read2scaffold_alignment_with_rmDup_and_fraglength_{sample}.log"
    benchmark:
        "logs/benchmark/Read2scaffold_alignment_with_rmDup_and_fraglength_{sample}.txt"
    threads: config["threads"]["Fragment_length_analysis"]
    params:
        remove_dups="-r", #! This will HARD remove the duplicates instead of only marking them, N.B. this is REQUIRED for the downstream bbtools' pileup.sh to work --> it ignores the DUP marker and counts the reads in its coverage metrics. Thus giving a false sense of confidence.
        markdup_mode="t",
        max_read_length="300", # This is the default value and also the max read length of Illumina in-house sequencing.
    shell:
        """
bwa index {input.fasta} > {log} 2>&1
bwa mem -t {threads} {input.fasta} \
{input.pR1} \
{input.pR2} 2>> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools collate -@ {threads} -O - 2>> {log} |\
samtools fixmate -@ {threads} -m - - 2>> {log} |\
samtools sort -@ {threads} - -o - 2>> {log} |\
samtools markdup -@ {threads} -l {params.max_read_length} -m {params.markdup_mode} {params.remove_dups} -f {output.dup_metrics} - {output.bam} >> {log} 2>&1
samtools index -@ {threads} {output.bam} >> {log} 2>&1

picard -Dpicard.useLegacyParser=false CollectInsertSizeMetrics \
-I {output.bam} \
-O {output.txt} \
-H {output.pdf} >> {log} 2>&1
        """

        