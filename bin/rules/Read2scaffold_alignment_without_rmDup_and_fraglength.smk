
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#TODO This rule is only still included for debugging purposes, this, and the corresponding part in the main snakefile, can be removed in version 1.0

rule Read2scaffold_alignment_without_rmDup_and_fraglength:
    input:
        fasta   =   f"{datadir + scf_filt}" + "{sample}_scaffolds_ge%snt.fasta" % config["Illumina_meta"]["minlen"],
        pR1     =   f"{datadir + scf_filt}" + "{sample}_pR1.fq",
        pR2     =   f"{datadir + scf_filt}" + "{sample}_pR2.fq",
    output:
        bam     =   f"{datadir + scf_filt}" + "{sample}_sorted.bam",
        bam_bai =   f"{datadir + scf_filt}" + "{sample}_sorted.bam.bai",
        txt     =   f"{datadir + scf_filt}" + "{sample}_insert_size_metrics.txt",
        pdf     =   f"{datadir + scf_filt}" + "{sample}_insert_size_histogram.pdf"
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}" + "Read2scaffold_alignment_without_rmDup_and_fraglength_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Read2scaffold_alignment_without_rmDup_and_fraglength_{sample}.txt"
    threads: config["threads"]["Fragment_length_analysis"]
    shell:
        """
bwa index {input.fasta} > {log} 2>&1
bwa mem -t {threads} {input.fasta} \
{input.pR1} \
{input.pR2} 2>> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.bam} >> {log} 2>&1
samtools index -@ {threads} {output.bam} >> {log} 2>&1
picard -Dpicard.useLegacyParser=false CollectInsertSizeMetrics \
-I {output.bam} \
-O {output.txt} \
-H {output.pdf} >> {log} 2>&1
        """