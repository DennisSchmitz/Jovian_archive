
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Generate_contigs_metrics:
    input:
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        ORF_NT_fasta="data/scaffolds_filtered/{sample}_ORF_NT.fa",
    output:
        summary="data/scaffolds_filtered/{sample}_MinLenFiltSummary.stats",
        perScaffold="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats",
        perORFcoverage="data/scaffolds_filtered/{sample}_perORFcoverage.stats",
    conda:
        "../envs/scaffold_analyses.yaml"
    log:
        "logs/Generate_contigs_metrics_{sample}.log"
    benchmark:
        "logs/benchmark/Generate_contigs_metrics_{sample}.txt"
    params:
        ""
    threads: 1
    shell: #! bbtools' pileup.sh counts every read, even those marked as duplicate upstream. Hence, upstream all duplicates are HARD removed.
        """
pileup.sh in={input.bam} \
ref={input.fasta} \
fastaorf={input.ORF_NT_fasta} \
outorf={output.perORFcoverage} \
out={output.perScaffold} \
secondary=f \
samstreamer=t \
2> {output.summary} 1> {log}
        """