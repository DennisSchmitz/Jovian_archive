
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Generate_contigs_metrics:
    input:
        bam             =   rules.Read2scaffold_alignment_with_rmDup_and_fraglength.output.bam,
        fasta           =   rules.De_novo_assembly.output.filt_scaffolds,
        ORF_NT_fasta    =   rules.ORF_analysis.output.ORF_NT_fasta,
    output:
        summary         =   "data/scaffolds_filtered/{sample}_MinLenFiltSummary.stats",
        perScaffold     =   "data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats",
        perORFcoverage  =   "data/scaffolds_filtered/{sample}_perORFcoverage.stats",
    conda:
        conda_envs + "Sequence_analysis.yaml"
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