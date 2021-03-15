
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
        summary         =   f"{datadir + scf_filt}" + "{sample}_MinLenFiltSummary.stats",
        perScaffold     =   f"{datadir + scf_filt}" + "{sample}_perMinLenFiltScaffold.stats",
        perORFcoverage  =   f"{datadir + scf_filt}" + "{sample}_perORFcoverage.stats",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}" + "Generate_contigs_metrics_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Generate_contigs_metrics_{sample}.txt"
    threads: 1
    resources:
        memory = 8
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