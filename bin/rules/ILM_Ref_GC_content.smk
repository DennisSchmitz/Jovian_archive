


rule Illumina_determine_GC_content:
    input:
        fasta   =   rules.Illumina_index_reference.output.reference_copy,
    output:
        fasta_fai   =   f"{datadir + refdir + reference_basename}.fasta.fai",
        fasta_sizes =   f"{datadir + refdir + reference_basename}.fasta.sizes",
        bed_windows =   f"{datadir + refdir + reference_basename}.windows",
        GC_bed      =   f"{datadir + refdir + reference_basename}_GC.bedgraph",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}Illumina_determine_GC_content.log"
    benchmark:
        f"{logdir + bench}Illumina_determine_GC_content.txt"
    threads: 1
    params:
        window_size =   config["Global"]["GC_window_size"]
    shell:
        """
samtools faidx -o {output.fasta_fai} {input.fasta} > {log} 2>&1
cut -f 1,2 {output.fasta_fai} 2> {log} 1> {output.fasta_sizes}
bedtools makewindows \
-g {output.fasta_sizes} \
-w {params.window_size} 2>> {log} 1> {output.bed_windows}
bedtools nuc \
-fi {input.fasta} \
-bed {output.bed_windows} 2>> {log} |\
cut -f 1-3,5 2>> {log} 1> {output.GC_bed}
        """