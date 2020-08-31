rule Illumina_determine_GC_content:
    input:
        fasta   =   rules.Illumina_extract_raw_consensus_it1.output.reference_copy_it2,
    output:
        fasta_fai   =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}.fasta.fai",
        fasta_sizes =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}.fasta.sizes",
        bed_windows =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}.windows",
        GC_bed      =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}_GC.bedgraph",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}" + "Illumina_determine_GC_content_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_determine_GC_content_{sample}.txt"
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