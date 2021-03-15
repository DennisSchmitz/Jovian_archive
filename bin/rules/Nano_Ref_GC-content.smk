

rule determine_GC_content:
    input:
        ref =   rules.Index_ref.output.refcopy,
    output:
        fai         =   f"{datadir + refdir + reference_basename}.fasta.fai",
        sizes       =   f"{datadir + refdir + reference_basename}.fasta.sizes",
        bed_windows =   f"{datadir + refdir + reference_basename}.windows",
        GC_bed      =   f"{datadir + refdir + reference_basename}_GC.bedgraph",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}Determine_GC_content.log"
    benchmark:
        f"{logdir + bench}Determine_GC_content.txt"
    threads: 1
    resources:
        memory = 8
    params:
        window_size =   config["Global"]["GC_window_size"]
    shell:
        """
samtools faidx -o {output.fai} {input.ref} > {log} 2>&1
cut -f 1,2 {output.fai} 2> {log} 1> {output.sizes}
bedtools makewindows \
-g {output.sizes} \
-w {params.window_size} 2>> {log} 1> {output.bed_windows}
bedtools nuc \
-fi {input.ref} \
-bed {output.bed_windows} 2>> {log} |\
cut -f 1-3,5 2>> {log} 1> {output.GC_bed}
        """ 