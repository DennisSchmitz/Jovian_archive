rule genomecoverage:
    input:
        bam         =   rules.Align_to_reference_pt1.output.bam
    output:
        bedgraph    =   f"{datadir + cons}" + "{sample}.bedgraph"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}" + "Determine_genome_coverage_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Determine_genome_coverage_{sample}.txt"
    threads: 1
    resources:
        memory = 8 * 1024
    shell:
        """
bedtools genomecov -bga -ibam {input.bam} > {output.bedgraph} 2>> {log}
        """