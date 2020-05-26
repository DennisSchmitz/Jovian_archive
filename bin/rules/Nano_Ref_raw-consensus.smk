

rule Create_raw_consensus:
    input:
        ref = reference,
        vcf = rules.Align_to_reference_pt2.output.vcf,
    output: 
        consensus = datadir + cons + raw + "{sample}.fasta"
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    benchmark:
        logdir + bench + "Create_raw_consensus_{sample}.txt"
    log:
        logdir + "Create_raw_consensus_{sample}.log"
    threads: 26
    shell:
        """
cat {input.ref} | bcftools consensus {input.vcf} 1> {output.consensus} 2> {log}
        """