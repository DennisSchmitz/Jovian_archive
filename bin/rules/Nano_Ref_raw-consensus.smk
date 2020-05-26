

rule Create_raw_consensus:
    input:
        ref =   reference,
        vcf =   rules.Align_to_reference_pt2.output.vcf
    output: 
        consensus   =   f"{datadir + cons + raw}" + "{sample}.fasta"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}" + "Create_raw_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Create_raw_consensus_{sample}.txt"
    threads: 26
    shell:
        """
cat {input.ref} | bcftools consensus {input.vcf} 1> {output.consensus} 2> {log}
        """