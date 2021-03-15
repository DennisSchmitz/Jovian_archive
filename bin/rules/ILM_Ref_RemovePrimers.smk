rule RemovePrimers_pt1:
    input:
        r1  = rules.Clean_the_data.output.r1,
        r2  = rules.Clean_the_data.output.r2,
        ur1 = rules.Clean_the_data.output.r1_unpaired,
        ur2 = rules.Clean_the_data.output.r2_unpaired
    output:
        fastq =   f"{datadir + cln + prdir}" + "{sample}_untrimmed.fastq"
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_RemovePrimers_pt1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_RemovePrimers_pt1_{sample}.txt"
    threads: 1
    resources:
        memory = 4 * 1024
    shell:
        """
cat {input} > {output}
        """

rule RemovePrimers_pt2:
    input: 
        fastq   = rules.RemovePrimers_pt1.output.fastq,
        ref     = rules.Illumina_index_reference.output.reference_copy,
        primers = primerfile
    output: f"{datadir + cln + prdir}" + "{sample}.fastq"
    conda:
        f"{conda_envs}Illumina_ref_CutPrimers.yaml"
    log:
        f"{logdir}" + "Illumina_RemovePrimers_pt2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_RemovePrimers_pt2_{sample}.txt"
    threads: config["threads"]["Illumina_RemovePrimers"]
    resources: 
        memory = (config["threads"]["Illumina_RemovePrimers"] * 12) * 1024
    params:
        primer_status = prstatus
    shell:
        """
if [ {params.primer_status} == "SET" ]; then
    python bin/scripts/RemoveIlluminaPrimers.py --input {input.fastq} --reference {input.ref} --primers {input.primers} --threads {threads} --output {output}
else
    cp {input.fastq} {output}
fi
        """