threadcounts = config["threads"]["Nanopore_reference_alignment"]
mappingthreads = threadcounts - 1 

rule Remove_Adapters_pt1:
    input:
        ref =   rules.Index_ref.output.refcopy, 
        fq  =   lambda wildcards: SAMPLES[wildcards.sample]
    output: 
        bam = f"{datadir + cln + trims}" + "{sample}.bam",
        index = f"{datadir + cln + trims}" + "{sample}.bam.bai"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    log:
        f"{logdir}" + "Remove_Adapters_pt1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Remove_Adapters_pt1_{sample}.txt"
    threads: config["threads"]["Nanopore_reference_alignment"]
    resources: 
        memory = config["threads"]["Nanopore_reference_alignment"] * 12
    params:
        mapthreads = mappingthreads
    shell: 
        """
minimap2 -ax map-ont -t {params.mapthreads} {input.ref} {input.fq} 2>> {log} |\
samtools view -@ {threads} -F 256 -F 512 -F 4 -F 2048 -uS 2>> {log} |\
samtools sort -o {output.bam} >> {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
        """

rule Remove_Adapters_pt2:
    input: rules.Remove_Adapters_pt1.output.bam
    output: f"{datadir + cln + trims}" + "{sample}.fastq"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    log:
        f"{logdir}" + "Remove_Adapters_pt2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Remove_Adapters_pt2_{sample}.txt"
    threads: config["threads"]["Nanopore_SoftClipper"]
    resources:
        memory = 12
    shell: 
        """
python bin/scripts/SoftClipper.py --input {input} --output {output} --threads {threads}
        """