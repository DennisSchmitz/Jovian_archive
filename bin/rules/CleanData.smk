
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Clean_the_data:
    input:
        lambda wildcards: (SAMPLES[wildcards.sample][i]
                            for i in ("R1", "R2")
                            )
    output:
        r1          =   f"{datadir + cln + hugo_no_rm}" + "{sample}_pR1.fastq",
        r2          =   f"{datadir + cln + hugo_no_rm}" + "{sample}_pR2.fastq",
        r1_unpaired =   f"{datadir + cln + hugo_no_rm}" + "{sample}_uR1.fastq",
        r2_unpaired =   f"{datadir + cln + hugo_no_rm}" + "{sample}_uR2.fastq",
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "Clean_the_data_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Clean_the_data_{sample}.txt"
    threads: config["threads"]["Clean_the_data"]
    resources:
        memory = (config["threads"]["Clean_the_data"] * 4) * 1024
    params:
        adapter_removal_config  =   config["Illumina"]["Clean"]["adapter_removal_config"],
        quality_trimming_config =   config["Illumina"]["Clean"]["quality_trimming_config"],
        minimum_length_config   =   config["Illumina"]["Clean"]["minimum_length_config"],
    shell:
        """
trimmomatic PE -threads {threads} \
{input[0]:q} {input[1]:q} \
{output.r1} {output.r1_unpaired} \
{output.r2} {output.r2_unpaired} \
{params.adapter_removal_config} \
{params.quality_trimming_config} \
{params.minimum_length_config} > {log} 2>&1
touch -r {output.r1} {output.r1_unpaired}
touch -r {output.r2} {output.r2_unpaired}
        """