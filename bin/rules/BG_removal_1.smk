
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HuGo_removal_pt1_alignment:
    input:
        background_ref  =   config["databases"]["background_ref"],
        r1              =   rules.Clean_the_data.output.r1,
        r2              =   rules.Clean_the_data.output.r2,
        r1_unpaired     =   rules.Clean_the_data.output.r1_unpaired,
        r2_unpaired     =   rules.Clean_the_data.output.r2_unpaired
    output:
        sorted_bam          =   f"{datadir + cln + hugo_no_rm}" + "{sample}_sorted.bam",
        sorted_bam_index    =   f"{datadir + cln + hugo_no_rm}" + "{sample}_sorted.bam.bai"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    log:
        f"{logdir}" + "HuGo_removal_pt1_alignment_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "HuGo_removal_pt1_alignment_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    resources: 
        memory = (config["threads"]["HuGo_removal"] * 12) * 1024
    params:
        alignment_type  =   config["Global"]["HuGo_removal_method"]
    shell:
        """
bowtie2 --time --threads {threads} {params.alignment_type} \
-x {input.background_ref} \
-1 {input.r1} \
-2 {input.r2} \
-U {input.r1_unpaired} \
-U {input.r2_unpaired} 2> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools sort -@ {threads} - -o {output.sorted_bam} >> {log} 2>&1
samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        """