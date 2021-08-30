
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HuGo_removal_pt3_extract_unpaired_unmapped_reads:
    input:
        bam         =   rules.HuGo_removal_pt1_alignment.output.sorted_bam,
        bam_index   =   rules.HuGo_removal_pt1_alignment.output.sorted_bam_index
    output:
        tbam = temp(f"{datadir + cln}" + "{sample}_tempbam_unpaired.bam"),
        fq = f"{datadir + cln}" + "{sample}_unpaired.fq"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    log:
        f"{logdir}" + "HuGo_removal_pt3_extract_unpaired_unmapped_reads_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "HuGo_removal_pt3_extract_unpaired_unmapped_reads_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    resources: 
        memory = (config["threads"]["HuGo_removal"] * 12) * 1024
    shell:
        """
samtools view -@ {threads} -b -F 1 -f 4 {input.bam} 2> {log} |\
samtools sort -@ {threads} -n -o {output.tbam} 2>> {log}
bedtools bamtofastq -i {output.tbam} -fq {output.fq} >> {log} 2>&1
        """