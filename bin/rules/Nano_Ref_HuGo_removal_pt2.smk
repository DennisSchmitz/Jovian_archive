

rule Hugo_removal_pt2:
    input:
        bam         =   rules.Hugo_removal_pt1.output.sorted_bam,
        bam_index   =   rules.Hugo_removal_pt1.output.sorted_bam_index
    output:
        cleanedfastq    =   f"{datadir + cln}" + "{sample}.fq"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    log:
        f"{logdir}" + "Hugo_removal_pt2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Hugo_removal_pt2_{sample}.txt"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
samtools view -b -F 1 -f 4 {input.bam} 2>> {log} |\
samtools sort -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output.cleanedfastq} >> {log} 2>&1
        """