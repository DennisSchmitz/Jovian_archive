

rule Hugo_removal_pt2:
    input:
        bam         =   rules.Hugo_removal_pt1.output.sorted_bam,
        bam_index   =   rules.Hugo_removal_pt1.output.sorted_bam_index
    output:
        cleanedfastq = datadir + cln + "{sample}.fq"
    conda:
        f"{conda_envs}HuGo_removal.yaml"
    benchmark:
        logdir + bench + "Hugo_removal_pt2_{sample}.txt"
    log:
        logdir + "Hugo_removal_pt2_{sample}.log"
    threads: config["threads"]["HuGo_removal"]
    shell:
        """
samtools view -b -F 1 -f 4 {input.bam} 2>> {log} |\
samtools sort -n - 2>> {log} |\
bedtools bamtofastq -i - -fq {output.cleanedfastq} >> {log} 2>&1
        """