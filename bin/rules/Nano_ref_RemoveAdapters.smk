rule Remove_Adapters_pt1:
    input:
        ref =   rules.Index_ref.output.refcopy 
        fq  =   rules.Cleanup.output.qc_fastq
    output: 
        bam = f"{datadir + cln + trims}" + "{sample}.bam"
        index = f"{datadir + cln + trims}" + "{sample}.bam.bai"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    shell: 
        """
minimap2 -ax map-ont {input.ref} {input.fq} 2>> {log} |\
samtools view -@ {threads} -F 256 -F 512 -F 4 -F 2048 -uS 2>> {log} |\
samtools sort -o {output.bam} >> {log} 2>&1
samtools index {output.bam} >> {log} 2>&1
        """

rule Remove_Adapters_pt2:
    input: rules.Remove_Adapters_pt1.output.bam
    output: f"{datadir + cln + trims}" + "{sample}.fastq"
    threads: 2
    shell: 
        """
python bin/scripts/SoftClipper.py --input {input} --output {output} --threads {threads}
        """