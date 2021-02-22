rule RemovePrimers_pt1:
    input:
        ref = rules.Illumina_index_reference.output.reference_copy,
        r1  = rules.Clean_the_data.output.r1,
        r2  = rules.Clean_the_data.output.r2,
        ur1 = rules.Clean_the_data.output.r1_unpaired,
        ur2 = rules.Clean_the_data.output.r2_unpaired
    output:
        bam =   f"{datadir + cln + prdir}" + "{sample}_sorted.bam",
        bai =   f"{datadir + cln + prdir}" + "{sample}_sorted.bam.bai"
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_RemovePrimers_pt1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_RemovePrimers_pt1_{sample}.txt"
    threads: config["threads"]["Illumina_align_to_reference"]
    params:
        aln_type    =   config["Illumina_ref"]["Alignment"]["Alignment_type"]
    shell:
        """
minimap2 -a -t {threads} {input} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools collate -@ {threads} -O - 2>> {log} |\
samtools fixmate -@ {threads} -m - - 2>> {log} |\
samtools sort -@ {threads} - -o {output.bam} 2>> {log}
samtools index -@ {threads} {output.bam} >> {log} 2>&1
        """
        #TODO: specify minimap input with -x flag. "-x sr" does not handle more than 2 input files
        
        # bowtie2 --time --threads {threads} {params.aln_type} \
        # -x {input.ref} \
        # -1 {input.r1} \
        # -2 {input.r2} \
        # -U {input.ur1},{input.ur2} 2> {log} |\

rule RemovePrimers_pt2:
    input: 
        bam = rules.RemovePrimers_pt1.output.bam,
        ref = rules.Illumina_index_reference.output.reference_copy,
        primers = primerfile
    output: f"{datadir + cln + prdir}" + "{sample}.fastq"
    conda:
        f"{conda_envs}Illumina_ref_CutPrimers.yaml"
    log:
        f"{logdir}" + "Illumina_RemovePrimers_pt2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_RemovePrimers_pt2_{sample}.txt"
    threads: config["threads"]["Illumina_RemovePrimers"]
    params:
        primer_status = prstatus
    shell:
        """
if [ {params.primer_status} == "SET" ]; then
    python3.8 bin/scripts/RemoveIlluminaPrimers.py --input {input.bam} --reference {input.ref} --primers {input.primers} --threads {threads} --output {output} 2>> {log}
else
    bedtools bamtofastq -i {input} -fq {output} 2>> {log}
fi
        """