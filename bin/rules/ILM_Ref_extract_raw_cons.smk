rule Illumina_extract_raw_consensus:
    input:
        bam                 = rules.Illumina_align_to_reference.output.sorted_bam,
        reference           = rules.Illumina_index_reference.output.reference_copy,
    output:
        unfiltered_vcf      = f"{datadir + cons + raw}" + "{sample}_unfiltered.vcf",
        majorSNP_vcf        = f"{datadir + cons + raw}" + "{sample}.vcf",
        majorSNP_vcf_gz     = f"{datadir + cons + raw}" + "{sample}.vcf.gz",
        raw_consensus_fasta = f"{datadir + cons + raw}" + "{sample}_raw_consensus.fa",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_extract_raw_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_extract_raw_consensus_{sample}.txt"
    threads: config["threads"]["Illumina_extract_raw_consensus"]
    params:
    shell:
        """
lofreq call-parallel -d 20000 --no-default-filter --pp-threads {threads} -f {input.reference} -o {output.unfiltered_vcf} {input.bam} >> {log} 2>&1
lofreq filter -a 0.50 -i {output.unfiltered_vcf} -o {output.majorSNP_vcf} >> {log} 2>&1
bgzip -c {output.majorSNP_vcf} 2>> {log} |\
bcftools norm -m -both -O z -f {input.reference} -o {output.majorSNP_vcf_gz} >> {log} 2>&1
tabix {output.majorSNP_vcf_gz} >> {log} 2>&1
cat {input.reference} 2>> {log} |\
bcftools consensus {output.majorSNP_vcf_gz} 2>> {log} |\
seqtk seq - > {output.raw_consensus_fasta} 2>> {log}
        """