rule Illumina_extract_raw_consensus:
    input:
        bam                 = rules.Illumina_align_to_reference.output.sorted_bam,
        reference           = rules.Illumina_index_reference.output.reference_copy,
    output:
        indelqual_bam       = f"{datadir + cons + raw}" + "{sample}_indelqual.bam",
        unfiltered_vcf      = f"{datadir + cons + raw}" + "{sample}_unfiltered.vcf",
        majorSNP_vcf        = f"{datadir + cons + raw}" + "{sample}.vcf",
        majorSNP_vcf_gz     = f"{datadir + cons + raw}" + "{sample}.vcf.gz",
        majorSNP_vcf_table  = f"{datadir + cons + raw}" + "{sample}.vcf.gz.tsv",
        raw_consensus_fasta = f"{datadir + cons + raw}" + "{sample}_raw_consensus.fa",
        minorSNP_vcf        = f"{datadir + cons + raw}" + "{sample}_minorSNPs.vcf",
        minorSNP_vcf_gz     = f"{datadir + cons + raw}" + "{sample}_minorSNPs.vcf.gz"
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_extract_raw_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_extract_raw_consensus_{sample}.txt"
    threads: config["threads"]["Illumina_extract_raw_consensus"]
    params:
        min_AF              = "0.05"
    shell: # First codeblock, set indelqual (requirement for LoFreq indel calling of Ill-data), index updated bam, SNP and indels calling. Second codeblock, filter majority SNPs from LoFreq output, zip/normalize/index it, call consensus based on these SNPs. Third codeblock, filter minority SNPs, flag them in IGVjs for manual inspection and/or tabular output.
        """
lofreq indelqual --dindel -f {input.reference} -o {output.indelqual_bam} {input.bam} >> {log} 2>&1
samtools index -@ {threads} {output.indelqual_bam} >> {log} 2>&1
lofreq call-parallel -d 20000 --no-default-filter --call-indels --use-orphan --pp-threads {threads} -f {input.reference} -o {output.unfiltered_vcf} {output.indelqual_bam} >> {log} 2>&1

lofreq filter -a 0.50 -v 0 -i {output.unfiltered_vcf} -o {output.majorSNP_vcf} >> {log} 2>&1
bgzip -c {output.majorSNP_vcf} 2>> {log} |\
bcftools norm -m -both -O z -f {input.reference} -o {output.majorSNP_vcf_gz} >> {log} 2>&1
tabix {output.majorSNP_vcf_gz} >> {log} 2>&1
cat {input.reference} 2>> {log} |\
bcftools consensus {output.majorSNP_vcf_gz} 2>> {log} |\
seqtk seq - > {output.raw_consensus_fasta} 2>> {log}

lofreq filter -a {params.min_AF} -A 0.50  -v 0 -i {output.unfiltered_vcf} -o {output.minorSNP_vcf} >> {log} 2>&1
bgzip -c {output.minorSNP_vcf} 2>> {log} |\
bcftools norm -m -both -O z -f {input.reference} -o {output.minorSNP_vcf_gz} - >> {log} 2>&1
tabix {output.minorSNP_vcf_gz} >> {log} 2>&1


bcftools query -f '{wildcards.sample}\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT{{0}}\t%QUAL\n' {output.majorSNP_vcf_gz} > {output.majorSNP_vcf_table}
        """