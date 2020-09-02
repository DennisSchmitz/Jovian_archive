rule Illumina_extract_raw_consensus_it2:
    input:
        bam                 = rules.Illumina_align_to_reference_it2.output.sorted_bam,
        reference           = rules.Illumina_extract_raw_consensus_it1.output.reference_copy_it2,
    output:
        unfiltered_vcf      = f"{datadir + it2 + cons}" + "{sample}_unfiltered.vcf",
        majorSNP_vcf        = f"{datadir + it2 + cons}" + "{sample}.vcf",
        majorSNP_vcf_gz     = f"{datadir + it2 + cons}" + "{sample}.vcf.gz",
        majorSNP_vcf_gz_tbi = f"{datadir + it2 + cons}" + "{sample}.vcf.gz.tbi", # Implicitly generated
        majorSNP_vcf_table  = f"{datadir + it2 + cons}" + "{sample}.tsv",
        raw_consensus_fasta = f"{datadir + it2 + cons}" + "{sample}_raw_consensus.fa",
        minorSNP_vcf        = f"{datadir + it2 + cons}" + "{sample}_minorSNPs.vcf",
        minorSNP_vcf_gz     = f"{datadir + it2 + cons}" + "{sample}_minorSNPs.vcf.gz",
        minorSNP_vcf_gz_tbi = f"{datadir + it2 + cons}" + "{sample}_minorSNPs.vcf.gz.tbi", # Implicitly generated
        minorSNP_vcf_table  = f"{datadir + it2 + cons}" + "{sample}_minorSNPs.tsv",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_extract_raw_consensus_it2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_extract_raw_consensus_it2_{sample}.txt"
    threads: config["threads"]["Illumina_extract_raw_consensus"]
    params:
        min_AF              = "0.05"
    shell: # First codeblock, SNP and indel calling, filter majority SNPs from LoFreq output, zip/normalize/index it, call consensus based on these SNPs. Third codeblock, filter minority SNPs, flag them in IGVjs for manual inspection and/or tabular output.
        """
lofreq call-parallel -d 20000 --no-default-filter --call-indels --use-orphan --pp-threads {threads} -f {input.reference} -o {output.unfiltered_vcf} {input.bam} >> {log} 2>&1

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

bcftools query -f '{wildcards.sample}\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT{{0}}\t%QUAL\n' {output.majorSNP_vcf_gz} > {output.majorSNP_vcf_table} 2>> {log}

echo -e "\nSample\tReference AccessionID\tPosition\tType\tReference\tAlternative\tQuality\tRaw Depth\tAllele Frequency\tPhred-scaled strand bias at this position\tCounts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases\tHomopolymer length to the right of report indel position" > {output.minorSNP_vcf_table} 2>> {log}
bcftools query -f '{wildcards.sample}\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT{{0}}\t%QUAL\t%DP\t%AF\t%SB\t%DP4\t%HRUN\n' {output.minorSNP_vcf_gz} >> {output.minorSNP_vcf_table} 2>> {log}
        """