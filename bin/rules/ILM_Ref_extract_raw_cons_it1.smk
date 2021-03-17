rule Illumina_extract_raw_consensus_it1:
    input:
        bam                 = rules.Illumina_align_to_reference_it1.output.sorted_bam,
        reference           = rules.Illumina_index_reference.output.reference_copy,
    output:
        unfiltered_vcf      = f"{datadir + it1 + cons}" + "{sample}_unfiltered.vcf",
        majorSNP_vcf        = f"{datadir + it1 + cons}" + "{sample}.vcf",
        majorSNP_vcf_gz     = f"{datadir + it1 + cons}" + "{sample}.vcf.gz",
        majorSNP_vcf_gz_tbi = f"{datadir + it1 + cons}" + "{sample}.vcf.gz.tbi", # Implicitly generated
        majorSNP_vcf_table  = f"{datadir + it1 + cons}" + "{sample}.tsv",
        raw_consensus_fasta = f"{datadir + it1 + cons}" + "{sample}_raw_consensus.fa",
        reference_copy_it2  = f"{datadir + it2 + refdir + reference_basename}" + "_{sample}.fasta",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_extract_raw_consensus_it1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_extract_raw_consensus_it1_{sample}.txt"
    threads: config["threads"]["Illumina_extract_raw_consensus"]
    resources: 
        memory = (config["threads"]["Illumina_extract_raw_consensus"] * 12) * 1024
    params:
        max_cov             = config["Illumina_ref"]["SNP"]["Max_coverage"],
    shell: # First codeblock, SNP and indel calling, filter majority SNPs from LoFreq output, zip/normalize/index it, call consensus based on these SNPs.
        """
lofreq call-parallel -d {params.max_cov} --no-default-filter --call-indels --use-orphan --pp-threads {threads} -f {input.reference} -o {output.unfiltered_vcf} {input.bam} >> {log} 2>&1

lofreq filter -a 0.50 -v 0 -i {output.unfiltered_vcf} -o {output.majorSNP_vcf} >> {log} 2>&1
bgzip -c {output.majorSNP_vcf} 2>> {log} |\
bcftools norm -m -both -O z -f {input.reference} -o {output.majorSNP_vcf_gz} >> {log} 2>&1
tabix {output.majorSNP_vcf_gz} >> {log} 2>&1
cat {input.reference} 2>> {log} |\
bcftools consensus {output.majorSNP_vcf_gz} 2>> {log} |\
seqtk seq - > {output.raw_consensus_fasta} 2>> {log}

cp {output.raw_consensus_fasta} {output.reference_copy_it2}

bcftools query -f '{wildcards.sample}\t%CHROM\t%POS\t%TYPE\t%REF\t%ALT{{0}}\t%QUAL\n' {output.majorSNP_vcf_gz} > {output.majorSNP_vcf_table} 2>> {log}
        """