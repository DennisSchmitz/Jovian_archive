
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule SNP_calling:
    input:
        fasta="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        bam="data/scaffolds_filtered/{sample}_sorted.bam",
        bam_bai="data/scaffolds_filtered/{sample}_sorted.bam.bai"
    output:
        fasta_fai="data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta.fai" % config["scaffold_minLen_filter"]["minlen"],
        unfilt_vcf="data/scaffolds_filtered/{sample}_unfiltered.vcf",
        filt_vcf="data/scaffolds_filtered/{sample}_filtered.vcf",
        zipped_vcf="data/scaffolds_filtered/{sample}_filtered.vcf.gz",
        zipped_vcf_index="data/scaffolds_filtered/{sample}_filtered.vcf.gz.tbi"
    conda:
        "../envs/scaffold_analyses.yaml"
    log:
        "logs/SNP_calling_{sample}.log"
    benchmark:
        "logs/benchmark/SNP_calling_{sample}.txt"
    threads: config["threads"]["SNP_calling"]
    params:
        max_cov=config["SNP_calling_params"]["max_cov"],
        minimum_AF=config["SNP_calling_params"]["minimum_AF"]
    shell:
        """
samtools faidx -o {output.fasta_fai} {input.fasta} > {log} 2>&1
lofreq call-parallel -d {params.max_cov} \
--no-default-filter \
--pp-threads {threads} \
-f {input.fasta} \
-o {output.unfilt_vcf} \
{input.bam} >> {log} 2>&1
lofreq filter -a {params.minimum_AF} \
-i {output.unfilt_vcf} \
-o {output.filt_vcf} >> {log} 2>&1
bgzip -c {output.filt_vcf} 2>> {log} 1> {output.zipped_vcf}
tabix -p vcf {output.zipped_vcf} >> {log} 2>&1
        """