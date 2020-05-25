rule Illumina_extract_raw_consensus:
    input:
        bam         =   rules.Illumina_align_to_reference.output.sorted_bam,
        reference   =   rules.Illumina_index_reference.output.reference_copy,
    output: #TODO check if it can use bcf output instead of vcf for downstream processing, saves diskspace, but not a huge differences for small viral genomes. Would require changes in the igvjs index
        gzipped_vcf         = f"{datadir + cons + raw}" + "{sample}_calls.vcf.gz",
        raw_consensus_fasta = f"{datadir + cons + raw}" + "{sample}_raw_consensus.fa",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_extract_raw_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_extract_raw_consensus_{sample}.txt"
    threads: 1 # Increasing this makes no differences for monopartite references/viruses. I think it splits different chromosomes to different threads #TODO check this in future version that is compatible with segmented viruses.
    params: #TODO move this param to pipeline_variables.yaml when we assess and optimize this value for different viral families.
        calling_prior   =   "1.1e-3" # From manual: mutation rate (use bigger for greater sensitivity), use with -m [1.1e-3]. Also see https://samtools.github.io/bcftools/howtos/variant-calling.html --> higher value is less strict and vice versa #TODO this can be (has to be?) adapted to the different virus mutation rate. Assess later and optimize for different viral families
    shell:
        """
bcftools mpileup --threads {threads} --ignore-RG -O u -d 10000 -f {input.reference} {input.bam} 2>> {log} |\
bcftools call --threads {threads} -m --prior {params.calling_prior} --ploidy 1 -mv -O z 2>> {log} |\
bcftools norm -m -both -O z -f {input.reference} -o {output.gzipped_vcf} >> {log} 2>&1
tabix {output.gzipped_vcf} >> {log} 2>&1
cat {input.reference} 2>> {log} |\
bcftools consensus {output.gzipped_vcf} | seqtk seq - > {output.raw_consensus_fasta} 2>> {log}
        """