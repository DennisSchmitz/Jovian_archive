#>##################################################################################################
#>#### Align to ref, mark and optionally remove duplicates, call SNPs, generate new consensus  #####
#>##################################################################################################
rule Illumina_align_to_reference_it2:
    input:
        reads       =   rules.RemovePrimers_pt2.output,
        reference   =   rules.Illumina_extract_raw_consensus_it1.output.reference_copy_it2
    output:
        reference_index     =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}.fasta.1.bt2", # I've only specified ".fasta.1.bt2", but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated. #TODO find a way to specify all output correctly (multiext snakemake syntax?)
        sorted_bam          =   f"{datadir + it2 + aln}" + "{sample}_sorted.bam",
        sorted_bam_index    =   f"{datadir + it2 + aln}" + "{sample}_sorted.bam.bai",
        dup_metrics         =   f"{datadir + it2 + aln}" + "{sample}_sorted.MarkDup_metrics", #TODO deze toevoegen aan MultiQC?
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_align_to_reference_it2_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_align_to_reference_it2_{sample}.txt"
    threads: config["threads"]["Illumina_align_to_reference"]
    resources: 
        memory = config["threads"]["Illumina_align_to_reference"] * 12
    params:
        aln_type        =   config["Illumina_ref"]["Alignment"]["Alignment_type"],
        remove_dups     =   config["Illumina_ref"]["Alignment"]["Duplicates"], #! Don't change this, see this gotcha with duplicate marked reads in bedtools genomecov (which is used downstream): https://groups.google.com/forum/#!msg/bedtools-discuss/wJNC2-icIb4/wflT6PnEHQAJ . bedtools genomecov is not able to filter them out and includes those dup-reads in it's coverage metrics. So the downstream BoC analysis and consensus at diff cov processes require dups to be HARD removed.
        markdup_mode    =   config["Illumina_ref"]["Alignment"]["Duplicate_marking"],
        max_read_length =   config["Illumina_ref"]["Alignment"]["Max_read_length"] # This is the default value and also the max read length of Illumina in-house sequencing.
    shell: # LoFreq dindel req for indel calling
        """
bowtie2-build --threads {threads} {input.reference} {input.reference} >> {log} 2>&1
bowtie2 --time --threads {threads} {params.aln_type} \
-x {input.reference} \
-U {input.reads} 2> {log} |\
samtools view -@ {threads} -uS - 2>> {log} |\
samtools collate -@ {threads} -O - 2>> {log} |\
samtools fixmate -@ {threads} -m - - 2>> {log} |\
samtools sort -@ {threads} - -o - 2>> {log} |\
lofreq indelqual --dindel -f {input.reference} -o - - 2>> {log} |\
samtools markdup -@ {threads} -l {params.max_read_length} -m {params.markdup_mode} {params.remove_dups} -f {output.dup_metrics} - {output.sorted_bam} >> {log} 2>&1
samtools index -@ {threads} {output.sorted_bam} >> {log} 2>&1
        """