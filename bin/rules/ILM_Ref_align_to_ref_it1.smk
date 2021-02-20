#>##################################################################################################
#>#### Align to ref, mark and optionally remove duplicates, call SNPs, generate new consensus  #####
#>##################################################################################################
rule Illumina_align_to_reference_it1:
    input:
        reads       =   rules.RemovePrimers_pt2.output,
        reference   =   rules.Illumina_index_reference.output.reference_copy
    output:
        sorted_bam          =   f"{datadir + it1 + aln}" + "{sample}_sorted.bam",
        sorted_bam_index    =   f"{datadir + it1 + aln}" + "{sample}_sorted.bam.bai",
        dup_metrics         =   f"{datadir + it1 + aln}" + "{sample}_sorted.MarkDup_metrics" #TODO deze toevoegen aan MultiQC?
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_align_to_reference_it1_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_align_to_reference_it1_{sample}.txt"
    threads: config["threads"]["Illumina_align_to_reference"]
    params:
        aln_type        =   config["Illumina_ref"]["Alignment"]["Alignment_type"],
        remove_dups     =   config["Illumina_ref"]["Alignment"]["Duplicates"], #! Don't change this, see this gotcha with duplicate marked reads in bedtools genomecov (which is used downstream): https://groups.google.com/forum/#!msg/bedtools-discuss/wJNC2-icIb4/wflT6PnEHQAJ . bedtools genomecov is not able to filter them out and includes those dup-reads in it's coverage metrics. So the downstream BoC analysis and consensus at diff cov processes require dups to be HARD removed.
        markdup_mode    =   config["Illumina_ref"]["Alignment"]["Duplicate_marking"],
        max_read_length =   config["Illumina_ref"]["Alignment"]["Max_read_length"] # This is the default value and also the max read length of Illumina in-house sequencing.
    shell: # LoFreq dindel is required for indel calling downstream
        """
bowtie2 --time --quiet --threads {threads} {params.aln_type} \
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