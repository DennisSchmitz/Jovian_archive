#TODO kijken of dit multithreaded kan worden.
rule Illumina_extract_clean_consensus:
    input:
        bam             =   rules.Illumina_align_to_reference_it2.output.sorted_bam,
        raw_consensus   =   rules.Illumina_extract_raw_consensus_it2.output.raw_consensus_fasta, # Only needed for when there are no positions in the bed with a coverage of 0; in that case the IlluminaW fasta is actually suitable for downstream processes and it is simply copied.
    output:
        bedgraph                            =   f"{datadir + it2 + cons}" + "{sample}.bedgraph",
        filt_consensus_N_filt_ge_1          =   f"{res + seqs}" + "{sample}_N-filt_cov_ge_1.fa",
        filt_consensus_N_filt_ge_5          =   f"{res + seqs}" + "{sample}_N-filt_cov_ge_5.fa",
        filt_consensus_N_filt_ge_10         =   f"{res + seqs}" + "{sample}_N-filt_cov_ge_10.fa",
        filt_consensus_N_filt_ge_30         =   f"{res + seqs}" + "{sample}_N-filt_cov_ge_30.fa",
        filt_consensus_N_filt_ge_100        =   f"{res + seqs}" + "{sample}_N-filt_cov_ge_100.fa",
        filt_consensus_minus_filt_ge_1      =   f"{res + seqs}" + "{sample}_minus-filt_cov_ge_1.fa",
        filt_consensus_minus_filt_ge_5      =   f"{res + seqs}" + "{sample}_minus-filt_cov_ge_5.fa",
        filt_consensus_minus_filt_ge_10     =   f"{res + seqs}" + "{sample}_minus-filt_cov_ge_10.fa",
        filt_consensus_minus_filt_ge_30     =   f"{res + seqs}" + "{sample}_minus-filt_cov_ge_30.fa",
        filt_consensus_minus_filt_ge_100    =   f"{res + seqs}" + "{sample}_minus-filt_cov_ge_100.fa",
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}" + "Illumina_extract_clean_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_extract_clean_consensus_{sample}.txt"
    threads: 1
    params:
        output_data_folder      =   f"{datadir + it2 + cons}",
        output_results_folder   =   f"{res + seqs}"
    shell:
        """
bash bin/scripts/consensus_at_diff_coverages.sh {wildcards.sample} {input.bam} {input.raw_consensus} \
{params.output_data_folder} {params.output_results_folder} {log} >> {log} 2>&1
        """