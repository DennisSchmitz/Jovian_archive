


rule extract_cleaned_consensus:
    input:
        raw_consensus = rules.Create_raw_consensus.output.consensus,
        bam = rules.Align_to_reference_pt1.output.bam,
    output:
        bedgraph                            = datadir + cons + filt + "{sample}.bedgraph",
        filt_consensus_N_filt_ge_1          = res + seqs + "{sample}_N-filt_cov_ge_1.fa",
        filt_consensus_N_filt_ge_5          = res + seqs + "{sample}_N-filt_cov_ge_5.fa",
        filt_consensus_N_filt_ge_10         = res + seqs + "{sample}_N-filt_cov_ge_10.fa",
        filt_consensus_N_filt_ge_30         = res + seqs + "{sample}_N-filt_cov_ge_30.fa",
        filt_consensus_N_filt_ge_100        = res + seqs + "{sample}_N-filt_cov_ge_100.fa",
        filt_consensus_minus_filt_ge_1      = res + seqs + "{sample}_minus-filt_cov_ge_1.fa",
        filt_consensus_minus_filt_ge_5      = res + seqs + "{sample}_minus-filt_cov_ge_5.fa",
        filt_consensus_minus_filt_ge_10     = res + seqs + "{sample}_minus-filt_cov_ge_10.fa",
        filt_consensus_minus_filt_ge_30     = res + seqs + "{sample}_minus-filt_cov_ge_30.fa",
        filt_consensus_minus_filt_ge_100    = res + seqs + "{sample}_minus-filt_cov_ge_100.fa",
    params:
        output_data_folder = datadir + cons + filt,
        output_results_folder = datadir + cons + seqs
    conda:
        conda_envs + "Nano_ref_alignment.yaml"
    threads: 26
    log:
        logdir + "Extract_cleaned_consensus_{sample}.log"
    benchmark:
        logdir + bench + "Extract_cleaned_consensus_{sample}.txt"
    shell:
        """
bash bin/scripts/RA_consensus_at_diff_coverages.sh {wildcards.sample} {input.bam} {input.raw_consensus} \
{params.output_data_folder} {params.output_results_folder} {log} >> {log} 2>&1
        """