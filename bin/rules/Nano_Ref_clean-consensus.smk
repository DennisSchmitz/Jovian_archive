


rule extract_cleaned_consensus:
    input:
        bam             =   rules.Align_to_reference_pt1.output.bam
    output:
        consensus_1     =   f"{res + seqs}" + "{sample}_cov_ge_1.fa",
        consensus_5     =   f"{res + seqs}" + "{sample}_cov_ge_5.fa",
        consensus_10    =   f"{res + seqs}" + "{sample}_cov_ge_10.fa",
        consensus_30    =   f"{res + seqs}" + "{sample}_cov_ge_30.fa",
        consensus_100   =   f"{res + seqs}" + "{sample}_cov_ge_100.fa"
    params:
        output_results_folder   =   f"{res + seqs}",
        coverages   =   [ '1', '5', '10', '30', '100' ],
        filename    =   "{sample}_cov_ge_"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}" + "Extract_cleaned_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Extract_cleaned_consensus_{sample}.txt"
    threads: 1
    shell:
        """
python bin/scripts/ConsensusFromBam.py -input {input.bam} -output {params.output_results_folder}{params.filename}{params.coverages} -cov {params.coverages} -name consensus_{params.filename}{params.coverages}
        """