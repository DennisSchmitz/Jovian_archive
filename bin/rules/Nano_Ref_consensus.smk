rule consensus_cov_1:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam
    output:
        cons    =   f"{res + seqs}" + "{sample}_cov_ge_1.fa",
    params:
        output_results_folder   =   f"{res + seqs}",
        coverage    =   "1",
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
python bin/scripts/ConsensusFromBam.py -input {input.bam} -output {output.cons} -cov {params.coverage} -name consensus_{params.filename}{params.coverage}
        """


rule consensus_cov_5:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam
    output:
        cons    =   f"{res + seqs}" + "{sample}_cov_ge_5.fa",
    params:
        output_results_folder   =   f"{res + seqs}",
        coverage    =   "5",
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
python bin/scripts/ConsensusFromBam.py -input {input.bam} -output {output.cons} -cov {params.coverage} -name consensus_{params.filename}{params.coverage}
        """

rule consensus_cov_10:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam
    output:
        cons    =   f"{res + seqs}" + "{sample}_cov_ge_10.fa",
    params:
        output_results_folder   =   f"{res + seqs}",
        coverage    =   "10",
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
python bin/scripts/ConsensusFromBam.py -input {input.bam} -output {output.cons} -cov {params.coverage} -name consensus_{params.filename}{params.coverage}
        """


rule consensus_cov_30:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam
    output:
        cons    =   f"{res + seqs}" + "{sample}_cov_ge_30.fa",
    params:
        output_results_folder   =   f"{res + seqs}",
        coverage    =   "30",
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
python bin/scripts/ConsensusFromBam.py -input {input.bam} -output {output.cons} -cov {params.coverage} -name consensus_{params.filename}{params.coverage}
        """

rule consensus_cov_100:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam
    output:
        cons    =   f"{res + seqs}" + "{sample}_cov_ge_100.fa",
    params:
        output_results_folder   =   f"{res + seqs}",
        coverage    =   "100",
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
python bin/scripts/ConsensusFromBam.py -input {input.bam} -output {output.cons} -cov {params.coverage} -name consensus_{params.filename}{params.coverage}
        """