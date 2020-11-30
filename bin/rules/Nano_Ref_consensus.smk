rule consensus_cov_1:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_1.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_1.fa"
    params:
        coverage    =   "1"
    conda:
        f"{conda_envs}Nano_ref_consensus.yaml"
    log:
        f"{logdir}" + "Generate_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Generate_consensus_{sample}.txt"
    threads: 2
    shell:
        """
python bin/scripts/Consensus.py \
--input {input.bam} \
--reference {input.ref} \
--gff {input.gff} \
--name {wildcards.sample} \
--mincov {params.coverage} \
--consensus {output.cons} \
--gapcorrected {output.gapcor} \
--threads {threads}
        """


rule consensus_cov_5:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_5.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_5.fa"
    params:
        coverage    =   "5"
    conda:
        f"{conda_envs}Nano_ref_consensus.yaml"
    log:
        f"{logdir}" + "Generate_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Generate_consensus_{sample}.txt"
    threads: 2
    shell:
        """
python bin/scripts/Consensus.py \
--input {input.bam} \
--reference {input.ref} \
--gff {input.gff} \
--name {wildcards.sample} \
--mincov {params.coverage} \
--consensus {output.cons} \
--gapcorrected {output.gapcor} \
--threads {threads}
        """

rule consensus_cov_10:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_10.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_10.fa"
    params:
        coverage    =   "10"
    conda:
        f"{conda_envs}Nano_ref_consensus.yaml"
    log:
        f"{logdir}" + "Generate_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Generate_consensus_{sample}.txt"
    threads: 2
    shell:
        """
python bin/scripts/Consensus.py \
--input {input.bam} \
--reference {input.ref} \
--gff {input.gff} \
--name {wildcards.sample} \
--mincov {params.coverage} \
--consensus {output.cons} \
--gapcorrected {output.gapcor} \
--threads {threads}
        """


rule consensus_cov_30:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_30.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_30.fa"
    params:
        coverage    =   30
    conda:
        f"{conda_envs}Nano_ref_consensus.yaml"
    log:
        f"{logdir}" + "Generate_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Generate_consensus_{sample}.txt"
    threads: 2
    shell:
        """
python bin/scripts/Consensus.py \
--input {input.bam} \
--reference {input.ref} \
--gff {input.gff} \
--name {wildcards.sample} \
--mincov {params.coverage} \
--consensus {output.cons} \
--gapcorrected {output.gapcor} \
--threads {threads}
        """

rule consensus_cov_100:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_100.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_100.fa"
    params:
        coverage    =   "100"
    conda:
        f"{conda_envs}Nano_ref_consensus.yaml"
    log:
        f"{logdir}" + "Generate_consensus_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Generate_consensus_{sample}.txt"
    threads: 2
    shell:
        """
python bin/scripts/Consensus.py \
--input {input.bam} \
--reference {input.ref} \
--gff {input.gff} \
--name {wildcards.sample} \
--mincov {params.coverage} \
--consensus {output.cons} \
--gapcorrected {output.gapcor} \
--threads {threads}
        """