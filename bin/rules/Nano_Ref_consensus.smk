#! the {output.cov} is redundant, only need one rule to make it. Don't have time to fix it now. But it's stupid and redundant, apologies.
rule consensus_cov_1:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_1.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_1.fa",
        cov     =   f"{res + covs}" + "{sample}_coverage1.tsv",
        ins     =   f"{res + insr}" + "{sample}_inserts_1.tsv",
        vcf     =   f"{res + vf}" + "{sample}_1.vcf",
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
--coverage {output.cov} \
--threads {threads} \
--vcf {output.vcf} \
--insertions {output.ins} 2>> {log}
        """


rule consensus_cov_5:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_5.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_5.fa",
        cov     =   f"{res + covs}" + "{sample}_coverage5.tsv",
        ins     =   f"{res + insr}" + "{sample}_inserts_5.tsv",
        vcf     =   f"{res + vf}" + "{sample}_5.vcf",
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
--coverage {output.cov} \
--threads {threads} \
--vcf {output.vcf} \
--insertions {output.ins} 2>> {log}
        """

rule consensus_cov_10:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_10.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_10.fa",
        cov     =   f"{res + covs}" + "{sample}_coverage10.tsv",
        ins     =   f"{res + insr}" + "{sample}_inserts_10.tsv",
        vcf     =   f"{res + vf}" + "{sample}_10.vcf",
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
--coverage {output.cov} \
--threads {threads} \
--vcf {output.vcf} \
--insertions {output.ins} 2>> {log}
        """


rule consensus_cov_30:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_30.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_30.fa",
        cov     =   f"{res + covs}" + "{sample}_coverage30.tsv",
        ins     =   f"{res + insr}" + "{sample}_inserts_30.tsv",
        vcf     =   f"{res + vf}" + "{sample}_30.vcf",
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
--coverage {output.cov} \
--threads {threads} \
--vcf {output.vcf} \
--insertions {output.ins} 2>> {log}
        """

rule consensus_cov_100:
    input:
        bam     =   rules.Align_to_reference_pt1.output.bam,
        ref     =   rules.Index_ref.output.refcopy,
        gff     =   rules.ORF_Analysis.output.ORF_gff
    output:
        cons    =   f"{res + seqs}" + "{sample}_standard_cov_ge_100.fa",
        gapcor  =   f"{res + seqs}" + "{sample}_gap_corrected_cov_ge_100.fa",
        cov     =   f"{res + covs}" + "{sample}_coverage100.tsv",
        ins     =   f"{res + insr}" + "{sample}_inserts_100.tsv",
        vcf     =   f"{res + vf}" + "{sample}_100.vcf",
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
--coverage {output.cov} \
--threads {threads} \
--vcf {output.vcf} \
--insertions {output.ins} 2>> {log}
        """