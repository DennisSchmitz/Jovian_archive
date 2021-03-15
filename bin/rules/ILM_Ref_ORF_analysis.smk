rule Illumina_reference_ORF_analysis:
    input:
        reference   =   rules.Illumina_extract_raw_consensus_it1.output.reference_copy_it2
    output: 
        ORF_AA_fasta        =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}_ORF_AA.fa",
        ORF_NT_fasta        =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}_ORF_NT.fa",
        ORF_annotation_gff  =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}_annotation.gff",
        zipped_gff3         =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}_annotation.gff.gz",
        index_zipped_gff3   =   f"{datadir + it2 + refdir + reference_basename}" + "_{sample}_annotation.gff.gz.tbi",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}" + "Illumina_reference_ORF_analysis_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_reference_ORF_analysis_{sample}.txt"
    threads: 1
    resources:
        memory = 8 * 1024
    params: #? Currently it's using the same prodigal settings as the main workflow, I see no problems with it since it's both foremost intended for viruses.
        procedure       =   config["Global"]["ORF_procedure"],
        output_format   =   config["Global"]["ORF_output_format"]
    shell:
        """
prodigal -q -i {input.reference} \
-a {output.ORF_AA_fasta} \
-d {output.ORF_NT_fasta} \
-o {output.ORF_annotation_gff} \
-p {params.procedure} \
-f {params.output_format} > {log} 2>&1
bgzip -c {output.ORF_annotation_gff} 2>> {log} 1> {output.zipped_gff3}
tabix -p gff {output.zipped_gff3} >> {log} 2>&1
        """