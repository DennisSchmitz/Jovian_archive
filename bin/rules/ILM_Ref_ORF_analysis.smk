
# Nuttig voor IGVjs vis. Gejat uit Jovian core met minor changes, kunnen we waarschijlijk efficienter doen. Bijvoorbeeld door gewoon een goed gecureerde ORF annotatie toe te voegen bij starten van analyse.
rule Illumina_reference_ORF_analysis:
    input:
        reference   =   rules.Illumina_index_reference.output.reference_copy
    output: 
        ORF_AA_fasta        =   f"{datadir + refdir + reference_basename}_ORF_AA.fa",
        ORF_NT_fasta        =   f"{datadir + refdir + reference_basename}_ORF_NT.fa",
        ORF_annotation_gff  =   f"{datadir + refdir + reference_basename}_annotation.gff",
        zipped_gff3         =   f"{datadir + refdir + reference_basename}_annotation.gff.gz",
        index_zipped_gff3   =   f"{datadir + refdir + reference_basename}_annotation.gff.gz.tbi",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}Illumina_reference_ORF_analysis.log"
    benchmark:
        f"{logdir + bench}Illumina_reference_ORF_analysis.txt"
    threads: 1
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