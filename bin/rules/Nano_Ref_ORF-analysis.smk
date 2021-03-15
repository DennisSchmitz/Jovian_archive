

rule ORF_Analysis:
    input: 
        ref =   rules.Index_ref.output.refcopy
    output:
        ORF_AA  =   f"{datadir + refdir + reference_basename}_ORF_AA.fa",
        ORF_NT  =   f"{datadir + refdir + reference_basename}_ORF_NT.fa",
        ORF_gff =   f"{datadir + refdir + reference_basename}_annotation.gff",
        gff_zip =   f"{datadir + refdir + reference_basename}_annotation.gff.gz",
        gff_ind =   f"{datadir + refdir + reference_basename}_annotation.gff.gz.tbi",
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}ORF_Analysis.log"
    benchmark:
        f"{logdir + bench}ORF_Analysis.txt"
    threads: 1
    resources:
        memory = 12
    params:
        procedure       =   config["Global"]["ORF_procedure"],
        output_format   =   config["Global"]["ORF_output_format"]
    shell:
        """
prodigal -q -i {input.ref} \
-a {output.ORF_AA} \
-d {output.ORF_NT} \
-o {output.ORF_gff} \
-p {params.procedure} \
-f {params.output_format} > {log} 2>&1
bgzip -c {output.ORF_gff} 2>> {log} 1> {output.gff_zip}
tabix -p gff {output.gff_zip} >> {log} 2>&1
        """ 