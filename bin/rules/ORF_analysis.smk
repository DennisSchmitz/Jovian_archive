
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule ORF_analysis:
    input:
        rules.De_novo_assembly.output.filt_scaffolds
    output:
        ORF_AA_fasta            =   f"{datadir + scf_filt}" + "{sample}_ORF_AA.fa",
        ORF_NT_fasta            =   f"{datadir + scf_filt}" + "{sample}_ORF_NT.fa",
        ORF_annotation_gff      =   f"{datadir + scf_filt}" + "{sample}_annotation.gff",
        zipped_gff3             =   f"{datadir + scf_filt}" + "{sample}_annotation.gff.gz",
        index_zipped_gff3       =   f"{datadir + scf_filt}" + "{sample}_annotation.gff.gz.tbi",
        contig_ORF_count_list   =   f"{datadir + scf_filt}" + "{sample}_contig_ORF_count_list.txt"
    conda:
        f"{conda_envs}Sequence_analysis.yaml"
    log:
        f"{logdir}" + "ORF_prediction_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "ORF_prediction_{sample}.txt"
    threads: 1
    resources:
        memory = 8 * 1024
    params:
        procedure       =   config["Global"]["ORF_procedure"],
        output_format   =   config["Global"]["ORF_output_format"]
    shell:
        """
prodigal -q -i {input} \
-a {output.ORF_AA_fasta} \
-d {output.ORF_NT_fasta} \
-o {output.ORF_annotation_gff} \
-p {params.procedure} \
-f {params.output_format} > {log} 2>&1
bgzip -c {output.ORF_annotation_gff} 2>> {log} 1> {output.zipped_gff3}
tabix -p gff {output.zipped_gff3} >> {log} 2>&1

egrep "^>" {output.ORF_NT_fasta} | sed 's/_/ /6' | tr -d ">" | cut -f 1 -d " " | uniq -c > {output.contig_ORF_count_list}
        """