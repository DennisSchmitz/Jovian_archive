
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule ORF_analysis:
    input:
        "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["Illumina_meta"]["minlen"]
    output:
        ORF_AA_fasta="data/scaffolds_filtered/{sample}_ORF_AA.fa",
        ORF_NT_fasta="data/scaffolds_filtered/{sample}_ORF_NT.fa",
        ORF_annotation_gff="data/scaffolds_filtered/{sample}_annotation.gff",
        zipped_gff3="data/scaffolds_filtered/{sample}_annotation.gff.gz",
        index_zipped_gff3="data/scaffolds_filtered/{sample}_annotation.gff.gz.tbi",
        contig_ORF_count_list="data/scaffolds_filtered/{sample}_contig_ORF_count_list.txt"
    conda:
        "../envs/Sequence_analysis.yaml"
    log:
        "logs/ORF_prediction_{sample}.log"
    benchmark:
        "logs/benchmark/ORF_prediction_{sample}.txt"
    threads: 1
    params:
        procedure=config["Global"]["ORF_procedure"],
        output_format=config["Global"]["ORF_output_format"]
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