
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule De_novo_assembly:
    input:
        fastq_pR1       =   rules.HuGo_removal_pt2_extract_paired_unmapped_reads.output.fastq_R1,
        fastq_pR2       =   rules.HuGo_removal_pt2_extract_paired_unmapped_reads.output.fastq_R2,
        fastq_unpaired  =   rules.HuGo_removal_pt3_extract_unpaired_unmapped_reads.output
    output:
        all_scaffolds   =   "data/scaffolds_raw/{sample}/scaffolds.fasta",
        filt_scaffolds  =   "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["Illumina_meta"]["minlen"],
    conda:
        conda_envs + "de_novo_assembly.yaml"
    benchmark:
        "logs/benchmark/De_novo_assembly_{sample}.txt"
    threads: config["threads"]["De_novo_assembly"]
    log:
        "logs/De_novo_assembly_{sample}.log"
    params:
        max_GB_RAM          =   "100",
        kmersizes           =   config["Illumina_meta"]["Spades"]["kmersizes"],
        outputfoldername    =   "data/scaffolds_raw/{sample}/",
        minlength           =   config["Illumina_meta"]["minlen"],
    shell:
        """
spades.py --only-assembler --meta \
-1 {input.fastq_pR1} \
-2 {input.fastq_pR2} \
-s {input.fastq_unpaired} \
-t {threads} \
-m {params.max_GB_RAM} \
-k {params.kmersizes} \
-o {params.outputfoldername} > {log} 2>&1
seqtk seq {output.all_scaffolds} 2>> {log} |\
gawk -F "_" '/^>/ {{if ($4 >= {params.minlength}) {{print $0; getline; print $0}};}}' 2>> {log} 1> {output.filt_scaffolds} 
        """