
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
        all_scaffolds   =   f"{datadir + scf_filt}" + "{sample}/scaffolds.fasta",
        filt_scaffolds  =   f"{datadir + scf_filt}" + "{sample}_scaffolds_ge%snt.fasta" % config["Illumina_meta"]["minlen"],
    conda:
        f"{conda_envs}de_novo_assembly.yaml"
    log:
        f"{logdir}" + "De_novo_assembly_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "De_novo_assembly_{sample}.txt"
    threads: config["threads"]["De_novo_assembly"]
    params:
        max_GB_RAM          =   config["Illumina_meta"]["Spades"]["Max_gb_ram"],
        kmersizes           =   config["Illumina_meta"]["Spades"]["kmersizes"],
        minlength           =   config["Illumina_meta"]["minlen"],
        outputfoldername    =   f"{datadir + scf_raw}" + "{sample}/"
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