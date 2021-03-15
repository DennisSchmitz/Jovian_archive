rule Index_ref:
    input:
        ref =   reference
    output:
        refcopy         =   f"{datadir + refdir + reference_basename}.fasta",
        refcopy_index   =   f"{datadir + refdir + reference_basename}.fasta.1.bt2"
    conda:
        f"{conda_envs}Nano_ref_alignment.yaml"
    log:
        f"{logdir}Index_ref.log"
    benchmark:
        f"{logdir + bench}Index_ref.txt"
    threads: 4
    resources:
        memory = 8 * 1024
    shell:
        """
cat {input.ref} | seqkit replace -p "\-" -s -r "N" > {output.refcopy}
bowtie2-build --threads {threads} {output.refcopy} {output.refcopy} >> {log} 2>&1
        """