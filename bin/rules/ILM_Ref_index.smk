

rule Illumina_index_reference:
    input:
        reference   =   f"{reference}"
    output:
        reference_copy      =   f"{datadir + refdir + reference_basename}.fasta",
        reference_index     =   f"{datadir + refdir + reference_basename}.fasta.1.bt2", # I've only specified ".fasta.1.bt2", but the "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2" and "rev.2.bt2" are implicitly generated. #TODO find a way to specify all output correctly (multiext snakemake syntax?)
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}Illumina_index_reference.log"
    benchmark:
        f"{logdir + bench}Illumina_index_reference.txt"
    threads: 4
    shell: # The reference is copied to the hardcoded subdir to make it standardized and easily logged. Convert it to a two-line fasta for easier downstream processing.
        """
cat {input.reference} | seqtk seq - > {output.reference_copy}
bowtie2-build --threads {threads} {output.reference_copy} {output.reference_copy} >> {log} 2>&1
        """