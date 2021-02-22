

rule Illumina_index_reference:
    input:
        reference   =   f"{reference}"
    output:
        reference_copy      =   f"{datadir + it1 + refdir + reference_basename}.fasta"
    conda:
        f"{conda_envs}Illumina_ref_alignment.yaml"
    log:
        f"{logdir}Illumina_index_reference.log"
    benchmark:
        f"{logdir + bench}Illumina_index_reference.txt"
    threads: 4
    shell: # The reference is copied to the hardcoded subdir to make it standardized and easily logged. Convert it to a two-line fasta for easier downstream processing, replace ambiguity nucleotides with an N to avoid errors downstream.
        """
cat {input.reference} | seqtk seq - | sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D]/N/g' > {output.reference_copy} 2>> {log}
        """