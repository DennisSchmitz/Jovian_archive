rule Cut_primers:
    input:
        fastq       =   rules.Cleanup.output.qc_fastq,
        primers     =   rules.Prepare_primers.output.primers,
        ref         =   rules.Index_ref.output.refcopy
    output: f"{datadir + cln + prdir}" + "{sample}.fastq"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    log:
        f"{logdir}" + "Primer_removal_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Primer_removal_{sample}.txt"
    threads: config["threads"]["Nanopore_primer_removal"]
    resources: 
        memory = (config["threads"]["Nanopore_primer_removal"] * 12) * 1024
    shell:
        """
python bin/scripts/RemoveONTPrimers.py --input {input.fastq} --reference {input.ref} --primers {input.primers} --threads {threads} --output {output}
sed -i "N;N;N;/\\n\\n/d" {output}
        """