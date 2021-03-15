rule Prepare_primers:
    input: primerfile
    output: 
        primers     = f"{datadir + prim}" + "primers.fasta",
        primers_5   = f"{datadir + prim}" + "primers_5.fasta",
        primers_3   = f"{datadir + prim}" + "primers_3.fasta"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    threads: 1
    resources:
        memory = 4 * 1024
    shell:
        """
cp {input} {output.primers}
python bin/scripts/prepare_primers.py --primers {output.primers} --three {output.primers_3} --five {output.primers_5}
        """ 