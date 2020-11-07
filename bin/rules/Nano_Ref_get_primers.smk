rule Prepare_primers:
    input: primerfile
    output: 
        primers     = f"{datadir + prim}" + "primers.fasta",
        primers_5   = f"{datadir + prim}" + "primers_5.fasta",
        primers_3   = f"{datadir + prim}" + "primers_3.fasta"
    shell:
        """
cp {input} {output.primers}
sed '2~2s/^/X/' {output.primers} > {output.primers_5}
sed '2~2s/$/X/' {output.primers} > {output.primers_3}
        """ 