rule Prepare_primers:
    input: primerfile
    output: 
        primers     = f"{datadir + prdir}" + "primers.fasta",
        primers_5   = f"{datadir + prdir}" + "primers_5.fasta",
        primers_3   = f"{datadir + prdir}" + "primers_3.fasta"
    shell:
        """
cp {input} {output.primers}
sed '2~2s/^/X/' {output.primers} > {output.primers_5}
sed -i '2~2s/$/X/' {output.primers} > {output.primers_3}
        """ 