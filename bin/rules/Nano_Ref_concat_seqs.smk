


rule concat_seqs:
    input: 
        cov_1   =   expand( "{p}{sample}_cov_ge_1.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_5   =   expand( "{p}{sample}_cov_ge_5.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_10  =   expand( "{p}{sample}_cov_ge_10.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_30  =   expand( "{p}{sample}_cov_ge_30.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_100 =   expand( "{p}{sample}_cov_ge_100.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            )
    output:
        conc_1      =   f"{res + seqs}" + "concat_cov_ge_1.fasta",
        conc_5      =   f"{res + seqs}" + "concat_cov_ge_5.fasta",
        conc_10     =   f"{res + seqs}" + "concat_cov_ge_10.fasta",
        conc_30     =   f"{res + seqs}" + "concat_cov_ge_30.fasta",
        conc_100    =   f"{res + seqs}" + "concat_cov_ge_100.fasta"
    shell:
        """
cat {input.cov_1} >> {output.conc_1}
cat {input.cov_5} >> {output.conc_5}
cat {input.cov_10} >> {output.conc_10}
cat {input.cov_30} >> {output.conc_30}
cat {input.cov_100} >> {output.conc_100}
        """