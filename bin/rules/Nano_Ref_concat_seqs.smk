rule concat_seqs_standard:
    input: 
        cov_1   =   expand( "{p}{sample}_standard_cov_ge_1.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_5   =   expand( "{p}{sample}_standard_cov_ge_5.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_10  =   expand( "{p}{sample}_standard_cov_ge_10.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_30  =   expand( "{p}{sample}_standard_cov_ge_30.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_100 =   expand( "{p}{sample}_standard_cov_ge_100.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            )
    output:
        conc_1      =   f"{res + seqs}" + "concat_standard_cov_ge_1.fasta",
        conc_5      =   f"{res + seqs}" + "concat_standard_cov_ge_5.fasta",
        conc_10     =   f"{res + seqs}" + "concat_standard_cov_ge_10.fasta",
        conc_30     =   f"{res + seqs}" + "concat_standard_cov_ge_30.fasta",
        conc_100    =   f"{res + seqs}" + "concat_standard_cov_ge_100.fasta"
    params:
        wildcard_1      =   f"{res + seqs}" + "*_standard_cov_ge_1.fa",
        wildcard_5      =   f"{res + seqs}" + "*_standard_cov_ge_5.fa",
        wildcard_10     =   f"{res + seqs}" + "*_standard_cov_ge_10.fa",
        wildcard_30     =   f"{res + seqs}" + "*_standard_cov_ge_30.fa",
        wildcard_100    =   f"{res + seqs}" + "*_standard_cov_ge_100.fa"
    threads: 1
    resources:
        memory = 6 * 1024
    shell:
        """
cat {params.wildcard_1} >> {output.conc_1}
cat {params.wildcard_5} >> {output.conc_5}
cat {params.wildcard_10} >> {output.conc_10}
cat {params.wildcard_30} >> {output.conc_30}
cat {params.wildcard_100} >> {output.conc_100}
        """


rule concat_seqs_gapcorrected:
    input: 
        cov_1   =   expand( "{p}{sample}_gap_corrected_cov_ge_1.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_5   =   expand( "{p}{sample}_gap_corrected_cov_ge_5.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_10  =   expand( "{p}{sample}_gap_corrected_cov_ge_10.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_30  =   expand( "{p}{sample}_gap_corrected_cov_ge_30.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            ),
        cov_100 =   expand( "{p}{sample}_gap_corrected_cov_ge_100.fa",
                            p       =   f"{res + seqs}",
                            sample  =   SAMPLES
                            )
    output:
        conc_1      =   f"{res + seqs}" + "concat_gap_corrected_cov_ge_1.fasta",
        conc_5      =   f"{res + seqs}" + "concat_gap_corrected_cov_ge_5.fasta",
        conc_10     =   f"{res + seqs}" + "concat_gap_corrected_cov_ge_10.fasta",
        conc_30     =   f"{res + seqs}" + "concat_gap_corrected_cov_ge_30.fasta",
        conc_100    =   f"{res + seqs}" + "concat_gap_corrected_cov_ge_100.fasta"
    params:
        wildcard_1      =   f"{res + seqs}" + "*_gap_corrected_cov_ge_1.fa",
        wildcard_5      =   f"{res + seqs}" + "*_gap_corrected_cov_ge_5.fa",
        wildcard_10     =   f"{res + seqs}" + "*_gap_corrected_cov_ge_10.fa",
        wildcard_30     =   f"{res + seqs}" + "*_gap_corrected_cov_ge_30.fa",
        wildcard_100    =   f"{res + seqs}" + "*_gap_corrected_cov_ge_100.fa"
    threads: 1
    resources:
        memory = 6 * 1024
    shell:
        """
cat {params.wildcard_1} >> {output.conc_1}
cat {params.wildcard_5} >> {output.conc_5}
cat {params.wildcard_10} >> {output.conc_10}
cat {params.wildcard_30} >> {output.conc_30}
cat {params.wildcard_100} >> {output.conc_100}
        """