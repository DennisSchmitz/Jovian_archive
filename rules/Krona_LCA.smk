rule Krona_chart_and_LCA:
    input:
        classification="data/taxonomic_classification/{sample}.blastn",
        blast="data/taxonomic_classification/{sample}.blastn",
        stats="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats"
    output:
        taxtab="data/taxonomic_classification/{sample}.taxtab",
        taxMagtab="data/taxonomic_classification/{sample}.taxMagtab",
    conda:
        "../envs/Krona_plot.yaml"
    benchmark:
        "logs/benchmark/Krona_chart_and_LCA_{sample}.txt"
    threads: 1
    log:
        "logs/Krona_chart_and_LCA_{sample}.log"
    params:
        bitscoreDeltaLCA=config["taxonomic_classification_LCA"]["Krona_LCA"]["bitscoreDeltaLCA"],
        krona_tax_db=config["databases"]["Krona_taxonomy"]
    shell:  # We rm the [conda_path]/opt/krona/taxonomy folder and replace that to our specified krona_taxonomy path, updated weekly via crontab, see scripts Robert
        """
if [ ! -L $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g') ] # If symlink to Krona db does not exist...
then # Clean and make symlink to Krona db from the current Conda env (which has a unique and unpredictable hash, therefore, the which command)
    rm -rf $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
    ln -s {params.krona_tax_db} $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
fi
ktClassifyBLAST -o {output.taxtab} -t {params.bitscoreDeltaLCA} {input.classification} > {log} 2>&1
python bin/krona_magnitudes.py {output.taxtab} {input.blast} {input.stats} {output.taxMagtab} >> {log} 2>&1
        """