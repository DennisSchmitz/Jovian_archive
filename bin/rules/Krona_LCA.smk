rule Krona_chart_and_LCA:
    input:
        classification  =   f"{datadir + taxclas}" + "{sample}.blastn",
        stats           =   f"{datadir + taxclas}" + "{sample}_perMinLenFiltScaffold.stats"
    output:
        taxtab      =   f"{datadir + taxclas}" + "{sample}.taxtab",
        taxMagtab   =   f"{datadir + taxclas}" + "{sample}.taxMagtab",
    conda:
        f"{conda_envs}Krona_plot.yaml"
    log:
        f"{logdir}" + "Krona_chart_and_LCA_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Krona_chart_and_LCA_{sample}.txt"
    threads: 1
    resources:
        memory = 8
    params:
        bitscoreDelta   =   config["Illumina_meta"]["LCA"]["bitscoreDelta"],
        krona_tax_db    =   config["databases"]["Krona_taxonomy"]
    shell:  # We rm the [conda_path]/opt/krona/taxonomy folder and replace that to our specified krona_taxonomy path, updated weekly via crontab, see scripts Robert
        """
if [ ! -L $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g') ] # If symlink to Krona db does not exist...
then # Clean and make symlink to Krona db from the current Conda env (which has a unique and unpredictable hash, therefore, the which command)
    rm -rf $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
    ln -s {params.krona_tax_db} $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
fi
ktClassifyBLAST -o {output.taxtab} -t {params.bitscoreDelta} {input.classification} > {log} 2>&1
python bin/scripts/krona_magnitudes.py {output.taxtab} {input.stats} {output.taxMagtab} >> {log} 2>&1
        """