
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule Krona_chart_combine:
    input:
        sorted(
            expand( rules.lca_mgkit.output.taxMagtab,
                    sample  =   set(SAMPLES)
                    )
            )
    output:
        f"{res}krona.html"
    conda:
        f"{conda_envs}Krona_plot.yaml"
    log:
        f"{logdir}Krona_chart_combine.log"
    benchmark:
        f"{logdir + bench}Krona_chart_combine.txt"
    threads: 1
    resources:
        memory = 8 * 1024
    params:
        krona_tax_db    =   config["databases"]["Krona_taxonomy"]
    shell:
        """
if [ ! -L $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g') ] # If symlink to Krona db does not exist...
then # Clean and make symlink to Krona db from the current Conda env (which has a unique and unpredictable hash, therefore, the which command)
    rm -rf $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
    ln -s {params.krona_tax_db} $(which ktClassifyBLAST | sed 's|/bin/ktClassifyBLAST|/opt/krona/taxonomy|g')
fi
ktImportTaxonomy {input} -i -k -m 4 -o {output} > {log} 2>&1
        """