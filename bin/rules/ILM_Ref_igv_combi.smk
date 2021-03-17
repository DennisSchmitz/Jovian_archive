rule Illumina_HTML_IGVJs_generate_final:
    input:
        expand( "{p}{chunk_name}_{sample}",
                p           = f"{datadir + html}",
                chunk_name  = [ '2_tab',
                                '4_html_divs',
                                '6_js_flex'
                                ],
                sample      = SAMPLES
                )
    output:
        f"{res}igv_ilr.html"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}Illumina_HTML_IGVJs_generate_final.log"
    benchmark:
        f"{logdir + bench}Illumina_HTML_IGVJs_generate_final.txt"
    threads: 1
    resources:
        memory = 4 * 1024
    params:
        chunkpath       =   f"{fls + chunks}", 
        tab_basename    =   f"{datadir + html}2_tab_",
        div_basename    =   f"{datadir + html}4_html_divs_",
        js_flex_output  =   f"{datadir + html}6_js_flex_"
    shell:
        """
cat {params.chunkpath}1_header.html > {output}
cat {params.tab_basename}* >> {output}
cat {params.chunkpath}3_tab_explanation_Illumina.html >> {output}
cat {params.div_basename}* >> {output}
cat {params.chunkpath}5_js_begin.html >> {output}
cat {params.js_flex_output}* >> {output}
cat {params.chunkpath}7_js_end.html >> {output}
        """