


rule HTML_IGVJs_variable_parts:
    input:
        ref     =   rules.Index_ref.output.refcopy,
        GC_bed  =   rules.determine_GC_content.output.GC_bed,
        ORF_gff =   rules.ORF_Analysis.output.gff_zip,
        vcf     =   rules.Align_to_reference_pt2.output.vcf,
        bam     =   rules.Align_to_reference_pt1.output.bam
    output:
        tab =   f"{datadir + chunks}" + "2_tab_{sample}",
        div =   f"{datadir + chunks}" + "4_html_divs_{sample}",
        js  =   f"{datadir + chunks}" + "6_js_flex_{sample}"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}" + "HTML_IGVJs_variable_parts_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_tabs.sh {wildcards.sample} {output.tab}

bash bin/html/igvjs_write_divs.sh {wildcards.sample} {output.div}

bash bin/html/nano_igvjs_write_flex_js_middle.sh {wildcards.sample} {output.js} \
{input.ref} {input.GC_bed} {input.ORF_gff} \
{input.vcf} {input.bam}
        """

rule HTML_IGVJs_generate_file:
    input:
        expand( "{p}{b}_{sample}",
                p       =   f"{datadir + chunks}",
                b       =   [   '2_tab',
                                '4_html_divs',
                                '6_js_flex'
                                ],
                sample  =   SAMPLES
                )
    output:
        html    =   f"{res}igv.html"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}HTML_IGVJs_generate_file.log"
    benchmark:
        f"{logdir + bench}HTML_IGVJs_generate_file.txt"
    threads: 1
    params:
        tab_basename    =   f"{datadir + chunks}2_tab_",
        div_basename    =   f"{datadir + chunks}4_html_divs_",
        js_flex_output  =   f"{datadir + chunks}6_js_flex_",
    shell:
        """
cat files/html_chunks/1_header.html > {output.html}
cat {params.tab_basename}* >> {output.html}
cat files/html_chunks/3_tab_explanation_Nano.html >> {output.html}
cat {params.div_basename}* >> {output.html}
cat files/html_chunks/5_js_begin.html >> {output.html}
cat {params.js_flex_output}* >> {output.html}
cat files/html_chunks/7_js_end.html >> {output.html}
        """ 