
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


rule HTML_IGVJs_variable_parts:
    input:
        fasta                   =   rules.De_novo_assembly.output.filt_scaffolds,
        ref_GC_bedgraph         =   rules.Determine_GC_content.output.GC_bed,
        ref_zipped_ORF_gff      =   rules.ORF_analysis.output.zipped_gff3,
        basepath_zipped_SNP_vcf =   rules.SNP_calling.output.zipped_vcf,
        basepath_sorted_bam     =   rules.Read2scaffold_alignment_with_rmDup_and_fraglength.output.bam
    output:
        tab_output      =   f"{datadir + html}" + "2_tab_{sample}",
        div_output      =   f"{datadir + html}" + "4_html_divs_{sample}",
        js_flex_output  =   f"{datadir + html}" + "6_js_flex_{sample}",
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}" + "HTML_IGVJs_variable_parts_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_tabs.sh {wildcards.sample} {output.tab_output}

bash bin/html/igvjs_write_divs.sh {wildcards.sample} {output.div_output}

bash bin/html/meta_igvjs_write_flex_js_middle.sh {wildcards.sample} {output.js_flex_output} \
{input.fasta} {input.ref_GC_bedgraph} {input.ref_zipped_ORF_gff} \
{input.basepath_zipped_SNP_vcf} {input.basepath_sorted_bam}
        """


rule HTML_IGVJs_generate_final:
    input:
        expand( "{p}{chunk_name}_{sample}",
                p           =   f"{datadir + html}",
                chunk_name  =   [   '2_tab',
                                    '4_html_divs',
                                    '6_js_flex'
                                    ],
                sample      =   SAMPLES
                )
    output:
        f"{res}igv.html"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}HTML_IGVJs_generate_final.log"
    benchmark:
        f"{logdir + bench}benchmark/HTML_IGVJs_generate_final.txt"
    threads: 1
    params:
        tab_basename    =   f"{datadir + html}2_tab_",
        div_basename    =   f"{datadir + html}4_html_divs_",
        js_flex_output  =   f"{datadir + html}6_js_flex_"
    shell:
        """
cat files/html_chunks/1_header.html > {output}
cat {params.tab_basename}* >> {output}
cat files/html_chunks/3_tab_explanation_meta.html >> {output}
cat {params.div_basename}* >> {output}
cat files/html_chunks/5_js_begin.html >> {output}
cat {params.js_flex_output}* >> {output}
cat files/html_chunks/7_js_end.html >> {output}
        """

