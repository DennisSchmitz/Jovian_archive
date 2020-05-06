
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


rule HTML_IGVJs_variable_parts:
    input:
        fasta= "data/scaffolds_filtered/{sample}_scaffolds_ge%snt.fasta" % config["scaffold_minLen_filter"]["minlen"],
        ref_GC_bedgraph= "data/scaffolds_filtered/{sample}_GC.bedgraph",
        ref_zipped_ORF_gff= "data/scaffolds_filtered/{sample}_annotation.gff.gz",
        basepath_zipped_SNP_vcf= "data/scaffolds_filtered/{sample}_filtered.vcf.gz",
        basepath_sorted_bam= "data/scaffolds_filtered/{sample}_sorted.bam",
    output:
        tab_output= "data/html/2_tab_{sample}",
        div_output= "data/html/4_html_divs_{sample}",
        js_flex_output= "data/html/6_js_flex_{sample}",
    conda:
        "../envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    log:
        "logs/HTML_IGVJs_variable_parts_{sample}.log"
    params:
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
        expand("data/html/{chunk_name}_{sample}", chunk_name = [ '2_tab', '4_html_divs', '6_js_flex' ], sample = SAMPLES)
    output:
        "results/igv.html"
    conda:
        "../envs/data_wrangling.yaml"
    benchmark:
        "logs/benchmark/HTML_IGVJs_generate_final.txt"
    threads: 1
    log:
        "logs/HTML_IGVJs_generate_final.log"
    params:
        tab_basename= "data/html/2_tab_",
        div_basename= "data/html/4_html_divs_",
        js_flex_output= "data/html/6_js_flex_",
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

