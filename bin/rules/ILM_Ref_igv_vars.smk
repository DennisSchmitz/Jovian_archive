rule Illumina_HTML_IGVJs_variable_parts:
    input:
        fasta                   =   rules.Illumina_index_reference.output.reference_copy,
        ref_GC_bedgraph         =   rules.Illumina_determine_GC_content.output.GC_bed,
        ref_zipped_ORF_gff      =   rules.Illumina_reference_ORF_analysis.output.zipped_gff3,
        basepath_zipped_SNP_vcf =   rules.Illumina_extract_raw_consensus.output.majorSNP_vcf_gz,
        basepath_sorted_bam     =   rules.Illumina_align_to_reference.output.sorted_bam
    output:
        tab_output      =   f"{datadir + html}" + "2_tab_{sample}",
        div_output      =   f"{datadir + html}" + "4_html_divs_{sample}",
        js_flex_output  =   f"{datadir + html}" + "6_js_flex_{sample}"
    conda:
        f"{conda_envs}data_wrangling.yaml"
    log:
        f"{logdir}" + "Illumina_HTML_IGVJs_variable_parts_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Illumina_HTML_IGVJs_variable_parts_{sample}.txt"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_tabs.sh {wildcards.sample} {output.tab_output}

bash bin/html/igvjs_write_divs.sh {wildcards.sample} {output.div_output}

bash bin/html/Illumina_igvjs_write_flex_js_middle.sh {wildcards.sample} {output.js_flex_output} \
{input.fasta} {input.ref_GC_bedgraph} {input.ref_zipped_ORF_gff} \
{input.basepath_zipped_SNP_vcf} {input.basepath_sorted_bam}
        """