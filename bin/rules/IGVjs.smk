
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule HTML_IGVJs_part1_static_head:
    output:
        "data/html/html_head.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_html_head.sh {output}
        """

rule HTML_IGVJs_part2_tabs:
    input:
        "data/html/html_head.ok"
    output:
        "data/html/html_tabs.{sample}.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_tabs.sh {wildcards.sample} {output} {input}
        """

rule HTML_IGVJs_part3_close_tabs:
    input:
        expand("data/html/html_tabs.{sample}.ok", sample = SAMPLES)
    output:
        "data/html/tabs_closed.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_close_tabs.sh {output} {input}
        """

rule HTML_IGVJs_part4_divs:
    input:
        "data/html/tabs_closed.ok"
    output:
        "data/html/html_divs.{sample}.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_divs.sh {wildcards.sample} {output} {input}
        """

rule HTML_IGVJs_part5_begin_js:
    input:
        expand("data/html/html_divs.{sample}.ok", sample = SAMPLES)
    output:
        "data/html/js-begin.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_static_js_begin.sh {output} {input}
        """

rule HTML_IGVJs_part6_middle_js:
    input:
        "data/html/js-begin.ok"
    output:
        "data/html/js-flex.{sample}.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_flex_js_middle.sh {wildcards.sample} {output} {input}
        """

rule HTML_IGVJs_part7_end_js:
    input:
        expand("data/html/js-flex.{sample}.ok", sample = SAMPLES)
    output:
        "data/html/js-end.ok"
    conda:
        "../envs/data_wrangling.yaml"
    threads: 1
    shell:
        """
bash bin/html/igvjs_write_static_js_end.sh {output} {input}
        """