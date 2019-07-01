rule make_gff: 
    input:
         "data/taxonomic_classification/{sample}.blastn"
    output:
         temp("data/taxonomic_classification/{sample}_lca_raw.gff")
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/make_gff_{sample}.txt"
    threads: 1
    log:
        "logs/make_gff_{sample}.log"
    params:
        bitscore_threshold=config["taxonomic_classification_LCA"]["mgkit_LCA"]["bitscore_threshold"] 
    shell:
        """
        sed -i "/vector\|construct\|synthetic/Id" {input}; 
        blast2gff blastdb -b {params.bitscore_threshold} -n {input} {output} > {log} 2>&1
        """
rule addtaxa_gff:
    input:
         "data/taxonomic_classification/{sample}_lca_raw.gff"
    output:
         temp("data/taxonomic_classification/{sample}_lca_tax.gff")
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/addtaxa_gff_{sample}.txt"
    threads: 1
    log:
        "logs/addtaxa_gff_{sample}.log"
    shell:
        """
         add-gff-info addtaxa -t <(gunzip -c /mnt/db/taxdb/nucl_gb.accession2taxid.gz | cut -f2,3) -e {input} {output} > {log} 2>&1  
        """
        
rule taxfilter_gff:
    input:
         "data/taxonomic_classification/{sample}_lca_tax.gff"
    output:
         temp("data/taxonomic_classification/{sample}_lca_taxfilt.gff")
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/taxfilter_gff_{sample}.txt"
    threads: 1
    log:
        "logs/taxfilter_gff_{sample}.log"
    shell:
        """
         taxon-utils filter -e 81077 -e 12908 -t /mnt/db/taxdb/taxonomy.pickle {input} {output} > {log} 2>&1
        """

rule qfilter_gff:
    input:
         "data/taxonomic_classification/{sample}_lca_taxfilt.gff"
    output:
         temp("data/taxonomic_classification/{sample}_lca_filt.gff")
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/qfilter_gff_{sample}.txt"
    threads: 1
    log:
        "logs/qfilter_gff_{sample}.log"
    params:
        quantile_threshold=config["taxonomic_classification_LCA"]["mgkit_LCA"]["quantile_threshold"],
    shell:
        """
         filter-gff sequence -t -f quantile -l {params.quantile_threshold} -c ge {input} {output} > {log} 2>&1
        """       
rule lca_mgkit:
    input:
        filtgff="data/taxonomic_classification/{sample}_lca_filt.gff",
        blast="data/taxonomic_classification/{sample}.blastn",
        stats="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats"
    output:
        no_lca=temp("data/taxonomic_classification/{sample}_nolca_filt.gff"),
        taxtab="data/taxonomic_classification/{sample}.taxtab",
        taxMagtab="data/taxonomic_classification/{sample}.taxMagtab"
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/lca_mgkit_{sample}.txt"
    threads: 1
    log:
        "logs/lca_mgkit_{sample}.log"
    params:
        bitscore_threshold=config["taxonomic_classification_LCA"]["mgkit_LCA"]["bitscore_threshold"]        
    shell:
        """
        taxon-utils lca -b {params.bitscore_threshold} -s -p -n {output.no_lca} -t /mnt/db/taxdb/taxonomy.pickle {input.filtgff} {output.taxtab} > {log} 2>&1;
        sed -i '1i #queryID\ttaxID' {output.taxtab} >> {log} 2>&1;
        if [[ ! -e {output.no_lca} ]]; then
        touch {output.no_lca}
        fi
        python bin/average_logevalue_no_lca.py {output.taxtab} {output.no_lca} {input.filtgff} {output.taxtab} >> {log} 2>&1;
        python bin/krona_magnitudes.py {output.taxtab} {input.blast} {input.stats} {output.taxMagtab} >> {log} 2>&1
        """