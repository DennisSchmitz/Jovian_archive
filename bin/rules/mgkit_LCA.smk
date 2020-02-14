# Reformat blast tsv output to gff	
## Remove all vector/construct/synthetic entries (because the filtering based on their specific taxid is not adequate)	
## Remove any entry with a lower bitscore than the user specified bitscore_threshold (i.e. filter short alignments since every match is a +2 bitscore)
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

# Reformat gff with accession id blasthit into taxid
rule addtaxa_gff:
    input:
         "data/taxonomic_classification/{sample}_lca_raw.gff"
    output:
         temp("data/taxonomic_classification/{sample}_lca_tax.gff")
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/addtaxa_gff_{sample}.txt"
    params:
        mgkit_tax_db=config["databases"]["MGKit_taxonomy"]
    threads: 1
    log:
        "logs/addtaxa_gff_{sample}.log"
    shell:
        """
         add-gff-info addtaxa -t <(gunzip -c {params.mgkit_tax_db}nucl_gb.accession2taxid.gz | cut -f2,3) -e {input} {output} > {log} 2>&1  
        """

# Filter taxid 81077 (https://www.ncbi.nlm.nih.gov/taxonomy/?term=81077 --> artificial sequences) and 12908 (https://www.ncbi.nlm.nih.gov/taxonomy/?term=12908 --> unclassified sequences)
rule taxfilter_gff:
    input:
         "data/taxonomic_classification/{sample}_lca_tax.gff"
    output:
         temp("data/taxonomic_classification/{sample}_lca_taxfilt.gff")
    conda:
        "../envs/mgkit_lca.yaml"
    benchmark:
        "logs/benchmark/taxfilter_gff_{sample}.txt"
    params:
        mgkit_tax_db=config["databases"]["MGKit_taxonomy"]
    threads: 1
    log:
        "logs/taxfilter_gff_{sample}.log"
    shell:
        """
         taxon-utils filter -e 81077 -e 12908 -t {params.mgkit_tax_db}taxonomy.pickle {input} {output} > {log} 2>&1
        """

# Filter gff on the user-specified bitscore-quantile settings.
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
        
# Perform the LCA analysis.	
## in `taxon-utils lca` the `-n {output.no_lca}` flag is important because these are the entries that reach no LCA and therefore have no taxid and are required later for the `bin/average_logevalue_no_lca.py` script.	
## in `taxon-utils lca` the `-b {params.bitscore_threshold} ` is redundant because this is already done in the first rule.	
## the `sed` rule creates a header for the output file	
## touch the {output.no_lca} if it hasn't been generated yet, otherwise you get an error downstream	
## `bin/average_logevalue_no_lca.py` add a `taxid=1` and `evalue=1` for all entries without an LCA result and also average the e-values (if e-value is 0 it is set to 10log evalue of -450).	
## `bin/krona_magnitudes.py` adds magnitude information for the Krona plot (same as default Krona method).
rule lca_mgkit:
    input:
        filtgff="data/taxonomic_classification/{sample}_lca_filt.gff",
        stats="data/scaffolds_filtered/{sample}_perMinLenFiltScaffold.stats"
    output:
        no_lca=temp("data/taxonomic_classification/{sample}_nolca_filt.gff"),
        taxtab="data/taxonomic_classification/{sample}.taxtab",
        taxMagtab="data/taxonomic_classification/{sample}.taxMagtab"
    conda:
        "../envs/mgkit_lca.yaml"
    params:
        mgkit_tax_db=config["databases"]["MGKit_taxonomy"],
        bitscore_threshold=config["taxonomic_classification_LCA"]["mgkit_LCA"]["bitscore_threshold"]   
    benchmark:
        "logs/benchmark/lca_mgkit_{sample}.txt"
    threads: 1
    log:
        "logs/lca_mgkit_{sample}.log"         
    shell:
        """
        taxon-utils lca -b {params.bitscore_threshold} -s -p -n {output.no_lca} -t {params.mgkit_tax_db}taxonomy.pickle {input.filtgff} {output.taxtab} > {log} 2>&1;
        sed -i '1i #queryID\ttaxID' {output.taxtab} >> {log} 2>&1;
        if [[ ! -e {output.no_lca} ]]; then
        touch {output.no_lca}
        fi
        python bin/average_logevalue_no_lca.py {output.taxtab} {output.no_lca} {input.filtgff} {output.taxtab} >> {log} 2>&1;
        python bin/krona_magnitudes.py {output.taxtab} {input.stats} {output.taxMagtab} >> {log} 2>&1
        """