# Reformat blast tsv output to gff	
## Remove all vector/construct/synthetic entries (because the filtering based on their specific taxid is not adequate)	
## Remove any entry with a lower bitscore than the user specified bitscore_threshold (i.e. filter short alignments since every match is a +2 bitscore)
rule make_gff: 
    input:
        rules.Scaffold_classification.output
    output:
        f"{datadir + taxclas}" + "{sample}_lca_raw.gff" #? This is a temp file, removed in the onSuccess//onError clause.
    conda:
        f"{conda_envs}mgkit_lca.yaml"
    log:
        f"{logdir}" + "make_gff_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "make_gff_{sample}.txt"
    threads: 1
    resources:
        memory = 12
    params:
        bitscore_threshold  =   config["Illumina_meta"]["LCA"]["bitscore_threshold"],
        filt_keywords       =   "/vector\|construct\|synthetic/Id" 
    shell:
        """
sed -i "{params.filt_keywords}" {input}; 
blast2gff blastdb -b {params.bitscore_threshold} -n {input} {output} > {log} 2>&1
        """

# Reformat gff with accession id blasthit into taxid
rule addtaxa_gff:
    input:
        rules.make_gff.output
    output:
        f"{datadir + taxclas}" + "{sample}_lca_tax.gff" #? This is a temp file, removed in the onSuccess//onError clause.
    conda:
        f"{conda_envs}mgkit_lca.yaml"
    log:
        f"{logdir}" + "addtaxa_gff_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "addtaxa_gff_{sample}.txt"
    threads: 1
    resources:
        memory = 12
    params:
        mgkit_tax_db    =   config["databases"]["MGKit_taxonomy"]
    shell:
        """
add-gff-info addtaxa -t <(gunzip -c {params.mgkit_tax_db}nucl_gb.accession2taxid.gz | cut -f2,3) -e {input} {output} > {log} 2>&1  
        """

# Filter taxid 81077 (https://www.ncbi.nlm.nih.gov/taxonomy/?term=81077 --> artificial sequences) and 12908 (https://www.ncbi.nlm.nih.gov/taxonomy/?term=12908 --> unclassified sequences)
rule taxfilter_gff:
    input:
        rules.addtaxa_gff.output
    output:
        f"{datadir + taxclas}" + "{sample}_lca_taxfilt.gff" #? This is a temp file, removed in the onSuccess//onError clause.
    conda:
        f"{conda_envs}mgkit_lca.yaml"
    log:
        f"{logdir}" + "taxfilter_gff_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "taxfilter_gff_{sample}.txt"
    threads: 1
    resources:
        memory = 12
    params:
        mgkit_tax_db    =   config["databases"]["MGKit_taxonomy"]
    shell:
        """
taxon-utils filter -e 81077 -e 12908 -t {params.mgkit_tax_db}taxonomy.pickle {input} {output} > {log} 2>&1
        """

# Filter gff on the user-specified bitscore-quantile settings.
rule qfilter_gff:
    input:
        rules.taxfilter_gff.output
    output:
        f"{datadir + taxclas}" + "{sample}_lca_filt.gff" #? This is a temp file, removed in the onSuccess//onError clause.
    conda:
        f"{conda_envs}mgkit_lca.yaml"
    log:
        f"{logdir}" + "qfilter_gff_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "qfilter_gff_{sample}.txt"
    threads: 1
    resources:
        memory = 12
    params:
        quantile_threshold  =   config["Illumina_meta"]["LCA"]["quantile_threshold"]
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
        filtgff =   rules.qfilter_gff.output,
        stats   =   rules.Generate_contigs_metrics.output.perScaffold
    output:
        no_lca      =   f"{datadir + taxclas}" + "{sample}_nolca_filt.gff", #? This is a temp file, removed in the onSuccess//onError clause.
        taxtab      =   f"{datadir + taxclas}" + "{sample}.taxtab",
        taxMagtab   =   f"{datadir + taxclas}" + "{sample}.taxMagtab"
    conda:
        f"{conda_envs}mgkit_lca.yaml"
    params:
        mgkit_tax_db        =   config["databases"]["MGKit_taxonomy"],
        bitscore_threshold  =   config["Illumina_meta"]["LCA"]["bitscore_threshold"]   
    log:
        f"{logdir}" + "lca_mgkit_{sample}.log" 
    benchmark:
        f"{logdir + bench}" + "lca_mgkit_{sample}.txt"
    threads: 1
    resources:
        memory = 12
    shell:
        """
taxon-utils lca -b {params.bitscore_threshold} -s -p -n {output.no_lca} -t {params.mgkit_tax_db}taxonomy.pickle {input.filtgff} {output.taxtab} > {log} 2>&1;
sed -i '1i #queryID\ttaxID' {output.taxtab} >> {log} 2>&1;
if [[ ! -e {output.no_lca} ]]; then
touch {output.no_lca}
fi
python bin/scripts/average_logevalue_no_lca.py {output.taxtab} {output.no_lca} {input.filtgff} {output.taxtab} >> {log} 2>&1;
python bin/scripts/krona_magnitudes.py {output.taxtab} {input.stats} {output.taxMagtab} >> {log} 2>&1
        """