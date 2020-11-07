rule Cut_primers:
    input:
        fastq       =   rules.Remove_Adapters_pt2.output,
        primers_5   =   rules.Prepare_primers.output.primers_5,
        primers_3   =   rules.Prepare_primers.output.primers_3
    output: f"{datadir + cln + prim}" + "{sample}.fastq"
    conda:
        f"{conda_envs}Nano_clean.yaml"
    log:
        f"{logdir}" + "Primer_removal_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Primer_removal_{sample}.txt"
    threads: config["threads"]["Nanopore_primer_removal"]
    params:
        primer_cutoff_plus  =   config["Nanopore_ref"]["Primer_cutoff_plus"],
        primer_cutoff_minus =   config["Nanopore_ref"]["Primer_cutoff_minus"],
        overlap             =   config["Nanopore_ref"]["Primer_min_overlap"],
        error_rate          =   config["Nanopore_ref"]["Primer_error_rate"],
        repeat_search       =   config["Nanopore_ref"]["Primer_repeat_search"]
    shell:
        """
cutadapt \
--cores={threads} \
--cut {params.primer_cutoff_plus} \
--cut {params.primer_cutoff_minus} \
-g file:{input.primers_5} \
-a file:{input.primers_3} \
--revcomp \
-O {params.overlap} \
-e {params.error_rate} \
--no-indels \
--times {params.repeat_search} \
-o {output} \
{input.fastq} > {log}
        """