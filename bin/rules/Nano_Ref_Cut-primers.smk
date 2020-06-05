


rule Cut_primers:
    input:
        fastq   =   rules.Adapter_trimming.output.trimmeddata,
        primers =   primerfile
    output:
        cleaneddata_pt1 =   f"{datadir + cln + prdir}" + "{sample}.fastq"
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "Primer_removal_{sample}.log"
    benchmark:
        f"{logdir + bench}" + "Primer_removal_{sample}.txt"
    threads: config["threads"]["Nanopore_primer_removal"]
    params:
        primer_cutoff_plus  =   config["Nanopore_ref"]["Primer_cutoff_plus"],
        primer_cutoff_minus =   config["Nanopore_ref"]["Primer_cutoff_minus"]
    shell:
        """
cutadapt \
--cores={threads} \
--cut {params.primer_cutoff_plus} \
--cut {params.primer_cutoff_minus} \
--revcomp -b file:{input.primers} \
-o {output.cleaneddata_pt1} {input.fastq} > {log} 2>&1
        """