
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#@ 
#@ 
#@ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

rule QC_raw_data:
    input:
        lambda wildcards: SAMPLES[wildcards.sample][wildcards.read]
    output:
        html    =   f"{datadir + qc_pre}" + "{sample}_{read}_fastqc.html",
        zip     =   f"{datadir + qc_pre}" + "{sample}_{read}_fastqc.zip"
    conda:
        f"{conda_envs}QC_and_clean.yaml"
    log:
        f"{logdir}" + "QC_raw_data_{sample}_{read}.log"
    benchmark:
        f"{logdir + bench}" + "QC_raw_data_{sample}_{read}.txt"
    threads: 6
    params:
        output_dir  =   f"{datadir + qc_pre}"
    shell:
        """
bash bin/scripts/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log}
        """