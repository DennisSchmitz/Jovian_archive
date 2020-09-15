#!/bin/bash

#####
# This script (part 3) this script writes the flexible part of the required JavaScript
# this script should be called for every sample.

INPUT="$1"
OUTPUT_HTML="$2"
INPUT_REF_FASTA="$3"
INPUT_REF_GC_BEDGRAPH="$4"
INPUT_REF_ZIPPED_ORF_GFF="$5"
INPUT_MAJOR_SNP_VCF_GZ="$6"
INPUT_MINOR_SNP_VCF_GZ="$7"
INPUT_SORTED_BAM="$8"


SAMPLE="sample_${INPUT//-/_}"


parse_yaml() {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/7;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

eval $(parse_yaml config/variables.yaml "vars_")
eval $(parse_yaml config/pipeline_parameters.yaml "params_")

cat << EOF >> ${OUTPUT_HTML}
        ${SAMPLE} = document.getElementById("${SAMPLE}");
        options =
            {
                reference:
                {
                    id: "${INPUT}",
                    fastaURL: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_REF_FASTA}",
                    wholeGenomeView: false,
                    tracks: [
                        {
                            type: "wig",
                            name: "GC contents",
                            format: "bedGraph",
                            color: "#fdae61",
                            url: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_REF_GC_BEDGRAPH}",
                            min: "0",
                            max: "1",
                            autoHeight: "true",
                            minHeight: "75",
                            maxHeight: "125",
                            order: Number.MAX_VALUE
                        },
                        {
                            name: "Majority SNPs",
                            type: "variant",
                            format: "vcf",
                            color: "#4dac26",
                            url: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_MAJOR_SNP_VCF_GZ}",
                            indexURL: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_MAJOR_SNP_VCF_GZ}.tbi",
                            displayMode: "EXPANDED",
                            order: 1
                        },
                        {
                            name: "Minority SNPs",
                            type: "variant",
                            format: "vcf",
                            color: "#d01c8b",
                            url: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_MINOR_SNP_VCF_GZ}",
                            indexURL: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_MINOR_SNP_VCF_GZ}.tbi",
                            displayMode: "EXPANDED",
                            order: 3
                        },
                        {
                            type: "alignment",
                            format: "bam",
                            colorBy: "strand",
                            url: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_SORTED_BAM}",
                            indexURL: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_SORTED_BAM}.bai",
                            indexed: "true",
                            name: "Alignment",
                            showSoftClips: false,
                            viewAsPairs: true,
                            order: 4
                        },
                        {
                            type: "annotation",
                            name: "ORF predictions",
                            format: "gff3",
                            color: "#2c7bb6",
                            url: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_REF_ZIPPED_ORF_GFF}",
                            indexURL: "${vars_Server_host_hostname}:${params_Global_server_port}/${vars_Jovian_run_identifier}/${INPUT_REF_ZIPPED_ORF_GFF}.tbi",
                            displayMode: "EXPANDED",
                            order: 2
                        }
                    ]
                },
            };
        igv.createBrowser(${SAMPLE}, options);

EOF
