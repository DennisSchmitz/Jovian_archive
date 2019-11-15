#!/bin/bash

#####
# This script (part 6) this script writes the flexible part of the required JavaScript
# this script should be called for every sample.

INPUT="$1"
OUTPUT="$2"

SAMPLE=${INPUT//_/-}


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

eval $(parse_yaml profile/variables.yaml "vars_")
eval $(parse_yaml profile/pipeline_parameters.yaml "params_")

cat << EOF >> results/igv.html
        ${SAMPLE} = document.getElementById("${SAMPLE}");
        options =
            {
                reference:
                {
                    id: "${SAMPLE}",
                    fastaURL: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_scaffolds_ge${params_scaffold_minLen_filter_minlen}nt.fasta",
                    wholeGenomeView: false,
                    tracks: [
                        {
                            type: "wig",
                            name: "GC contents",
                            format: "bedGraph",
                            url: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_GC.bedgraph",
                            min: "0",
                            max: "1",
                            order: Number.MAX_VALUE
                        },
                        {
                            name:"SNPs",
                            type:"variant",
                            format:"vcf",
                            url: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_filtered.vcf.gz",
                            indexURL: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_filtered.vcf.gz.tbi",
                            displayMode: "SQUISHED",
                            order: 2
                        },
                        {
                            type: "alignment",
                            format: "bam",
                            colorBy: "strand",
                            url: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_sorted.bam",
                            indexURL: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_sorted.bam.bai",
                            indexed: "true",
                            name: "Alignment",
                            showSoftClips: true,
                            order: 3
                        },
                        {
                            type: "annotation",
                            name: "ORF predictions",
                            format: "gff3",
                            url: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_annotation.gff.gz",
                            indexURL: "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/data/scaffolds_filtered/${SAMPLE}_annotation.gff.gz.tbi",
                            displayMode: "EXPANDED",
                            order: 1
                        }
                    ]
                },
            };
        igv.createBrowser(${SAMPLE}, options);

EOF

touch "${OUTPUT}"