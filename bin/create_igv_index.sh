#!/bin/bash

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

tree -H "${vars_Server_host_hostname}:${params_server_info_port}/${vars_Jovian_run_identifier}/Jovian/data/scaffolds_filtered" -L 1 -T "${params_HTML_index_titles_IGVjs_title}" --noreport --charset utf-8 -P "*.html" -o results/IGVjs_index.html data/scaffolds_filtered/
#mv IGVjs_index.html results/IGVjs_index.html