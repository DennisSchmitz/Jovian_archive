#!/bin/bash
source bin/functions.sh

eval $(parse_yaml profile/variables.yaml "vars_")


if [ ! -h $CONDA_PREFIX/etc/nginx/default-site/$vars_Jovian_run_identifier ]
   then
      ln -s "$(pwd)/bin/software/igv.js" "$CONDA_PREFIX/etc/nginx/default-site/$vars_Jovian_run_identifier"
fi