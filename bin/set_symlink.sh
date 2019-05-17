#!/bin/bash
source bin/functions.sh

eval $(parse_yaml profile/variables.yaml "vars_")


if [ ! -h /tmp/etc/nginx/default-site/$vars_Jovian_run_identifier ]
   then
      mkdir -p /tmp/etc/nginx/default-site/
      ln -s "$(pwd)/bin/software/igv.js" "/tmp/etc/nginx/default-site/$vars_Jovian_run_identifier"
fi

chmod 777 /tmp/etc/nginx/default-site/*