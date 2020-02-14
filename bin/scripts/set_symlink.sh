#!/bin/bash
source bin/includes/functions

eval $(parse_yaml config/variables.yaml "vars_")

currentuser=$(whoami)

find /tmp/etc/nginx/default-site/* -user ${currentuser} -xtype l -exec unlink {} \;

if [ ! -h /tmp/etc/nginx/default-site/$vars_Jovian_run_identifier ]
   then
      mkdir -p /tmp/etc/nginx/default-site/
      ln -s "$(pwd)" "/tmp/etc/nginx/default-site/$vars_Jovian_run_identifier"
fi