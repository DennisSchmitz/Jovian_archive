#!/bin/bash

source bin/functions.sh

SET_HOSTNAME=$(bin/gethostname.sh)

minispacer
line
ARCHIVE_NAME=$(ls archive_*.tar.gz)
echo -e "Unpacking ${ARCHIVE_NAME}"
(tar -xzf archive_*.tar.gz) &
spinner

minispacer
line
minispacer

if [ -e profile/variables.yaml ]; then
    echo -e "Resetting some of the variables..."
    sed '3,4d' -i profile/variables.yaml
    echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml
fi

minispacer
line
minispacer

eval $(parse_yaml profile/variables.yaml "vars_")


if [ ! -h $CONDA_PREFIX/etc/nginx/default-site/$vars_Jovian_run_identifier ]; then
    echo -e "Setting symlink"
    bin/set_symlink.sh
else
    echo "Symlink has already been set"
fi

