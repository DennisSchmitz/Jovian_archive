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

eval $(parse_yaml profile/variables.yaml "vars_")
oldhostname="$vars_Server_host_hostname"
newhostname="http://${SET_HOSTNAME}"


###### Modifying the HTML index file for IGVjs
if [ -e results/IGVjs_index.html ]; then
    echo -e "Modifying IGVjs index file..."
    sed -i -e "s@${oldhostname}@${newhostname}@g" "results/IGVjs_index.html"
fi


if [ -e profile/variables.yaml ]; then
    echo -e "Resetting some of the variables..."
    sed '3,4d' -i profile/variables.yaml
    echo -e "Server_host:\n    hostname: http://${SET_HOSTNAME}" >> profile/variables.yaml
fi

minispacer
line
minispacer



if [ ! -h $CONDA_PREFIX/etc/nginx/default-site/$vars_Jovian_run_identifier ]; then
    echo -e "Setting symlink..."
    bin/set_symlink.sh
else
    echo "Symlink has already been set"
fi

echo -e "Done"
echo -e ":)"
