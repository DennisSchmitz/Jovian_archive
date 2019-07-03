#!/bin/bash

### functions

line() {
    printf "%-70s\n" "#" | sed 's/\s/#/g'
}

thinline() {
    printf "%-70s\n" "-" | sed 's/\s/-/g'
}

spacer() {
    printf "\n\n"
}

minispacer() {
    printf "\n"
}

spinner() {
    PID=$!
    i=1
    sp="/-\|"
    echo -n ' '
    while [ -d /proc/$PID ]
    do
        printf "\b${sp:i++%${#sp}:1}"
        sleep 0.5
    done
    printf "\b"
}

installer_intro() {
    tput reset
    line
    echo -e "Welcome to the interactive installation process of Jovian"
    echo -e "You are seeing this because one or multiple dependencies of Jovian are missing on your system."
    echo -e "You are using Jovian $VERSION"
    spacer
    echo -e "Please keep in mind that this pipeline is still a work-in-progress"
    echo -e "It is expected that several portions are unstable until release 1.0.0"
    line
    spacer
}

database_installer_intro () {
    tput reset
    line
    echo -e "Welcome to the interactive database installation process of Jovian"
    echo -e "This process will ask you for the locations where to place the specific databases and will download/install them for you"
    minispacer
    echo -e "\e[1mWe strongly advise you to contact your local system administrators before continuing. \nIf you're inexperienced with the commandline, or linux computers in general, it is best to leave this installation process to your system administrators.\e[0m"
    minispacer
    line
    minispacer
}

database_installer () {
    tput reset
    line
    echo -e "Welcome to the interactive database installation process of Jovian"
    echo -e "This process will ask you for the locations where to place the specific databases and will download/install them for you"
    minispacer
    line
    minispacer
}

ready_for_start() {
    spacer
    line
    echo -e "All pre-flight checks of Jovian have been completed"
    echo -e "You are using Jovian version ${VERSION}"
    echo -e "Please keep in mind that this pipeline is still a work-in-progress"
    echo -e "It is expected that several portions are unstable until release 1.0.0"
    line
    spacer
}

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
         printf("%s%s%s=\"%s\"\n", "'${prefix}'",vn, $2, $3);
      }
   }'
}
