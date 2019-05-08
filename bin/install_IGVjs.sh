#!/bin/bash
# script downloads and build IGVjs
#
# https://github.com/igvteam/igv.js?files=1
#

# set -x

tools="git npm nginx node"

home="bin/"
igvdir=${home}/software/igv.js

function test_npm_internet {
# test internet for npm
echo -n "test npm ping, " 
npm ping > /dev/null 2>&1
echo ok
let err=$?
if [ $err -ne 0 ]
then
   echo "can not connect to internet, try to set npm proxy..."
   if [ ! -z "${http_proxy}" ]
   then
       npm config set proxy http://${http_proxy}
   fi
   if [ ! -z "${https_proxy}" ]
   then
       npm config set https-proxy http://${https_proxy}
   fi
fi

if [ ${err} -ne 0 ]
then
   echo -n "test npm ping again" 
   npm ping > /dev/null 2>&1
   if [ $? -ne 0 ] 
   then
       echo
       echo "npm can not connect to internet"
       echo "try to set proxy"
       echo "npm config set proxy http://proxy.company.com:8080"
       echo "npm config set https-proxy http://proxy.company.com:8080"
       echo
       echo "proxy env settings:"
       set| grep -i proxy
       exit 1
   else
   echo ok
   fi
fi

}

# test if tools exists
echo "check tools:"
for name in $(echo $tools)
do
    which $name > /dev/null 2>&1
    if [ $? -ne 0 ]
    then
        echo $name not installed
        exit 1
    else 
       echo "`which $name` OK"
    fi
done
echo

# download IGVjs
if [ ! -d ${home}/software/ ]
then
   mkdir -p ${home}/software/
fi


if [ ! -f ${igvdir}/examples/igvjs.html ] || [ ! -f ${igvdir}/css/igv.css ] || [ ! -f ${igvdir}/js/igv-utils.js ]
then
   cd ${home}/software/
   git clone https://github.com/igvteam/igv.js.git
else
   echo igv.js already downloaded
   cd ${home}/software/
fi
echo

cd igv.js

if [ -d node_modules ] && [ -d node_modules/grunt ] && [ -f node_modules/grunt/bin/grunt ]
then
    echo "node modules already installed"
else
   echo "Install all Node.js modules using npm (this can take a while...) output in npm_install_output and npm_install_output_err"
   test_npm_internet
   npm install > npm_install_output 2> npm_install_output_err
   if [ $? -ne 0 ]
   then
       echo "npm gave an error:"
       cat npm_install_output_err
       echo "still continuing..."
       echo
   fi
fi
echo

if [ -f dist/igv.min.js ]
then
    echo "dist/igv.min.js already exist"
else
   echo "build dist/igv.min.jst (npm run grunt) output in npm_run_grunt_output and npm_run_grunt_output_err"
   test_npm_internet
   npm run grunt > npm_run_grunt_output 2> npm_run_grunt_output_err
   if [ $? -ne 0 ]
   then
       echo "npm run grunt gave an error:"
       cat npm_run_grunt_output_err
       echo
   fi
fi
echo

# maken link to generated htmls
if [ ! -h Jovian ]
then
   ln -s ../../.. Jovian
fi
