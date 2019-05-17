#!/bin/bash
# make CORS enabled nginx config (https://enable-cors.org/server_nginx.html)
# and start nginx in background
# parameter start or stop
#
# changes
# 22-11-2018, rv, from 8080 to port 8083 because jupyter notebook is using this port

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

eval $(parse_yaml profile/pipeline_parameters.yaml "config_")

echo "Selected port for NGINX is: $config_server_info_port "
port=$config_server_info_port

#set -x

# check parameter
if [ -z "${1}" ]
then
   echo "start or stop?"
   exit 1
fi

#CONDA_PREFIX=$2

# check if nginx is installed
if [ -z "${CONDA_PREFIX}" ] || [ ! -d "${CONDA_PREFIX}" ]
then
   echo "CONDA_PREFIX env not set"
   exit 1
fi

mkdir -p /tmp/etc/nginx/sites.d/
cp $CONDA_PREFIX/etc/nginx/sites.d/default-site.conf /tmp/etc/nginx/sites.d/default-site.conf

if [ ! -f "/tmp/etc/nginx/sites.d/default-site.conf" ] || [ ! -f "$CONDA_PREFIX/bin/nginx" ]
then
   echo "nginx/cfg not found"
   echo "conda install nginx?"
   exit 1
fi

# make config and start nginx
if [ "${1}" == "start" ]
then
   # check if nginx already runs
   if [ ! -z "$(ps -ef| grep nginx| grep worker)" ]
   then
       echo "nginx is already running"
       exit 1
   fi

   cat << EOF > /tmp/etc/nginx/sites.d/default-site.conf
server {
    listen       ${port};
    server_name  localhost;
    root   /tmp/etc/nginx/default-site/;
    location / {
        index  index.html index.htm;

        if (\$request_method = 'OPTIONS') {
           add_header 'Access-Control-Allow-Origin' '*';
           add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS';
           #
           # Custom headers and headers various browsers *should* be OK with but aren't
           #
           add_header 'Access-Control-Allow-Headers' 'DNT,User-Agent,X-Requested-With,If-Modified-Since,Cache-Control,Content-Type,Range';
           #
           # Tell client that this pre-flight info is valid for 20 days
           #
           add_header 'Access-Control-Max-Age' 1728000;
           add_header 'Content-Type' 'text/plain; charset=utf-8';
           add_header 'Content-Length' 0;
           return 204;
        }

        if (\$request_method = 'POST') {
           add_header 'Access-Control-Allow-Origin' '*';
           add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS';
           add_header 'Access-Control-Allow-Headers' 'DNT,User-Agent,X-Requested-With,If-Modified-Since,Cache-Control,Content-Type,Range';
           add_header 'Access-Control-Expose-Headers' 'Content-Length,Content-Range';
        }

        if (\$request_method = 'GET') {
           add_header 'Access-Control-Allow-Origin' '*';
           add_header 'Access-Control-Allow-Methods' 'GET, POST, OPTIONS';
           add_header 'Access-Control-Allow-Headers' 'DNT,User-Agent,X-Requested-With,If-Modified-Since,Cache-Control,Content-Type,Range';
           add_header 'Access-Control-Expose-Headers' 'Content-Length,Content-Range';
        }
    }

    #error_page  404              /404.html;

    # redirect server error pages to the static page /50x.html
    error_page   500 502 503 504  /50x.html;
    location = /50x.html {
        root   /tmp/etc/nginx/default-site/;
    }

}
EOF

cp /tmp/etc/nginx/sites.d/default-site.conf $CONDA_PREFIX/etc/nginx/sites.d/default-site.conf
chmod 777 /tmp/etc/nginx/sites.d/
chmod 777 /tmp/etc/nginx/sites.d/*
chmod 777 /tmp/etc/nginx/default-site/
chmod 777 /tmp/etc/nginx/default-site/*
   nginx&
   sleep 3

   # check if nginx runs
   if [ -z "$(ps -ef| grep nginx| grep worker)" ]
   then
       echo "nginx is not running"
       exit 1
   fi
fi


if [ "${1}" == "stop" ]
then
   # check if nginx  runs
   if [ -z "$(ps -ef| grep nginx| grep worker)" ]
   then
       echo "nginx is not running"
       exit 1
   else
      nginx -s quit
   fi
fi
