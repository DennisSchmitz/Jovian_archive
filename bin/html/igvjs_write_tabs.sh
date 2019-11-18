#!/bin/bash

#####
# This script (part 2) writes the various tabs for the IGVJS html.
# This script should be called for every individual sample

INPUT="$1"
OUTPUT="$2"

SAMPLE="sample_${INPUT//-/_}"

cat << EOF >> results/igv.html
	<li><a href="#${SAMPLE}">${SAMPLE}</a></li>
EOF
touch "${OUTPUT}"