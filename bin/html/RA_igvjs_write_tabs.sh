#!/bin/bash

#####
# This script (part 2) writes the various tabs for the IGVJS html.
# This script should be called for every individual sample

INPUT="$1"
OUTPUT="$2"
# $3 is input file from smk rule, it's only to trick the smk into making a DAG it's not used in this script
OUTPUT_HTML="$4"

SAMPLE="sample_${INPUT//-/_}"

cat << EOF >> ${OUTPUT_HTML}
	<li><a href="#${SAMPLE}">${SAMPLE}</a></li>
EOF
touch "${OUTPUT}"