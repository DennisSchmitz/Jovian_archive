#!/bin/bash

#####
# This script (part 2) writes the various tabs for the IGVJS html.
# This script should be called for every individual sample

INPUT="$1"
OUTPUT_HTML="$2"

SAMPLE="sample_${INPUT//-/_}"

cat << EOF >> ${OUTPUT_HTML}
	<li><a href="#${SAMPLE}">${SAMPLE}</a></li>
EOF
