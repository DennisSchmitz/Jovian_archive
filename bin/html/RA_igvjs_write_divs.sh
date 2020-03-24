#!/bin/bash

#####
# This script (part 4) this script writes the various div elements for the samples
# this script should be called for every sample individually.

INPUT="$1"
OUTPUT="$2"
OUTPUT_HTML="$3"

SAMPLE="sample_${INPUT//-/_}"

cat << EOF >> ${OUTPUT_HTML}
<div id="${SAMPLE}">
    <div id="${SAMPLE}"></div>
</div>
EOF
touch "${OUTPUT}"