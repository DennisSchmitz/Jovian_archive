#!/bin/bash

#####
# This script (part 4) this script writes the various div elements for the samples
# this script should be called for every sample individually.

INPUT="$1"
OUTPUT_HTML="$2"

SAMPLE="sample_${INPUT//-/_}"

cat << EOF >> ${OUTPUT_HTML}
<div id="${SAMPLE}">
    <div id="${SAMPLE}"></div>
</div>
EOF