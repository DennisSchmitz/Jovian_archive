#!/bin/bash

#####
# This script (part 4) this script writes the various div elements for the samples
# this script should be called for every sample individually.

INPUT="$1"
OUTPUT="$2"

SAMPLE=${INPUT//_/-}

cat << EOF >> results/igv.html
<div id="${SAMPLE}">
    <div id="${SAMPLE}"></div>
</div>
EOF
touch "${OUTPUT}"