#!/bin/bash

#####
# This script (part 3) closes the html tabs-list and writes the basic explanation
# this script should be called only once

OUTPUT="$1"
OUTPUT_HTML="$2"

cat << EOF >> ${OUTPUT_HTML}
</ul>

<div id="standardTab">
    <ul>
        <li> Pick a sample to open the Interactive Genome Viewer </li>
        <li> The first "Node" within a sample is automatically loaded, you can view a specific node (scaffold) with the dropdown menu. The scaffolds are ordered from large to small. </li>
        <li> Be aware that loading might take a very long time depending on the size of the scaffold. This usually is a non-issue for scaffolds smaller than 28k nucleotides.
    </ul>
</div>
EOF

touch "${OUTPUT}"