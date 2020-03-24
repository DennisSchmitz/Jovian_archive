#!/bin/bash

#####
# This script (part 3) closes the html tabs-list and writes the basic explanation
# this script should be called only once

OUTPUT_HTML="$1"

cat << EOF >> ${OUTPUT_HTML}
</ul>

<div id="standardTab">
    <ul>
        <li> Pick a sample to open the Interactive Genome Viewer </li>
        <li> The first "segment" within the given reference is automatically loaded, you can view a specific segment with the dropdown menu. </li>
        <li> Be aware that loading might take a very long time depending on the size of the scaffold and the sequencing depth. 
    </ul>
</div>
EOF
