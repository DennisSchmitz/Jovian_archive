#!/bin/bash

#####
# This script (part 3) closes the html tabs-list and writes the basic explanation
# this script should be called only once

OUTPUT="$1"

cat << EOF >> results/igv.html
</ul>

<div id="standardTab">
    <ul>
        <li> When this tab (really, just a div) comes first in the document, the igv.js div has display:none, width == 0
        <li> Clicking the IGV Tab does not generate any events that igv can listen to. So to inform IGV of the
            visibility
            change clients must call igv.visibilityChange(), as shown above in the showTabs function.
    </ul>
</div>
EOF

touch "${OUTPUT}"