#!/bin/bash

#####
# This script (part 5) this script writes the first (static) part of the required JavaScript
# this script should be called only once.

OUTPUT="$1"
OUTPUT_HTML="$2"

cat << EOF >> ${OUTPUT_HTML}
<script type="text/javascript">

    document.addEventListener("DOMContentLoaded", function () {
        initTabs();
        initIGV();
    });

    function initIGV() {

        var igvDiv, options;
EOF

touch "${OUTPUT}"