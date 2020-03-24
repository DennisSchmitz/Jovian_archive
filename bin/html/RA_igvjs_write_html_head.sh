#!/bin/bash

#####
# This script (part 1) writes the standard header information
# it should only be called once

OUTPUT_HTML="$1"

cat << EOF > ${OUTPUT_HTML}
<!doctype html>
<head>

    <!-- Adapted from https://www.elated.com/res/File/articles/development/javascript/document-object-model/javascript-tabs/javascript-tabs.html -->

    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no">
    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="shortcut icon" href=../img/favicon.ico>
    <title>igv.js</title>

    <!-- IGV JS-->
    <script src="https://cdn.jsdelivr.net/npm/igv@2.3.5/dist/igv.min.js"></script>

    <style type="text/css">
        body {
            font-size: 80%;
            font-family: 'Lucida Grande', Verdana, Arial, Sans-Serif;
        }

        ul#tabs {
            list-style-type: none;
            margin: 30px 0 0 0;
            padding: 0 0 0.3em 0;
        }

        ul#tabs li {
            display: inline-block;
            margin-bottom: 15px;
        }

        ul#tabs li a {
            color: #42454a;
            background-color: #dedbde;
            padding: 0.3em;
            text-decoration: none;
            border: 1px solid #333333;
        }

        ul#tabs li a:hover {
            background-color: #f1f0ee;
        }

        ul#tabs li a.selected {
            color: #fff;
            background-color: #333;
            font-weight: bold;
            padding: 0.3em 0.3em 0.3em 0.3em;
        }

        div.tabContent {
            border: 1px solid #c9c3ba;
            padding: 0.5em;
        }

        div.tabContent.hide {
            display: none;
        }
    </style>

</head>

<body>

<h1>Alignments from a BAM file</h1>

<ul id="tabs">
    <li><a href="#standardTab">Explanation</a></li>
EOF

