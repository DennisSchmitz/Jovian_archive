#!/bin/bash

#####################################################################################################################
### To do, ask Robert            ###
###     Example:                                                   ###
#####################################################################################################################

#set -x
#
#   id="7_20688_AGGCAGAAGCGTAAGA_L001_contigs"
#   fastaURL="../robert/7_20688_AGGCAGAAGCGTAAGA_L001_scaffolds_ge500nt.fasta"
#   url="../robert/7_20688_AGGCAGAAGCGTAAGA_L001_sorted.bam"
#   indexURL="../robert/7_20688_AGGCAGAAGCGTAAGA_L001_sorted.bam.bai"
#   indexed="true"
#   name="7_20688_AGGCAGAAGCGTAAGA_L001_sorted"

parse_yaml() {
   local prefix=$2
   local s='[[:space:]]*' w='[a-zA-Z0-9_]*' fs=$(echo @|tr @ '\034')
   sed -ne "s|^\($s\)\($w\)$s:$s\"\(.*\)\"$s\$|\1$fs\2$fs\3|p" \
        -e "s|^\($s\)\($w\)$s:$s\(.*\)$s\$|\1$fs\2$fs\3|p"  $1 |
   awk -F$fs '{
      indent = length($1)/7;
      vname[indent] = $2;
      for (i in vname) {if (i > indent) {delete vname[i]}}
      if (length($3) > 0) {
         vn=""; for (i=0; i<indent; i++) {vn=(vn)(vname[i])("_")}
         printf("%s%s%s=\"%s\"\n", "'$prefix'",vn, $2, $3);
      }
   }'
}

eval $(parse_yaml profile/variables.yaml "vars_")

if [ -z "${1}" ] || [ $# -ne 11 ]
then
   echo "give the following parameters:"
   echo "id(1) fastaURL(2) url_bam_file(3) indexURL_bam_file(bai)(4) htmlnamei(5) dir_script(pwd)(6) url_vcf_file(7) indexURL_vcf(8) url_gff3(9) indexURL_gff3(10) GC_content_bedgraph(11)"
   exit 1
else
   pwd="${6}"
   id="${1}"
   fastaURL="/${vars_Jovian_run_identifier}/${pwd}/${2}"
   url_bam="/${vars_Jovian_run_identifier}/${pwd}/${3}"
   indexURL_bam="/${vars_Jovian_run_identifier}/${pwd}/${4}"
   htmlname="${5}"
   url_vcf="/${vars_Jovian_run_identifier}/${pwd}/${7}"
   indexURL_vcf="/${vars_Jovian_run_identifier}/${pwd}/${8}"
   url_gff3="/${vars_Jovian_run_identifier}/${pwd}/${9}"
   indexURL_gff3="/${vars_Jovian_run_identifier}/${pwd}/${10}"
   GCcontentBedGraph="/${vars_Jovian_run_identifier}/${pwd}/${11}"
fi

a="$(cat << EOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1, maximum-scale=1, user-scalable=no">
    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="shortcut icon" href="https://igv.org/web/img/favicon.ico">

    <title>igv.js</title>

    <!-- IGV JS-->
    <script src="/${vars_Jovian_run_identifier}/dist/igv.min.js"></script>

</head>

<body>

<h1>Alignments from a BAM file</h1>

<div id="igv-div" style="padding-top: 10px;padding-bottom: 10px; border:1px solid lightgray"></div>

<script type="text/javascript">

    document.addEventListener("DOMContentLoaded", function () {

        var options =
        {
                reference:
                    {
                        id: "<<<id>>>",
                        fastaURL: "<<<fastaURL>>>"
                    },

            tracks: [
                {
                    type: "wig",
                    name: "GC contents",
                    format: "bedGraph",
                    url: "<<<GCcontentBedGraph>>>",
                    min: "0",
                    max: "1"
                },
                {
                    name:"SNPs",
                    type:"variant",
                    format:"vcf",
                    url: "<<<url_vcf>>>",
                    indexURL: "<<<indexURL_vcf>>>",
                },
                {
                    type: "alignment",
                    format: "bam",
                    colorBy: "strand",
                    url: "<<<url_bam>>>",
                    indexURL: "<<<indexURL_bam>>>",
                    indexed: "true",
                    name: "Alignment"
                },
                {
                    type: "annotation",
                    name: "ORF predictions",
                    format: "gff3",
                    url: "<<<url_gff3>>>",
                    indexURL: "<<<indexURL_gff3>>>",
                    displayMode: "EXPANDED"
                }
            ]
        };

        var igvDiv = document.getElementById("igv-div");

        igv.createBrowser(igvDiv, options)
                .then(function (browser) {
                    console.log("Created IGV browser");
                })

    });

</script>

</body>

</html>
EOF
)"

echo "${a}" |awk '{
       gsub ("<<<id>>>",id)
       gsub ("<<<fastaURL>>>",fastaURL)
       gsub ("<<<locus>>>",locus)
       gsub ("<<<url_bam>>>",url_bam)
       gsub ("<<<indexURL_bam>>>",indexURL_bam)
       gsub ("<<<url_vcf>>>",url_vcf)
       gsub ("<<<indexURL_vcf>>>",indexURL_vcf)
       gsub ("<<<url_gff3>>>",url_gff3)
       gsub ("<<<indexURL_gff3>>>",indexURL_gff3)
       gsub ("<<<GCcontentBedGraph>>>", GCcontentBedGraph)
       print $0
      }' id="$id" fastaURL="${fastaURL}" locus="${locus}" url_bam="${url_bam}" indexURL_bam="${indexURL_bam}" url_vcf="${url_vcf}" indexURL_vcf="${indexURL_vcf}" url_gff3="${url_gff3}" indexURL_gff3="${indexURL_gff3}" GCcontentBedGraph="${GCcontentBedGraph}" > ${htmlname}

#cp ${htmlname} /data/software/igv.js/generated/
