#!/usr/bin/env python3
"""Parser for the XML output of the papillomaviridae (PV) typing tool (TT).

Usage:
  typingtool_PV_XML_to_csv_parser.py <sample_name> <input_PV_XML> <output_PV_csv>

<sample_name> is how you want to call this sample (in plain text).
<input_PV_XML> is the input XML file that you got from the PV TT.
<output_PV_csv> is the output csv file.

Example:
  typingtool_PV_XML_to_csv_parser.py [sample_name] \
    data/virus_typing_tables/[sample_name]_PV.xml \
    data/virus_typing_tables/[sample_name]_PV.csv
"""

import csv
import xml.etree.cElementTree as et
from collections import OrderedDict
from sys import argv

SCRIPT, SAMPLE_NAME, INPUT_XML, OUTPUT_CSV = argv

parsedXML = et.parse(INPUT_XML)

csv_data = []

input_fields = ["result", "start", "end", "nucleotides", "conclusion"]

output_fields = [
    "Sample_name",
    "Query_name",
    "start",
    "length",
    "end",
    "blast_concluded-name",
    "blast_absolute-similarity",
    "blast_refseq",
    "blast_reverse-complement",
    "cluster_result",
    "cluster_support",
    "nucleotides",
]

for elem in parsedXML.findall(".//sequence"):
    inner_dict = OrderedDict({k: None for k in output_fields})

    inner_dict["Sample_name"] = SAMPLE_NAME

    inner_dict["Query_name"] = elem.attrib["name"]
    inner_dict["length"] = elem.attrib["length"]

    for item in elem.findall(".//*"):
        if item.tag in input_fields:
            if item.tag == "nucleotides":
                inner_dict["nucleotides"] = item.text.strip()
            elif item.tag == "result":
                tmp_result_id = item.attrib.get("id")
                if tmp_result_id == "blast":
                    inner_dict["blast_concluded-name"] = item.find(
                        ".//concluded-name"
                    ).text
                    inner_dict["blast_absolute-similarity"] = item.find(
                        ".//absolute-similarity"
                    ).text
                    try:
                        inner_dict["blast_refseq"] = item.find(".//refseq").text
                    except:
                        inner_dict["blast_refseq"] = item.find(".//refseq")
                    inner_dict["blast_reverse-complement"] = item.find(
                        ".//reverse-compliment"
                    ).text
            elif item.tag == "conclusion":
                tmp_conclusion_type = item.attrib.get("type")
                tmp_conclusion_id = item.attrib.get("id")
                tmp_conclusion_region = item.attrib.get("region")
                if tmp_conclusion_type == "simple" and tmp_conclusion_id == "type":
                    inner_dict["cluster_result"] = item.find(".//name").text
                    inner_dict["cluster_support"] = item.find(".//support").text
            else:
                inner_dict[item.tag] = item.text
    csv_data.append(inner_dict)

with open(OUTPUT_CSV, "w", newline="") as f:
    w = csv.DictWriter(f, csv_data[0].keys())
    w.writeheader()
    w.writerows(csv_data)
