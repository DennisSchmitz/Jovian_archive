#!/bin/bash
#####################################################################################################################
### This script interacts with the public web-based virus typingtools hosted by the RIVM:                         ###
###    Norovirus:   https://www.rivm.nl/mpf/typingservice/norovirus/                                              ###
###    Enterovirus: https://www.rivm.nl/mpf/typingservice/enterovirus/                                            ###
###    Hepatitis A: https://www.rivm.nl/mpf/typingservice/hav/                                                    ###
###    Hepatitis E: https://www.rivm.nl/mpf/typingservice/hev/                                                    ###
###    Rotavirus:   [work in progress, coming soon]                                                               ###
###                                                                                                               ###
### Usage: bin/fastqc_wrapper.sh {NoV|EV|HAV|HEV}                                                                 ###
#####################################################################################################################

# Check commandline argument, throw error if wrong, parse argument if right
if [ -z "${1}" ] || [ $# -ne 1 ]
then
   echo "Please specify which typingtool you want to use by entering the following parameter:"
   echo -e "\t'NoV' for Norovirus"
   echo -e "\t'EV' for Enterovirus"
   echo -e "\t'HAV' for Hepatitis A"
   echo -e "\t'HEV' for Hepatitis E"
   echo -e "Please note, these parameters are case-sensitive."
   exit 1
else
   WHICH_TT="${1}"
fi

# Setup
OUTPUT_FOLDER="data/virus_typing_tables/"
mkdir -p ${OUTPUT_FOLDER}

# Functions
extract_fasta() {
    local input="${1}"
    local output="${2}"
    local capture_name="${3}"
    local capture_field="${4}"
    # Extract the scaffold name and sequence of a certain taxonomic rank from the complete Jovian taxonomic output and write it as a fasta
    gawk -F "\t" -v name="${capture_name}" -v field="${capture_field}" '$field == name {print ">" $2 "\n" $24}' < ${input} > ${output}
}
submit_query_fasta() {
    local input="${1}"
    local output="${2}"
    local url="${3}"
    # Send the extracted taxonomic slice fasta to the specified public typing tool and wait for the XML results
    curl -s --data-urlencode fasta-sequence@${input} ${url} > ${output}
}
typingtool() {
    local file_path="${1}"
    local basename="${2}"
    local which_tt="${3}"
    local sample_name=${basename/_taxClassified.tsv/}   # Base sample name without path and suffixes

    # Set proper variables depending on chosen typingtool (either 'NoV', 'EV', 'HAV' or 'HEV')
    if [ "${which_tt}" == "NoV" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/norovirus/"
        local parser_py="bin/typingtool_NoV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_NoV.fa}
        local extract_name="Caliciviridae" # Family
        local extract_field="8" # Family
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with species == Norwalk virus found."
    elif [ "${which_tt}" == "EV" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/enterovirus/"
        local parser_py="bin/typingtool_EV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_EV.fa}
        local extract_name="Picornaviridae" # Family
        local extract_field="8" # Family
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with family == Picornaviridae found."
    elif [ "${which_tt}" == "HAV" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/hav/"
        local parser_py="bin/typingtool_HAV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_HAV.fa}
        local extract_name="Hepatovirus" # Genus
        local extract_field="7" # Genus
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with genus == Hepatovirus found."
    elif [ "${which_tt}" == "HEV" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/hev/"
        local parser_py="bin/typingtool_HEV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_HEV.fa}
        local extract_name="Orthohepevirus" # Genus
        local extract_field="7" # Genus
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with genus == Orthohepevirus found."
    else
        echo -e "Unknown typingtool specified, please specify either 'NoV', 'EV', 'HAV' or 'HEV'"
        exit 1
    fi

    local tt_xml=${query_fasta/.fa/.xml}
    local tt_csv=${tt_xml/.xml/.csv}

    # Extract taxonomic slice fasta, send to TT, parse the results XML into csv
    extract_fasta "${file_path}" "${query_fasta}" "${extract_name}" "${extract_field}"
    if [ -s "${query_fasta}" ]
    then
        echo -e "Sample:\t${sample_name}\t${which_tt} scaffolds found and sent to typingtool, waiting for results... This may take a while..."
        submit_query_fasta "${query_fasta}" "${tt_xml}" "${tt_url}"
        python ${parser_py} "${sample_name}" "${tt_xml}" "${tt_csv}"
    else
        echo -e "${nothing_found_message}"
    fi
}

# Perform all typingtool functions for each input file in the glob below (standard Jovian output)
echo -e "\nStarting with ${WHICH_TT} typingtool analysis.\nN.B. depending on the size of your dataset, and the load of the virus typingtool webservice, this might take some time...\n"
for FILE in data/tables/*_taxClassified.tsv
do
    BASENAME=${FILE##*/}   # Filename without path but WITH suffixes
    typingtool "${FILE}" "${BASENAME}" "${WHICH_TT}"
done

if test -n "$(find data/virus_typing_tables/ -maxdepth 1 -name "*_${WHICH_TT}.csv" -print -quit)"
then
    # If any files were created in the first place; concat individual outputs into one combined output, the awk magic is to not repeat headers
    gawk 'FNR==1 && NR!=1 { next; } { print }' data/virus_typing_tables/*_${WHICH_TT}.csv > results/all_${WHICH_TT}-TT.csv
fi

# Cleanup
find data/virus_typing_tables/ -type f -empty -delete
#rm -f data/virus_typing_tables/*_${WHICH_TT}.fa # Commented this out for debugging purposes, should be activated in v.1.0
#rm -f data/virus_typing_tables/*_${WHICH_TT}.xml # Commented this out for debugging purposes, should be activated in v.1.0
echo -e "\nFinished"