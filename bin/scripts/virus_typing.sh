#!/bin/bash
#####################################################################################################################
### This script interacts with the public web-based virus typingtools hosted by the RIVM:                         ###
###    Norovirus:      https://www.rivm.nl/mpf/[typingservice|typingtool]/norovirus/                              ###
###    Enterovirus:    https://www.rivm.nl/mpf/[typingservice|typingtool]/enterovirus/                            ###
###    Hepatitis A:    https://www.rivm.nl/mpf/[typingservice|typingtool]/hav/                                    ###
###    Hepatitis E:    https://www.rivm.nl/mpf/[typingservice|typingtool]/hev/                                    ###
###    Rotavirus:      https://www.rivm.nl/mpf/[typingservice|typingtool]/rotavirusa/                             ###
###    Papillomavirus: https://www.rivm.nl/mpf/[typingservice|typingtool]/papillomavirus/                         ###
###    Flavivirus:     https://www.rivm.nl/mpf/[typingservice|typingtool]/flavivirus/                             ###
###                                                                                                               ###
### Usage: bin/virus_typing.sh {NoV|EV|HAV|HEV|RVA|PV|Flavi} (--force)                                            ###
###     --force     Will redo and force-overwrite previously generated results.                                   ###
###     --help      Will print the help message and exit.                                                         ###
#####################################################################################################################

# Setup
INPUT_FOLDER="data/tables/" # Default Jovian output directory ("data/tables/")
INPUT_FILES="${INPUT_FOLDER}*_taxClassified.tsv" # Combine the INPUT_FOLDER with the default output file suffix ("*_taxClassified.tsv") into a file glob
OUTPUT_FOLDER="data/virus_typing_tables/"
mkdir -p ${OUTPUT_FOLDER}

WHICH_TT=${1,,} # The ${var,,} syntax converts it all to lowercase, this makes it easier to regex/match the keywords later and makes the keywords case-insensitive.

usage_msg() {
    cat <<HELP_USAGE
Usage: bash jovian -vt [virus-keyword]

Mandatory arguments:
[virus-keyword]                         A virus keyword as specified below.

Optional arguments:
bash jovian -vt-force [virus-keyword]   Force overwrite existing output.
bash jovian -vt-help                    Print this help message.

The first argument should always be one of the keyword listed below:
    ---------------------------------------------------
    | KEYWORD | TYPEABLE VIRUSSES                     |
    |-------------------------------------------------|
    | all     | perform all the different typing-tools|
    |         | listed below.                         |
    |-------------------------------------------------|
    | NoV     | Caliciviridae family;                 |
    |         | - Norwalk virus (GI & GII),           |
    |         | - Sapporo virus                       |
    |-------------------------------------------------|
    | EV      | Picornaviridae family;                |
    |         | - Aichivirus A-C                      |
    |         | - Cosavirus A,B,D,E                   |
    |         | - Enterovirus A-H, J                  |
    |         |   - Incl. Enterovirus, Aichivirus,    |
    |         |     Human Poliovirus, Echovirus       |
    |         | - Human Parechovirus A,B              |
    |         | - Rhinovirus A,B,C                    |
    |-------------------------------------------------|
    | RVA     | Rotavirus genus;                      |
    |         | - Rotavirus A                         |
    |-------------------------------------------------|
    | HAV     | Hepatovirus genus;                    |
    |         | - Hepatovirus A                       |
    |-------------------------------------------------|
    | HEV     | Orthohepevirus genus;                 |
    |         | - Orthohepevirus A (Hepatitis E)      |
    |-------------------------------------------------|
    | PV      | Papillomaviridae family;              |
    |         | - Alphapapillomaviruses,              |
    |         | - Betapapillomaviruses,               |
    |         | - Gammapapillomaviruses               |
    |-------------------------------------------------|
    | Flavi   | Flaviviridae family;                  |
    |         | - Dengue virus 1-4                    |
    |         | - Hepacivirus C (Hepatitis C)         |
    |         | - Pegivirus C (GB-C)                  |
    |         | - Tick-borne encephalitis virus       |
    |         | - Zika virus                          |
    ---------------------------------------------------
HELP_USAGE
}

#################################################################################
### Parse CLI argument, throw error if wrong, set proper variables if right #####
#################################################################################

wrong_tt_keyword_err_msg() {
    echo -e "\nUnknown parameter '"${1}"' given.\n"
    usage_msg
}

validate_input_tt_keyword() {
   if [[ "${1}" =~ ^(nov|ev|hav|hev|rva|pv|flavi|all)$ ]]; then
        if [ ! -d "${INPUT_FOLDER}" ]; then
            echo -e "No '"${INPUT_FOLDER}"' folder found. Virus typing can only be performed after a completed Jovian analysis."
            exit 1
        fi
    else
        wrong_tt_keyword_err_msg "${1}"
        exit 1
    fi
}

#! If on any positional argument `--help` is given, print help message
if [[ "$*" == *--help* ]]; then
    usage_msg
    exit 0
#! If $1 is not empty, and $2 is not empty, and number of arguments is equal to 2 --> so the check if a virus keyword was given with the --force keyword
elif [ ! -z "${1}" -a ! -z "${2}" -a $# -eq 2 ]; then
    validate_input_tt_keyword "${WHICH_TT}"
    #! If $2 equal the literal --force
    if [ "${2}" == "--force" ]; then
        FORCE_FLAG="TRUE"
    else
        echo -e "\nUnknown parameter '"${2}"' given. Did you mean '"--force"'?"
        exit 1
    fi
#! If $1 is not empty, and number of arguments is equal to 1 --> so the check if only a virus keyword was given without --force keyword
elif [ ! -z "${1}" -a $# -eq 1 ]; then
    validate_input_tt_keyword "${WHICH_TT}"
else
    echo -e "Invalid input parameters, please use one of these parameters:\n"
    usage_msg
    exit 1
fi

#################################################################################
### Virus typing                                                            #####
#################################################################################

extract_fasta() {
    local input="${1}"
    local output="${2}"
    local capture_name="${3}"
    local capture_field="${4}"
    # Extract the scaffold name and sequence of a certain taxonomic rank from the complete Jovian taxonomic output and write it as a fasta
    gawk -F "\t" -v name="${capture_name}" -v field="${capture_field}" '$field == name {print ">" $2 "\n" $25}' < ${input} > ${output}
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
    if [ "${which_tt}" == "nov" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/norovirus/"
        local parser_py="bin/scripts/typingtool_NoV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_NoV.fa}
        local extract_name="Caliciviridae" # Family
        local extract_field="8" # Family
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with family == Caliciviridae found."
    elif [ "${which_tt}" == "ev" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/enterovirus/"
        local parser_py="bin/scripts/typingtool_EV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_EV.fa}
        local extract_name="Picornaviridae" # Family
        local extract_field="8" # Family
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with family == Picornaviridae found."
    elif [ "${which_tt}" == "hav" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/hav/"
        local parser_py="bin/scripts/typingtool_HAV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_HAV.fa}
        local extract_name="Hepatovirus" # Genus
        local extract_field="7" # Genus
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with genus == Hepatovirus found."
    elif [ "${which_tt}" == "hev" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/hev/"
        local parser_py="bin/scripts/typingtool_HEV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_HEV.fa}
        local extract_name="Orthohepevirus" # Genus
        local extract_field="7" # Genus
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with genus == Orthohepevirus found."
    elif [ "${which_tt}" == "rva" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/rotavirusa/"
        local parser_py="bin/scripts/typingtool_RVA_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_RVA.fa}
        local extract_name="Rotavirus" # Genus
        local extract_field="7" # Genus
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with genus == Rotavirus found."
    elif [ "${which_tt}" == "pv" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/papillomavirus/"
        local parser_py="bin/scripts/typingtool_PV_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_PV.fa}
        local extract_name="Papillomaviridae" # Family
        local extract_field="8" # Family
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with family == Papillomaviridae found."
    elif [ "${which_tt}" == "flavi" ]; then
        local tt_url="https://www.rivm.nl/mpf/typingservice/flavivirus/"
        local parser_py="bin/scripts/typingtool_Flavi_XML_to_csv_parser.py"
        local query_fasta=${OUTPUT_FOLDER}${basename/_taxClassified.tsv/_Flavi.fa}
        local extract_name="Flaviviridae" # Family
        local extract_field="8" # Family
        local nothing_found_message="Sample:\t${sample_name}\tNo scaffolds with family == Flaviviridae found."
    else
        echo -e "Unknown typingtool specified."
        wrong_tt_keyword_err_msg "${1}"
        exit 1
    fi

    local tt_xml=${query_fasta/.fa/.xml}
    local tt_csv=${tt_xml/.xml/.csv}

    #! If tt_csv doesn't exist, OR, FORCE_FLAG is NOT empty (i.e. force overwrite previously generated output)
    #TODO this will need to be changed if we remove the tt_csv output in the snakemake script (remove temp chunk/onsuccess)
    if [ ! -e "${tt_csv}" ] || [ ! -z "${FORCE_FLAG}" ]; then
        # Extract taxonomic slice fasta, send to TT, parse the results XML into csv
        extract_fasta "${file_path}" "${query_fasta}" "${extract_name}" "${extract_field}"
        if [ -s "${query_fasta}" ]
        then
            echo -e "Sample:\t${sample_name}\tScaffolds compatible with the ${which_tt} tool found, sent to typingtool service, waiting for results... This may take a while..."
            submit_query_fasta "${query_fasta}" "${tt_xml}" "${tt_url}"

            # Sadly, the current version of the typingtool service has some issues, resulting in errors because it can't handle the request.
            ### One of two things can happen; you get a terse xml output that states "502 Proxy Error" or you get a verbose html output that states "Error reading from remote server". Hence, the double grep OR statement...
            if grep -q -e "502 Proxy Error" -e "Error reading from remote server" ${tt_xml}; then
                echo -e "Sample:\t${sample_name}\tQuery cannot currently be handled by typingtool... Please try again later, for further information, see: https://github.com/DennisSchmitz/Jovian/issues/51"
            else
                # If error code is not found; parse the XML
                python ${parser_py} "${sample_name}" "${tt_xml}" "${tt_csv}"
            fi
        else
            echo -e "${nothing_found_message}"
        fi
    #! If tt_csv is not empty (i.e. it exists and has content), AND, FORCE_FLAG is empty (i.e. don't force overwrite previously generated output)
    elif [ -s "${tt_csv}" ] && [ -z "${FORCE_FLAG}" ]; then
        echo -e "Sample:\t${sample_name}\tScaffolds compatible with the ${which_tt} tool were already found and analyzed in earlier analysis. Skipping..."
    else
        echo -e "${nothing_found_message}"
    fi
}

# Perform all typingtool functions for each input file specified in the glob below, based on the INPUT_FILES variable
echo -e "\nStarting with ${WHICH_TT} typingtool analysis.\nN.B. depending on the size of your dataset, and the load of the virus typingtool webservice, this might take some time...\n"
for FILE in ${INPUT_FILES}
do
    BASENAME=${FILE##*/}   # Filename without path but WITH suffixes
    if [ "${WHICH_TT}" == "all" ]; then
        #TODO make this more elegant with a list of keywords
        typingtool "${FILE}" "${BASENAME}" "nov"
        typingtool "${FILE}" "${BASENAME}" "ev"
        typingtool "${FILE}" "${BASENAME}" "hav"
        typingtool "${FILE}" "${BASENAME}" "hev"
        typingtool "${FILE}" "${BASENAME}" "rva"
        typingtool "${FILE}" "${BASENAME}" "pv"
        typingtool "${FILE}" "${BASENAME}" "flavi"
    else
        typingtool "${FILE}" "${BASENAME}" "${WHICH_TT}"
    fi
done

#################################################################################
### Concatenate all indivual files together into one big file               #####
#################################################################################

if [ -n "$( find data/virus_typing_tables/ -maxdepth 1 -name "*_${WHICH_TT}.csv" -print -quit )" ]
then
    # If any files were created in the first place; concat individual outputs into one combined output, the awk magic is to not repeat headers
    gawk 'FNR==1 && NR!=1 { next; } { print }' data/virus_typing_tables/*_${WHICH_TT}.csv > results/all_${WHICH_TT}-TT.csv
fi

#################################################################################
### Cleanup #####
#################################################################################
find data/virus_typing_tables/ -type f -empty -delete
#TODO rm -f data/virus_typing_tables/*_${WHICH_TT}.fa # Commented this out for debugging purposes, should be activated in v.1.0 (but please first see other #TODO above)
#TODO rm -f data/virus_typing_tables/*_${WHICH_TT}.xml # Commented this out for debugging purposes, should be activated in v.1.0 (but please first see other #TODO above)
echo -e "\nFinished"