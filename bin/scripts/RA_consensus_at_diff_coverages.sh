#!/bin/bash

#####################################################################################################################
### Script to extract the consensus greater or equal to different coverage thresholds (1, 5, 10, 30, 100)         ###
###     Usage: bin/scripts/RA_consensus_at_diff_coverages.sh sample_name sample_name_sorted.bam reference.fasta \ ###
###            sample_name_raw_consensus.fa output/dir/ sample_name.log                                           ###
###     Output: Per coverage threshold (1, 5, 10, 30, 100), two fasta files are generated. One where the failing  ###
###             coordinates are replaced with N nucleotides (ambiguous) and another replaced by the "-" character ###
###             (gap). This is to let the end-user decide which fasta best suits their downstream program.        ###
###             N.B. in the "-" filtered fasta you cannot distinguish real gaps from filtered regions!            ###
#####################################################################################################################

#####################################################################################################################
### Import the positional arguments (specified in the Snakefile)                                                  ###
#####################################################################################################################
INPUT_SAMPLE_NAME="$1"
INPUT_BAM="$2"
INPUT_REFERENCE="$3"
INPUT_RAW_CONSENSUS="$4"
OUTPUT_FOLDER_DATA="$5"
OUTPUT_FOLDER_RESULTS="$6"
LOG_FILE="$7"

OUTPUT_BASENAME_DATA="${OUTPUT_FOLDER_DATA}/${INPUT_SAMPLE_NAME}"
OUTPUT_BEDGRAPH="${OUTPUT_BASENAME_DATA}.bedgraph"
OUTPUT_BASENAME_RESULTS="${OUTPUT_FOLDER_RESULTS}/${INPUT_SAMPLE_NAME}"

# Make the empty output dir:
mkdir -p ${OUTPUT_FOLDER_RESULTS}

#####################################################################################################################
### Define function                                                                                               ###
#####################################################################################################################

function consensus_at_variable_cov_threshold {
    local coverage_threshold="$1"
    local output_bed="${OUTPUT_BASENAME_DATA}_cov_ge_${coverage_threshold}.bed"
    local consensus_N_filt_tmp="${OUTPUT_BASENAME_RESULTS}_N-filt_cov_ge_${coverage_threshold}.tmp"
    local consensus_N_filt="${OUTPUT_BASENAME_RESULTS}_N-filt_cov_ge_${coverage_threshold}.fa"
    local consensus_minus_filt_tmp="${OUTPUT_BASENAME_RESULTS}_minus-filt_cov_ge_${coverage_threshold}.tmp"
    local consensus_minus_filt="${OUTPUT_BASENAME_RESULTS}_minus-filt_cov_ge_${coverage_threshold}.fa"

    # From the bedgraph, extract the regions that do NOT match the coverage threshold. These will be masked/filtered later on.
    gawk -v threshold="${coverage_threshold}" '$4 < threshold' ${OUTPUT_BEDGRAPH} > ${output_bed}

    if [ -s "${output_bed}" ]
    then #? bedfile is not empty, hence there are regions to be filtered
        # Mask the nucleotide positions in the genome that do not reach the coverage threshold based on the previously made .bed file which contains the failing nucleotide positions
        #### Importantly, the default setting is to replace the bed coordinates with N nucleotides, hence this generates the N nucleotide filtered consensus
        bedtools maskfasta -fullHeader -fi ${INPUT_REFERENCE} -bed ${output_bed} -fo ${consensus_N_filt_tmp} 2>> ${LOG_FILE}
        # Convert the multi-line fasta output of bedtools maskfasta to the normal two-line fasta output.
        seqtk seq ${consensus_N_filt_tmp} > ${consensus_N_filt} 2>> ${LOG_FILE}
        rm ${consensus_N_filt_tmp}
        # Prefix the sample name, filtering character and cov-threshold information to the fasta header
        sed -i "s/^>/>${INPUT_SAMPLE_NAME}_N-filt_cov-ge-${coverage_threshold}_/g" ${consensus_N_filt} 2>> ${LOG_FILE}

        # Mask the nucleotide positions in the genome that do not reach the coverage threshold based on the previously made .bed file which contains the failing nucleotide positions
        #### Importantly, here the bed coordinates are replaced with the minus ("-") character which represents a gap in alignment software. While it is NOT an actual gap, downstream
        #### alignment software is usually more compatible with gap characters than the default N characters. This was implemented on end-user request, they can decide which version
        #### of the output they want to use for downstream processing.
        bedtools maskfasta -fullHeader -mc "-" -fi ${INPUT_REFERENCE} -bed ${output_bed} -fo ${consensus_minus_filt_tmp} 2>> ${LOG_FILE}
        # Convert the multi-line fasta output of bedtools maskfasta to the normal two-line fasta output.
        seqtk seq ${consensus_minus_filt_tmp} > ${consensus_minus_filt} 2>> ${LOG_FILE}
        rm ${consensus_minus_filt_tmp}
        # Prefix the sample name, filtering character and cov-threshold information to the fasta header
        sed -i "s/^>/>${INPUT_SAMPLE_NAME}_minus-filt_cov-ge-${coverage_threshold}_/g" ${consensus_minus_filt} 2>> ${LOG_FILE}
    else #? bedfile is empty, hence there are no regions to be filtered and the RAW consensus can be considered suitable for downstream processes
        cp ${INPUT_RAW_CONSENSUS} ${consensus_N_filt} >> ${LOG_FILE} 2>&1
        # Prefix the sample name, filtering character and cov-threshold information to the fasta header
        sed -i "s/^>/>${INPUT_SAMPLE_NAME}_N-filt_cov-ge-${coverage_threshold}_/g" ${consensus_N_filt} 2>> ${LOG_FILE}

        cp ${INPUT_RAW_CONSENSUS} ${consensus_minus_filt} >> ${LOG_FILE} 2>&1
        # Prefix the sample name, filtering character and cov-threshold information to the fasta header
        sed -i "s/^>/>${INPUT_SAMPLE_NAME}_minus-filt_cov-ge-${coverage_threshold}_/g" ${consensus_minus_filt} 2>> ${LOG_FILE}
    fi
}

#####################################################################################################################
### Actual code                                                                                                   ###
#####################################################################################################################

# First create a genome-coverage bedgraph that will be used as input for the downstream processes.
bedtools genomecov -bga -ibam ${INPUT_BAM} > ${OUTPUT_BEDGRAPH} 2>> ${LOG_FILE}

# Generate consensuses at different coverage thresholds
consensus_at_variable_cov_threshold "1"
consensus_at_variable_cov_threshold "5"
consensus_at_variable_cov_threshold "10"
consensus_at_variable_cov_threshold "30"
consensus_at_variable_cov_threshold "100"