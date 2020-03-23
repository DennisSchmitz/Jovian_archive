#!/bin/bash

#####################################################################################################################
### Script to determine the breadth of coverage (BoC) of the reference alignment at different coverage            ###
### thresholds (1, 5, 10, 30, 100)                                                                                ###
###     Usage: bin/scripts/RA_BoC_analysis.sh sample_name sample_name_sorted.bam reference.fa \                   ###
###            sample_name_BoC.pileup sample_name_BoC_pct.tsv sample_name_BoC_int.tsv                             ###
###     Output:                                                                                                   ###
###            Sample_name     Total_ref_size  BoC_at_coverage_threshold_1     BoC_at_coverage_threshold_5 \      ###
###            BoC_at_coverage_threshold_10    BoC_at_coverage_threshold_30    BoC_at_coverage_threshold_100      ###
###            RUN25-01_S1     29904           .9998   .9993   .9993   .9990   .9989                              ###
###            RUN25-02_S6     29904           .9182   .8351   .8070   .4027   .4406                              ###
###            RUN25-03_S11    29904           .8636   .7477   .5843   .0388   0                                  ###
#####################################################################################################################
### TODO can probably be made prettier via a function.
### TODO make the coverage thresholds a CLI argument

# Import the positional arguments (specified in the Snakefile)
INPUT_SAMPLE_NAME="$1"
INPUT_BEDGRAPH="$2"
INPUT_REF="$3"
OUTPUT_PCT_BoC="$4"
OUTPUT_INT_BoC="$5"

#! N.B. this script only works if there are no duplicate genomes given as reference, since the sum of all reference sequences is used as $ref_total_size. I.e. adding multiple genomes of the same virus will ruin the analysis.
# Determine the total reference size:
#? Bowtie2-inspect reports a 0-start ref length, but since the bed file is also 0-start it doesn't matter.
ref_total_size=$(bowtie2-inspect -s ${INPUT_REF} | gawk '{ FS = "\t" } ; BEGIN{L=1}; {L=L+$3}; END{print L}') # The `BEGIN{L=1}` is because the bowtie2-inspect length value starts counting at position 0 instead of 1, will be confusing for end-users.

# Determine the BoC at different coverage thresholds (1, 5, 10, 30 and 100) both as integer and as percentage
BoC_at_threshold_1=$(gawk -v threshold="1" 'BEGIN {covered_nt=1} {if ($4>=threshold) block_length=$3-$2; covered_nt+=block_length;} END { print covered_nt }' ${INPUT_BEDGRAPH}) # The `BEGIN {covered_nt=1} ` is because the bedgraph format starts counting at position 0 instead of 1, will be confusing for end-users.
BoC_at_threshold_1_pct=$(echo "scale=4 ; ${BoC_at_threshold_1} / ${ref_total_size} " | bc )
BoC_at_threshold_5=$(gawk -v threshold="5" 'BEGIN {covered_nt=1} {if ($4>=threshold) block_length=$3-$2; covered_nt+=block_length;} END { print covered_nt }' ${INPUT_BEDGRAPH}) # The `BEGIN {covered_nt=1} ` is because the bedgraph format starts counting at position 0 instead of 1, will be confusing for end-users.
BoC_at_threshold_5_pct=$(echo "scale=4 ; ${BoC_at_threshold_5} / ${ref_total_size} " | bc )
BoC_at_threshold_10=$(gawk -v threshold="10" 'BEGIN {covered_nt=1} {if ($4>=threshold) block_length=$3-$2; covered_nt+=block_length;} END { print covered_nt }' ${INPUT_BEDGRAPH}) # The `BEGIN {covered_nt=1} ` is because the bedgraph format starts counting at position 0 instead of 1, will be confusing for end-users.
BoC_at_threshold_10_pct=$(echo "scale=4 ; ${BoC_at_threshold_10} / ${ref_total_size} " | bc )
BoC_at_threshold_30=$(gawk -v threshold="30" 'BEGIN {covered_nt=1} {if ($4>=threshold) block_length=$3-$2; covered_nt+=block_length;} END { print covered_nt }' ${INPUT_BEDGRAPH}) # The `BEGIN {covered_nt=1} ` is because the bedgraph format starts counting at position 0 instead of 1, will be confusing for end-users.
BoC_at_threshold_30_pct=$(echo "scale=4 ; ${BoC_at_threshold_30} / ${ref_total_size} " | bc )
BoC_at_threshold_100=$(gawk -v threshold="100" 'BEGIN {covered_nt=1} {if ($4>=threshold) block_length=$3-$2; covered_nt+=block_length;} END { print covered_nt }' ${INPUT_BEDGRAPH}) # The `BEGIN {covered_nt=1} ` is because the bedgraph format starts counting at position 0 instead of 1, will be confusing for end-users.
BoC_at_threshold_100_pct=$(echo "scale=4 ; ${BoC_at_threshold_100} / ${ref_total_size} " | bc )

# Print to the integer and percentage BoC output files.
echo -e "${INPUT_SAMPLE_NAME}\t${ref_total_size}\t${BoC_at_threshold_1}\t${BoC_at_threshold_5}\t${BoC_at_threshold_10}\t${BoC_at_threshold_30}\t${BoC_at_threshold_100}" > ${OUTPUT_INT_BoC}
echo -e "${INPUT_SAMPLE_NAME}\t${ref_total_size}\t${BoC_at_threshold_1_pct}\t${BoC_at_threshold_5_pct}\t${BoC_at_threshold_10_pct}\t${BoC_at_threshold_30_pct}\t${BoC_at_threshold_100_pct}" > ${OUTPUT_PCT_BoC}