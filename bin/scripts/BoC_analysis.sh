#!/bin/bash

#####################################################################################################################
### Script to determine the breadth of coverage (BoC) of the reference alignment at different coverage            ###
### thresholds (1, 5, 10, 30, 100)                                                                                ###
###     Usage: bin/scripts/RA_BoC_analysis.sh sample_name sample_name_sorted.bam reference.fa \                   ###
###            sample_name_BoC.pileup sample_name_BoC_pct.tsv sample_name_BoC_int.tsv                             ###
###     Output:                                                                                                   ###
###            Sample_name     Total_ref_size  BoC_at_coverage_threshold_1     BoC_at_coverage_threshold_5 \      ###
###            BoC_at_coverage_threshold_10    BoC_at_coverage_threshold_30    BoC_at_coverage_threshold_100      ###
###            RUN25-01_S1     29904           .9999   .6666   .3333   .1000   0                                  ###
###            RUN25-02_S6     29904           .9999   .6666   .3333   .1000   0                                  ###
###            RUN25-03_S11    29904           .9999   .6666   .3333   .1000   0                                  ###
#####################################################################################################################
### TODO can probably be made prettier via a function.
### TODO parallize via parallel function (would need to be included into the env)
### TODO after the step above, make this rule a localrule in the snakefile
### TODO make the coverage thresholds a CLI argument

# Import the positional arguments (specified in the Snakefile)
INPUT_SAMPLE_NAME="$1"
INPUT_BEDGRAPH="$2"
INPUT_REF="$3"
OUTPUT_PCT_BoC="$4"
OUTPUT_INT_BoC="$5"

#? Awk 1 explanation:
####? BEGIN: set variable
####? the body: if field 4 which contains te coverage is greater than or equal to the threshold subtract the start coordinate (field 2) from the end coordinate (field 3) and sum these for all lines that abide the threshold
#? Awk 2 explanation:
####? simply a calculation to convert the integer value into a percentage of the reference. Force output with 2 decimal points via the "printf" function with "%.2f" (=2 decimal points).


#! N.B. this script only works if there are no duplicate genomes given as reference, since the sum of all reference sequences is used as $ref_total_size. I.e. adding multiple genomes of the same virus will ruin the analysis.
# Determine the total reference size:
ref_total_size=$(bowtie2-inspect -s ${INPUT_REF} | gawk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}') # The `BEGIN{L=0}` is only to switch easily between 0- and 1-count later, if needed.

# Determine the BoC at different coverage thresholds (1, 5, 10, 30 and 100) both as integer and as percentage
BoC_at_threshold_1=$(gawk -v threshold="1" 'BEGIN {covered_nt=0} {if ($4>=threshold) {block_length=$3-$2; covered_nt+=block_length;}} END { print covered_nt;}' ${INPUT_BEDGRAPH})
BoC_at_threshold_1_pct=$(gawk -v total_length="${ref_total_size}" -v counted_length="${BoC_at_threshold_1}" 'BEGIN {result = counted_length / total_length * 100; printf ("%.2f\n", result);}')
BoC_at_threshold_5=$(gawk -v threshold="5" 'BEGIN {covered_nt=0} {if ($4>=threshold) {block_length=$3-$2; covered_nt+=block_length;}} END { print covered_nt;}' ${INPUT_BEDGRAPH})
BoC_at_threshold_5_pct=$(gawk -v total_length="${ref_total_size}" -v counted_length="${BoC_at_threshold_5}" 'BEGIN {result = counted_length / total_length * 100; printf ("%.2f\n", result);}')
BoC_at_threshold_10=$(gawk -v threshold="10" 'BEGIN {covered_nt=0} {if ($4>=threshold) {block_length=$3-$2; covered_nt+=block_length;}} END { print covered_nt;}' ${INPUT_BEDGRAPH})
BoC_at_threshold_10_pct=$(gawk -v total_length="${ref_total_size}" -v counted_length="${BoC_at_threshold_10}" 'BEGIN {result = counted_length / total_length * 100; printf ("%.2f\n", result);}')
BoC_at_threshold_30=$(gawk -v threshold="30" 'BEGIN {covered_nt=0} {if ($4>=threshold) {block_length=$3-$2; covered_nt+=block_length;}} END { print covered_nt;}' ${INPUT_BEDGRAPH})
BoC_at_threshold_30_pct=$(gawk -v total_length="${ref_total_size}" -v counted_length="${BoC_at_threshold_30}" 'BEGIN {result = counted_length / total_length * 100; printf ("%.2f\n", result);}')
BoC_at_threshold_100=$(gawk -v threshold="100" 'BEGIN {covered_nt=0} {if ($4>=threshold) {block_length=$3-$2; covered_nt+=block_length;}} END { print covered_nt;}' ${INPUT_BEDGRAPH})
BoC_at_threshold_100_pct=$(gawk -v total_length="${ref_total_size}" -v counted_length="${BoC_at_threshold_100}" 'BEGIN {result = counted_length / total_length * 100; printf ("%.2f\n", result);}')

# Print to the integer and percentage BoC output files.
echo -e "${INPUT_SAMPLE_NAME}\t${ref_total_size}\t${BoC_at_threshold_1}\t${BoC_at_threshold_5}\t${BoC_at_threshold_10}\t${BoC_at_threshold_30}\t${BoC_at_threshold_100}" > ${OUTPUT_INT_BoC}
echo -e "${INPUT_SAMPLE_NAME}\t${ref_total_size}\t${BoC_at_threshold_1_pct}\t${BoC_at_threshold_5_pct}\t${BoC_at_threshold_10_pct}\t${BoC_at_threshold_30_pct}\t${BoC_at_threshold_100_pct}" > ${OUTPUT_PCT_BoC}