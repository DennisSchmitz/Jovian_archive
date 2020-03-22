#!/bin/bash

#####################################################################################################################
### Script to determine the breadth of coverage (BoC) of the reference alignment at different coverage            ###
### thresholds (1, 5, 10, 30, 100)                                                                                ###
###     Usage: bin/scripts/RA_BoC_analysis.sh sample_name sample_name_sorted.bam reference.fa \                   ###
###            sample_name_BoC.pileup sample_name_BoC_pct.tsv sample_name_BoC_int.tsv                             ###
###     Output:                                                                                                   ###
###            Sample_name     Total_ref_size  BoC_at_coverage_threshold_1     BoC_at_coverage_threshold_5 \      ###
###            BoC_at_coverage_threshold_10    BoC_at_coverage_threshold_30    BoC_at_coverage_threshold_100      ###
###            RUN25-01_S1     29903           .9990   .9987   .9981   .9977                                      ###
###            RUN25-02_S6     29903           .1248   .0198   0       0                                          ###
###            RUN25-03_S11    29903           .0029   0       0       0                                          ###
###     Source: http://www.metagenomics.wiki/tools/samtools/breadth-of-coverage                                   ###
#####################################################################################################################
### TODO can probably be made prettier via a function.
### TODO make the coverage thresholds a CLI argument

# Import the positional arguments (specified in the Snakefile)
INPUT_SAMPLE_NAME="$1"
INPUT_BAM="$2"
INPUT_REF="$3"
OUTPUT_PILEUP="$4"
OUTPUT_PCT_BoC="$5"
OUTPUT_INT_BoC="$6"

# Generate a pileup file (https://en.wikipedia.org/wiki/Pileup_format)
samtools mpileup -f ${INPUT_REF}  ${INPUT_BAM} -o ${OUTPUT_PILEUP}


#! N.B. this script only works if there are no duplicate genomes given as reference, since the sum of all reference sequences is used as $ref_total_size. I.e. adding multiple genomes of the same virus will ruin the analysis.
# Determine the total reference size:
ref_total_size=$(bowtie2-inspect -s ${INPUT_REF} | gawk '{ FS = "\t" } ; BEGIN{L=0}; {L=L+$3}; END{print L}')

# Determine the BoC at different coverage thresholds (1, 5, 10, 30 and 100) both as integer and as percentage
BoC_at_threshold_1=$(cat ${OUTPUT_PILEUP} | gawk -v X="1" '$4>=X' | wc -l)
BoC_at_threshold_1_pct=$(echo "scale=4 ; ${BoC_at_threshold_1} / ${ref_total_size} " | bc )
BoC_at_threshold_5=$(cat ${OUTPUT_PILEUP} | gawk -v X="5" '$4>=X' | wc -l)
BoC_at_threshold_5_pct=$(echo "scale=4 ; ${BoC_at_threshold_5} / ${ref_total_size} " | bc )
BoC_at_threshold_10=$(cat ${OUTPUT_PILEUP} | gawk -v X="10" '$4>=X' | wc -l)
BoC_at_threshold_10_pct=$(echo "scale=4 ; ${BoC_at_threshold_10} / ${ref_total_size} " | bc )
BoC_at_threshold_30=$(cat ${OUTPUT_PILEUP} | gawk -v X="30" '$4>=X' | wc -l)
BoC_at_threshold_30_pct=$(echo "scale=4 ; ${BoC_at_threshold_30} / ${ref_total_size} " | bc )
BoC_at_threshold_100=$(cat ${OUTPUT_PILEUP} | gawk -v X="100" '$4>=X' | wc -l)
BoC_at_threshold_100_pct=$(echo "scale=4 ; ${BoC_at_threshold_100} / ${ref_total_size} " | bc )

# Print to the integer and percentage BoC output files.
echo -e "${INPUT_SAMPLE_NAME}\t${ref_total_size}\t${BoC_at_threshold_1}\t${BoC_at_threshold_5}\t${BoC_at_threshold_10}\t${BoC_at_threshold_30}\t${BoC_at_threshold_100}" > ${OUTPUT_INT_BoC}
echo -e "${INPUT_SAMPLE_NAME}\t${ref_total_size}\t${BoC_at_threshold_1_pct}\t${BoC_at_threshold_5_pct}\t${BoC_at_threshold_10_pct}\t${BoC_at_threshold_30_pct}\t${BoC_at_threshold_100_pct}" > ${OUTPUT_PCT_BoC}