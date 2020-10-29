#!/bin/bash
#####################################################################################################################
### This is a wrapper for fastqc used in the snakefile.                                                           ###
###     Usage: bin/fastqc_wrapper.sh {input} {params.output_dir} {output.html} {output.zip} {log}                 ###
###     Reason: Fastqc implicitly generates output files based on the input file name. Therefore, if these are    ###
###             not the same, this output is moved to the explicitly Snakefile defined output location.           ###
###     I had to put this simple bash code into a script because of the Snakemake-bash strict settings, see:      ###
### https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-with-errors-about-an-unbound-variable-what-s-wrong
#####################################################################################################################

# Import the positional arguments (specified in the Snakefile)
INPUT_FILE="$1"
OUTPUT_DIR="$2"
DESIRED_OUTPUT_HTML="$3"
DESIRED_OUTPUT_ZIP="$4"
LOG="$5"

# Generate sample basename (i.e. remove everything before the last '/' and then remove everything after the '.')
EXTENSION1=".gz"
EXTENSION2=".fastq"
EXTENSION3=".fq"


SAMPLE_TEMP="${INPUT_FILE##*/}"
TEMP1=${SAMPLE_TEMP//$EXTENSION1/}
TEMP2=${TEMP1//$EXTENSION2/}
TEMP3=${TEMP2//$EXTENSION3/}

SAMPLE_NAME="${TEMP3}"

# These are the output files that fastqc implicitly generates (based on the input file basename)
REAL_OUTPUT_HTML="${OUTPUT_DIR}${SAMPLE_NAME}_fastqc.html"
REAL_OUTPUT_ZIP="${OUTPUT_DIR}${SAMPLE_NAME}_fastqc.zip"

# Run fastqc
fastqc -t 6 --quiet --outdir ${OUTPUT_DIR} ${INPUT_FILE} > ${LOG} 2>&1

# If the implicit output file are not equal to the explicitly defined output in the Snakefile, rename them to the explicitly defined output name.
if [ "$DESIRED_OUTPUT_HTML" != "$REAL_OUTPUT_HTML" ]
then
    mv $REAL_OUTPUT_HTML $DESIRED_OUTPUT_HTML
fi

if [ "$DESIRED_OUTPUT_ZIP" != "$REAL_OUTPUT_ZIP" ]
then
    mv $REAL_OUTPUT_ZIP $DESIRED_OUTPUT_ZIP
fi