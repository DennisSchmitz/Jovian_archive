#! /bin/bash

# Script to count reads mapped to contigs,
# counting ONLY PRIMARY alignments and 
# NO SUPPLEMENTARY alignments.

# Returns counts as a table with two columns:
# 'scaffold_name' and 'mapped_reads'
# (for merging with all_taxClassified.tsv and
# all_taxUnclassified.tsv in bin/quantify_output.py).

INPUT_FOLDER=$1
BAM_FILES="${INPUT_FOLDER}/*_sorted.bam"
RESULT_FILE=$2

printf "mapped_reads\tscaffold_name\n" #> $RESULT_FILE
# have the complete script redirect ouput to a file (and log file?)

for file in $BAM_FILES
do
    paired_mapped_reads=$(samtools view -F 4 -F 8 -F 256 -F 2048 $file | cut -f 3 | sort | uniq -c | sed $"s/NODE_/\tNODE_/g")
    #cut the contig names, sort and count, and place tabs between numbers and contigs names
    unpaired_mapped_reads=$(samtools view -F 4 -f 8 -F 256 -F 2048 $file | cut -f 3 | sort | uniq -c | sed $"s/NODE_/\tNODE_/g")
    mapped_reads=$(samtools view -F 4 -F 256 -F 2048 $file | cut -f 3 | sort | uniq -c | sed $"s/NODE_/\tNODE_/g")
    #mapped_reads = both paired and unpaired! (do not filter for paired or not = 8)
    samtools view -F 4 -F 256 -F 2048 $file | cut -f 3 | sort | uniq -c | sed $"s/NODE_/\tNODE_/g"
done