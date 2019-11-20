#! /bin/bash

# Script to count reads mapped to scaffolds,
# counting ONLY PRIMARY alignments and 
# NO SUPPLEMENTARY alignments.

# Returns counts as a table with two columns:
# 'scaffold_name' and 'mapped_reads'
# (for merging with all_taxClassified.tsv and
# all_taxUnclassified.tsv in bin/quantify_output.py).

# Requires an input file as command-line argument.
# Recommended to redirect output to a (tsv) file.
# Example:
# $ bash count_mapped_reads.sh data/mapped_reads.bam > results/mapped_read_counts.tsv

INPUT_FILE=$1

#Write a header line
printf "mapped_reads\tscaffold_name\n"

#Append read counts and scaffold names
samtools view -F 4 -F 256 -F 2048 $INPUT_FILE | cut -f 3 | sort | uniq -c | sed $"s/NODE_/\tNODE_/g"

#Pipe explained:
# 1. use samtools to parse the bam file (samtools view)
# 1.1 exclude unmapped reads (-F 4)
# 1.2 exclude 'not primary alignments' (-F 256)
# 1.3 exclude supplementary alignments (-F 2048)
# (mates may be extracted separately by using the flag -F 8/-f 8
#  to exclude or include reads for which the mate is unmapped)
# 2. cut only the scaffold names from the sam file output (cut)
# 3. count the number of occurrences of each scaffold name (sort | uniq -c)
# 4. insert a tab to make a tab-separated table (sed)
