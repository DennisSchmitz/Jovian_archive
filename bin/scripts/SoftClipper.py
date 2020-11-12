import pysam
import re
import argparse

arg = argparse.ArgumentParser()

arg.add_argument(
    "--input",
    metavar="File",
    help="BAM file of the first minimap2 alignment run, sorted and indexed.",
    type=str,
    required=True,
)

arg.add_argument(
    "--output",
    metavar="File",
    help="File with cleaned fastq reads",
    type=argparse.FileType("w"),
    required=True,
)

arg.add_argument(
    "--threads",
    metavar="Number",
    help="Number of threads that can be used for decompressing/compressing the BAM file",
    default=1,
    type=int,
    required=False,
)

flags = arg.parse_args()

with flags.output as fileout:
    bamfile = pysam.AlignmentFile(flags.input, "rb", threads=flags.threads)
    for read in bamfile:
        try:
            s_clipped_left = re.match("^(\d*S)", read.cigarstring).group().split("S")[0]
        except:
            s_clipped_left = "0"

        try:
            s_clipped_right = (
                re.search("(\d*S)$", read.cigarstring).group().split("S")[0]
            )
        except:
            s_clipped_right = "0"

        trim_seq_left = read.seq[int(s_clipped_left) :]
        trim_seq_right = trim_seq_left[: -int(s_clipped_right)]

        trimmed_seq = trim_seq_right

        trim_qual_left = read.qual[int(s_clipped_left) :]
        trim_qual_right = trim_qual_left[: -int(s_clipped_right)]

        trimmed_qual = trim_qual_right
        
        if read.is_reverse == True:
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
            bases = list(trimmed_seq)
            bases = [complement[base] for base in bases]
            trimmed_seq = ''.join(bases)
            trimmed_seq = trimmed_seq[::-1]
            trimmed_qual = trimmed_qual[::-1]

        fileout.write(
            "@"
            + str(read.query_name)
            + "\n"
            + str(trimmed_seq)
            + "\n"
            + "+"
            + "\n"
            + str(trimmed_qual)
            + "\n"
        )

fileout.close()
