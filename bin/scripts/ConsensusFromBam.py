#!/usr/bin/env python


# @ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# @ This script is a modified version of the "bam2consensus.py" script
# @ written by David Nieuwenhuijse (https://github.com/dnieuw)
# @ The "logic" is unchanged from the original script.
# @
# @ The original script can be found at https://github.com/dnieuw/ENA_SARS_Cov2_nanopore/blob/master/bin/bam2consensus.py
# @
# @ # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

from argparse import FileType
import argparse
from typing import Counter
import pysam


par = argparse.ArgumentParser()

par.add_argument(
    "-input", metavar="File", help="Sorted + Indexed BAM file", type=str, required=True
)

par.add_argument("-output", metavar="File", type=argparse.FileType("w"), required=True)

par.add_argument(
    "-cov", help="Minimal coverage threshold", default=1, type=int, required=True
)

par.add_argument(
    "-name", help="Sequence name to write to the fasta file", type=str, required=True
)

flags = par.parse_args()


def buildCons(BamAln, mincov):
    with pysam.AlignmentFile(BamAln, "rb") as bam:
        all_references = []
        for ref in bam.references:
            cons = []
            prev = -1
            for n, pileupcolumn in enumerate(
                bam.pileup(
                    ref, ignore_orphans=False, min_mapping_quality=0, min_base_quality=0
                )
            ):
                if pileupcolumn.pos != prev + 1:
                    gapregion = pileupcolumn.pos - prev - 1
                    cons.append("N" * gapregion)
                prev = pileupcolumn.pos

                coverage = pileupcolumn.get_num_aligned()
                if coverage < mincov:
                    cons.append("N")
                    continue

                pos = Counter()
                for read in pileupcolumn.pileups:
                    if read.is_del:
                        pos["*"] += 1
                    elif not read.is_refskip:
                        nuc = read.alignment.query_sequence[read.query_position]
                        pos[nuc] += 1

                match = max(pos, key=pos.get)

                if match == "*":
                    cons.append("-")
                    continue

                cons.append(match)

            if bam.get_reference_length(ref) != prev:
                gapregion = bam.get_reference_length(ref) - prev - 1
                cons.append("N" * gapregion)

            if len(cons) == 0:
                continue
            if flags.name != None:
                all_references.append(">" + flags.name + "\n" + "".join(cons))

        return "\n".join(all_references)


if __name__ == "__main__":
    cons = buildCons(flags.input, flags.cov)
    with flags.output as consfile:
        print(cons, file=consfile)
