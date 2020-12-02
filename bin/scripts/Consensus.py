from os import system
import pysam
import pysamstats
import pandas as pd
import argparse
import gffpandas.gffpandas as gffpd

arg = argparse.ArgumentParser()

arg.add_argument(
    "--input",
    "-i",
    metavar="File",
    help="BAM file, sorted and indexed.",
    type=str,
    required=True,
)

arg.add_argument(
    "--reference",
    "-ref",
    metavar="File",
    help="Reference fasta file",
    type=str,
    required=True,
)

arg.add_argument(
    "--mincov",
    "-mc",
    metavar="Number",
    help="Minimum coverage threshold",
    default=1,
    type=int,
    required=True,
)

arg.add_argument(
    "--name",
    help="Name that will be used to make the fasta files",
    type=str,
    required=True,
)

arg.add_argument(
    "--consensus",
    metavar="File",
    help="File with the non-corrected consensus sequence",
    type=argparse.FileType("w"),
    required=True,
)

arg.add_argument(
    "--gapcorrected",
    metavar="File",
    help="File with the corrected consensus sequence",
    type=argparse.FileType("w"),
    required=True,
)

arg.add_argument(
    "--gff",
    metavar="File",
    help="GFF3 file produced by Prodigal",
    type=str,
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

arg.add_argument(
    "--coverage",
    "-cov",
    metavar="File",
    help="Output file listing the coverage per position",
    type=argparse.FileType("w"),
    required=True,
)

flags = arg.parse_args()


def MakeGFFindex(gff3file):
    GFFindex = gffpd.read_gff3(gff3file)

    return GFFindex.df


def BuildIndex(aln, fasta):
    bam = pysam.AlignmentFile(aln, "rb", threads=flags.threads)

    columns = ["coverage", "A", "T", "C", "G", "D"]
    pileup_index = pd.DataFrame(columns=columns)

    for rec in pysamstats.stat_pileup(
        type="variation",
        alignmentfile=bam,
        fafile=fasta,
        pad=True,
        one_based=True,
        max_depth=1000000000,
    ):
        pileup_index.loc[rec["pos"]] = (
            [rec["reads_all"]]
            + [rec["A"]]
            + [rec["T"]]
            + [rec["C"]]
            + [rec["G"]]
            + [rec["deletions"]]
        )

    return pileup_index


def ORFfinder(location, GFFindex):
    exists_in_orf = []
    for index, orf in GFFindex.iterrows():
        in_orf = location in range(orf.start, orf.end)
        exists_in_orf.append(in_orf)

    if any(exists_in_orf) == True:
        in_orf = True
    elif any(exists_in_orf) == False:
        in_orf = False
    else:
        in_orf = False
        system.exit(
            """
                    Something went wrong during the processing of this GFF.
                    Please check or re-generate the GFF file and try again.
                    """
        )

    return in_orf


def slices(mintwo, minone, zero, plusone, plustwo):
    dist_mintwo = {
        "A": mintwo[1],
        "T": mintwo[2],
        "C": mintwo[3],
        "G": mintwo[4],
        "D": mintwo[5],
    }
    dist_minone = {
        "A": minone[1],
        "T": minone[2],
        "C": minone[3],
        "G": minone[4],
        "D": minone[5],
    }
    dist_zero = {
        "A": zero[1],
        "T": zero[2],
        "C": zero[3],
        "G": zero[4],
        "D": zero[5],
    }
    dist_plusone = {
        "A": plusone[1],
        "T": plusone[2],
        "C": plusone[3],
        "G": plusone[4],
        "D": plusone[5],
    }
    distplustwo = {
        "A": plustwo[1],
        "T": plustwo[2],
        "C": plustwo[3],
        "G": plustwo[4],
        "D": plustwo[5],
    }

    return dist_mintwo, dist_minone, dist_zero, dist_plusone, distplustwo


def BuildCoverage(pileupindex):
    with flags.coverage as coverage_output:
        for index, rows in pileupindex.iterrows():
            slice_c = []
            for items in pileupindex.loc[index]:
                slice_c.append(items)
            coverage = slice_c[0]
            coverage_output.write(
                str(index)
                + "\t"
                + str(coverage)
                + "\n"
            )

def BuildCons(pileupindex, IndexedGFF, mincov):
    standard_cons = []
    corrected_cons = []

    for index, rows in pileupindex.iterrows():

        currentloc = index
        secondback = currentloc - 2
        firstback = currentloc - 1
        firstnext = currentloc + 1
        secondnext = currentloc + 2

        lastposition = pileupindex.tail(1).index.item()

        if currentloc == 1:
            firstback = 1
            secondback = 1
        if currentloc == 2:
            secondback = 1

        ORFPosition = ORFfinder(currentloc, IndexedGFF)

        # currentloc = slice_c
        # firstnext = slice_n1
        # secondnext = slice_n2
        # firstback = slice_p1
        # secondback = slice_np2

        slice_c = []
        slice_n1 = []
        slice_n2 = []
        slice_p1 = []
        slice_p2 = []

        for items in pileupindex.loc[currentloc]:
            slice_c.append(items)

        if firstnext < lastposition:
            for items in pileupindex.loc[firstnext]:
                slice_n1.append(items)
        else:
            for items in pileupindex.loc[lastposition]:
                slice_n1.append(items)

        if secondnext < lastposition:
            for items in pileupindex.loc[secondnext]:
                slice_n2.append(items)
        else:
            for items in pileupindex.loc[lastposition]:
                slice_n2.append(items)

        for items in pileupindex.loc[firstback]:
            slice_p1.append(items)
        for items in pileupindex.loc[secondback]:
            slice_p2.append(items)

        ## get the actual nucleotide distributions for every position

        prv2_nuc_dist, prv_nuc_dist, cur_nuc_dist, nxt_nuc_dist, nxt2_nuc_dist = slices(
            slice_p2, slice_p1, slice_c, slice_n1, slice_n2
        )


        # get the most primary nucleotide and secondary nucleotide for every position
        ## >> sort the distribution of the earlier made dict based on the values, return the keys with the highest and secondary highest values
        # > current pos
        cur_sorted_dist = sorted(((value, key) for key, value in cur_nuc_dist.items()))
        # Most abundant nuc at current pos: set it to lower case when that nuc is < mincov, else, uppercase.
        # (locations are filtered based on cur_cov (i.e. sum count of all nucs, incl dels), however, when the most
        # abundant nuc is < mincov (but sum total count of all nucs still >= mincov) it makes it a lowercase letter to
        # signify uncertainty.
        if cur_sorted_dist[-1][0] < mincov:
            cur_primary_nuc = cur_sorted_dist[-1][1].lower()
        else:
            cur_primary_nuc = cur_sorted_dist[-1][1].upper()

        # 2nd most abundant nuc at current pos, used later for gap-filling when del is inframe:
        # set it to N if count of that nuc == 0, this to assure that random nucs aren't inserted when count is zero.
        # Set to lower case when that nuc is < mincov, see explanation above, else set uppercase.
        if cur_sorted_dist[-2][0] == 0:
            cur_second_nuc = "N"
        elif cur_sorted_dist[-2][0] < mincov:
            cur_second_nuc = cur_sorted_dist[-2][1].lower()
        else:
            cur_second_nuc = cur_sorted_dist[-2][1].upper()
        
        # 3rd most abundant nuc at current pos, used later for gap-filling when del is inframe:
        # set it to N if count of that nuc == 0, this to assure that random nucs aren't inserted when count is zero.
        # Set to lower case when that nuc is < mincov, see explanation above, else set uppercase.
        if cur_sorted_dist[-3][0] == 0:
            cur_third_nuc = "N"
        elif cur_sorted_dist[-3][0] < mincov:
            cur_third_nuc = cur_sorted_dist[-3][1].lower()
        else:
            cur_third_nuc = cur_sorted_dist[-3][1].upper()

        # > next pos
        nxt_sorted_dist = sorted(((value, key) for key, value in nxt_nuc_dist.items()))
        nxt_primary_nuc = nxt_sorted_dist[-1][1]

        # > nextnext pos
        nxt2_sorted_dist = sorted(
            ((value, key) for key, value in nxt2_nuc_dist.items())
        )
        nxt2_primary_nuc = nxt2_sorted_dist[-1][1]

        # > prev pos
        prv_sorted_dist = sorted(((value, key) for key, value in prv_nuc_dist.items()))
        prv_primary_nuc = prv_sorted_dist[-1][1]

        # > prevprev pos
        prv2_sorted_dist = sorted(
            ((value, key) for key, value in prv2_nuc_dist.items())
        )
        prv2_primary_nuc = prv2_sorted_dist[-1][1]

        # Get the current coverage
        cur_cov = slice_c[0]

        # if 2nd most abundant nuc is not "N" (see above, i.e. its non-zero) and most and 2nd most abundant nucs are tied in counts, and DoC at cur-pos is >= mincov, then its an ambigious call and warning is thrown to stdout.
        # TODO exclude the "D"-nuc from the comparison, not informative and almost always artefacts IMHO, but we can discuss it later. Currently, warnings are thrown that an equal "D" count is ambigious
        if (
            (cur_second_nuc != "N")
            and (cur_sorted_dist[-1][0] == cur_sorted_dist[-2][0])
            and (cur_cov >= mincov)
        ):
            # TODO hier later nog een optie van maken om abiguity nuc-code te outputten op user request
            print(
                "File: ",
                flags.input,
                ". Ambigious call during consensus-calling at position ",
                currentloc,
                ', a "',
                cur_primary_nuc,
                '" was called (N=',
                cur_sorted_dist[-1][0],
                ') but could also be a "',
                cur_second_nuc,
                '" with (N=',
                cur_sorted_dist[-2][0],
                ").",
                sep="",
            )

        # als de coverage op de "currentposition" lager is dan de minimale coverage, plak dan een "N"
        # Zoniet, ga door met de daadwerkelijke nucleotides checken
        if cur_cov < mincov:
            standard_cons.append("N")
            corrected_cons.append("N")
        else:
            # als de "currentposition" géén deletie is, plak dan de meerderheid (A/T/C/G/D) in de index van deze positie (komt uit alignment)
            if cur_primary_nuc.upper() != "D":
                standard_cons.append(cur_primary_nuc)
                corrected_cons.append(cur_primary_nuc)
            # als de "currentposition" wél een deletie is, ga dan kijken naar de status van omliggende nucleotides
            elif cur_primary_nuc.upper() == "D":
                standard_cons.append(
                    "-"
                )  ## <-- deze is hier om beide een "standaard consensus" te maken naast de "gap-corrected consensus"
                if ORFPosition == False:
                    corrected_cons.append("-")
                elif ORFPosition == True:

                    is_del = False

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "-1" en "+1" (t.o.v. de "currentposition", dit is positie 0) beide géén deletie hebben
                    # vul dan de positie met de tweede meest voorkomende gerapporteerde nucleotide (A/C/T/G/D)
                    if (
                        nxt_primary_nuc.upper() != "D"
                        and prv_primary_nuc.upper() != "D"
                    ):
                        is_del = False

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "-1" en "+1" hebben beide wél een deletie gerapporteerd
                    # keur dan de gerapporteerde deletie goed
                    if (
                        nxt_primary_nuc.upper() == "D"
                        and prv_primary_nuc.upper() == "D"
                    ):
                        is_del = True

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "+1" en "+2" hebben beide wél een deletie gerapporteerd
                    # keur dan de gerapporteerde deletie goed
                    if (
                        nxt_primary_nuc.upper() == "D"
                        and nxt2_primary_nuc.upper() == "D"
                    ):
                        is_del = True

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "-1" en "-2" hebben beide wél een deletie gerapporteerd
                    # keur dan de gerapporteerde deletie goed
                    if (
                        prv_primary_nuc.upper() == "D"
                        and prv2_primary_nuc.upper() == "D"
                    ):
                        is_del = True

                    if is_del == False:
                        # if cur_second_nuc and cur_third_nuc are not "N" (see above, i.e. its non-zero) and 2nd and 3rd most abundant nucs are tied in counts,
                        # its an ambigious call and throw a warning to std. out.
                        if (cur_second_nuc and cur_third_nuc != "N") and (
                            cur_sorted_dist[-2][0] == cur_sorted_dist[-3][0]
                        ):
                            print(
                                "File: ",
                                flags.input,
                                ". Ambigious call during gap-filling at position ",
                                currentloc,
                                ', a "',
                                cur_second_nuc,
                                '" was called (N=',
                                cur_sorted_dist[-2][0],
                                ') but could also be a "',
                                cur_third_nuc,
                                '" with (N=',
                                cur_sorted_dist[-3][0],
                                ").",
                                sep="",
                            )
                        corrected_cons.append(cur_second_nuc)
                    elif is_del == True:
                        corrected_cons.append("-")

    sequences = "".join(standard_cons) + "," + "".join(corrected_cons)
    return sequences


if __name__ == "__main__":
    GFF_index = MakeGFFindex(flags.gff)
    pileindex = BuildIndex(flags.input, flags.reference)
    BuildCoverage(pileindex)
    sequences = BuildCons(pileindex, GFF_index, flags.mincov)

    standard_seq = sequences.split(",")[0]
    corrected_seq = sequences.split(",")[1]
    with flags.consensus as raw_consensus_seq:
        raw_consensus_seq.write(
            ">"
            + flags.name
            + "_standard_consensus_cov_"
            + str(flags.mincov)
            + "\n"
            + standard_seq
            + "\n"
        )
    with flags.gapcorrected as corrected_consensus_seq:
        corrected_consensus_seq.write(
            ">"
            + flags.name
            + "_gap-corrected_consensus_cov_"
            + str(flags.mincov)
            + "\n"
            + corrected_seq
            + "\n"
        )
