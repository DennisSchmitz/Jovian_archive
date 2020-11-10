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
        max_depth=0,
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

        exists_in_orf = []
        for index, orf in IndexedGFF.iterrows():
            in_orf = currentloc in range(orf.start, orf.end)
            exists_in_orf.append(in_orf)

        if any(exists_in_orf) == True:
            ORFPosition = True
        elif any(exists_in_orf) == False:
            ORFPosition = False

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
        cur_nuc_dist = {
            "A": slice_c[1],
            "T": slice_c[2],
            "C": slice_c[3],
            "G": slice_c[4],
            "D": slice_c[5],
        }
        prv_nuc_dist = {
            "A": slice_p1[1],
            "T": slice_p1[2],
            "C": slice_p1[3],
            "G": slice_p1[4],
            "D": slice_p1[5],
        }
        prv2_nuc_dist = {
            "A": slice_p2[1],
            "T": slice_p2[2],
            "C": slice_p2[3],
            "G": slice_p2[4],
            "D": slice_p2[5],
        }
        nxt_nuc_dist = {
            "A": slice_n1[1],
            "T": slice_n1[2],
            "C": slice_n1[3],
            "G": slice_n1[4],
            "D": slice_n1[5],
        }
        nxt2_nuc_dist = {
            "A": slice_n2[1],
            "T": slice_n2[2],
            "C": slice_n2[3],
            "G": slice_n2[4],
            "D": slice_n2[5],
        }

        # get the most primary nucleotide and secondary nucleotide for every position
        ## >> sort the distribution of the earlier made dict based on the values, return the keys with the highest and secondary highest values
        # > current pos
        cur_sorted_dist = sorted(((value, key) for key, value in cur_nuc_dist.items()))
        cur_primary_nuc = cur_sorted_dist[-1][1]
        cur_second_nuc = cur_sorted_dist[-2][1]

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

        # als de coverage op de "currentposition" lager is dan de minimale coverage, plak dan een "N"
        # Zoniet, ga door met de daadwerkelijke nucleotides checken
        if cur_cov < mincov:
            standard_cons.append("N")
            corrected_cons.append("N")
        else:
            # als de "currentposition" géén deletie is, plak dan de meerderheid (A/T/C/G/D) in de index van deze positie (komt uit alignment)
            if cur_primary_nuc != "D":
                standard_cons.append(cur_primary_nuc)
                corrected_cons.append(cur_primary_nuc)
            # als de "currentposition" wél een deletie is, ga dan kijken naar de status van omliggende nucleotides
            elif cur_primary_nuc == "D":
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
                    if nxt_primary_nuc != "D" and prv_primary_nuc != "D":
                        is_del = False

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "-1" en "+1" hebben beide wél een deletie gerapporteerd
                    # keur dan de gerapporteerde deletie goed
                    if nxt_primary_nuc == "D" and prv_primary_nuc == "D":
                        is_del = True

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "+1" en "+2" hebben beide wél een deletie gerapporteerd
                    # keur dan de gerapporteerde deletie goed
                    if nxt_primary_nuc == "D" and nxt2_primary_nuc == "D":
                        is_del = True

                    # In het geval er een deletie is gerapporteerd als meerderheid op de "currentposition"
                    # en zowel nucleotide positie "-1" en "-2" hebben beide wél een deletie gerapporteerd
                    # keur dan de gerapporteerde deletie goed
                    if prv_primary_nuc == "D" and prv2_primary_nuc == "D":
                        is_del = True

                    if is_del == False:
                        corrected_cons.append(cur_second_nuc)
                    elif is_del == True:
                        corrected_cons.append("-")

    sequences = "".join(standard_cons) + "," + "".join(corrected_cons)
    return sequences


if __name__ == "__main__":
    GFF_index = MakeGFFindex(flags.gff)
    pileindex = BuildIndex(flags.input, flags.reference)
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
