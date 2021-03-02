from Bio import SeqIO
import re
import mappy as mp
from io import StringIO
import argparse
import os
import numpy as np

##* import van modin is later (na argparse sectie) zodat het aantal bruikbare threads ingesteld kan worden

arg = argparse.ArgumentParser()

arg.add_argument(
    "--input", "-i", metavar="File", help="Input FastQ file", type=str, required=True
)

arg.add_argument(
    "--reference",
    "-ref",
    metavar="File",
    help="Input reference fasta",
    type=str,
    required=True,
)

arg.add_argument(
    "--primers",
    "-pr",
    metavar="File",
    help="Fasta file with used primers (no ambiguity codes!)",
    type=str,
    required=True,
)

arg.add_argument(
    "--output", "-o", metavar="File", help="Output FastQ File", type=str, required=True
)

arg.add_argument(
    "--threads",
    "-t",
    metavar="N",
    help="Number of threads that can be used in parallel",
    type=int,
    default=4,
    required=False,
)

flags = arg.parse_args()

os.environ["MODIN_CPUS"] = str(flags.threads)
import modin.pandas as pd

### zoek de coordinaten van primersequenties
def search_primers(pattern, reference, id):
    for record in SeqIO.parse(reference, "fasta"):
        chrom = record.id

        for match in re.finditer(str(pattern), str(record.seq)):
            start_pos = match.start()
            end_pos = match.end()

            return chrom, start_pos, end_pos, id


"""
Maak een "BED" dataframe van de primercoordinaten
vereiste is een vorm van een orientatie in de naam van de primer (zoals in de fasta)
als voorbeeld:

>primer_2_LEFT
ACTGAGTATCG

De primers in de opgegeven fasta moeten altijd  5' --> 3' zijn.
Aan de hand van de orientatie opgegeven in de primernaam kunnen we terugvinden welke primers revcomp zijn
"""


def PrimerCoordinates(primers, ref):
    left = ["LEFT", "PLUS", "POSITIVE"]
    right = ["RIGHT", "MINUS", "NEGATIVE"]

    BedLeft = pd.DataFrame([])
    BedRight = pd.DataFrame([])

    for record in SeqIO.parse(primers, "fasta"):
        if any(orient in record.id for orient in left) is True:
            try:
                ref, start, end, id = search_primers(record.seq, reference, record.id)
                BedLeft = BedLeft.append(
                    pd.DataFrame(
                        {"chrom": ref, "start": start, "stop": end, "name": id},
                        index=[0],
                    ),
                    ignore_index=True,
                )
            except:
                ref, start, end, id = search_primers(
                    record.seq.reverse_complement(), reference, record.id
                )
                BedLeft = BedLeft.append(
                    pd.DataFrame(
                        {"chrom": ref, "start": start, "stop": end, "name": id},
                        index=[0],
                    ),
                    ignore_index=True,
                )
        if any(orient in record.id for orient in right) is True:
            try:
                ref, start, end, id = search_primers(record.seq, reference, record.id)
                BedRight = BedRight.append(
                    pd.DataFrame(
                        {"chrom": ref, "start": start, "stop": end, "name": id},
                        index=[0],
                    ),
                    ignore_index=True,
                )
            except:
                ref, start, end, id = search_primers(
                    record.seq.reverse_complement(), reference, record.id
                )
                BedRight = BedRight.append(
                    pd.DataFrame(
                        {"chrom": ref, "start": start, "stop": end, "name": id},
                        index=[0],
                    ),
                    ignore_index=True,
                )

    BedAll = pd.concat([BedLeft, BedRight])

    ForwardList = PrimerCoordinates_forward(BedLeft)
    ReverseList = PrimerCoordinates_reverse(BedRight)
    AllList = PrimerCoordinates_all(BedAll)

    ForwardList.sort()
    ReverseList.sort()
    AllList.sort()

    return ForwardList, ReverseList, AllList, BedLeft, BedRight


### Itereren over een dataframe per operatie duurt te lang.
## Maken dus een lijst van de coordinaten per primer-range en gebruiken dat
## nodig om de verschillende coordinaten binnen een range van de primers terug te vinden.
def PrimerCoordinates_forward(bed):
    coordlist = []
    for index, chrom in bed.iterrows():
        list = [*range(chrom.start, chrom.stop, 1)]
        for i in list:
            coordlist.append(i)

    return coordlist


def PrimerCoordinates_reverse(bed):
    coordlist = []
    for index, chrom in bed.iterrows():
        list = [*range(chrom.start, chrom.stop, 1)]
        for i in list:
            coordlist.append(i)

    return coordlist


def PrimerCoordinates_all(bed):
    coordlist = []
    for index, chrom in bed.iterrows():
        list = [*range(chrom.start, chrom.stop, 1)]
        for i in list:
            coordlist.append(i)

    return coordlist


def IndexReads(fastqfile):
    #### Big thanks to this guy who suggested using a giant dict instead of df.append for speed improvements
    #### https://stackoverflow.com/a/50105723

    ## Laad het FastQ bestand in memory zodat we geen operaties hoeven te doen op een filehandle
    ## vermijden van filesystem-latency
    with open(fastqfile, "r") as f:
        line = f.read()

    ReadDict = {}
    i = 0

    ## Biopython wil enkel en alleen handelen op filehandles (raar) dus is het nodig om via StringIO een filehandle te maken vanuit een memory-stream
    fastq_io = StringIO(line)
    for record in SeqIO.parse(fastq_io, "fastq"):

        RecordQualities = "".join(
            map(lambda x: chr(x + 33), record.letter_annotations["phred_quality"])
        )
        ReadDict[i] = {
            "Readname": str(record.id),
            "Sequence": str(record.seq),
            "Qualities": str(RecordQualities),
        }
        i = i + 1

    fastq_io.close()

    ## Maak van de dict een dataframe
    ReadIndex = pd.DataFrame.from_dict(ReadDict, "index")

    return ReadIndex


### de verschillende cutting-functies om de sequenties en bijbehorende qualities weg te halen
def slice_forward_left(readstart, seq, qual):
    readstart = readstart + 1
    trimmedseq = seq[1:]
    trimmedqual = qual[1:]

    return trimmedseq, trimmedqual, readstart


def slice_forward_right(readend, seq, qual):
    readend = readend - 1
    trimmedseq = seq[:-1]
    trimmedqual = qual[:-1]

    return trimmedseq, trimmedqual, readend


def slice_reverse_left(readstart, seq, qual):
    readstart = readstart + 1
    trimmedseq = seq[:-1]
    trimmedqual = qual[:-1]

    return trimmedseq, trimmedqual, readstart


def slice_reverse_right(readend, seq, qual):
    readend = readend - 1
    trimmedseq = seq[1:]
    trimmedqual = qual[1:]

    return trimmedseq, trimmedqual, readend


## basis van de alignerfunctie
def InitAligner(reference):
    AlignObject = mp.Aligner(reference, preset="map-ont")

    return AlignObject


def ReadBeforePrimer_FW(readstart, coordlist):
    diff = lambda val: abs(val - readstart)
    nearest_int = min(coordlist, key=diff)

    if readstart <= nearest_int:
        return True
    else:
        return False


def ReadAfterPrimer_FW(readend, coordlist):
    diff = lambda val: abs(val - readend)
    nearest_int = min(coordlist, key=diff)

    if readend >= nearest_int:
        return True
    else:
        return False


def ReadBeforePrimer_RV(readstart, coordlist):
    diff = lambda val: abs(val - readstart)
    nearest_int = min(coordlist, key=diff)

    if readstart <= nearest_int:
        return True
    else:
        return False


def ReadAfterPrimer_RV(readend, coordlist):
    diff = lambda val: abs(val - readend)
    nearest_int = min(coordlist, key=diff)

    if readend >= nearest_int:
        return True
    else:
        return False


def Cut_FastQ(input, reference, ForwardList, ReverseList, AllList):

    seq = input[1]
    qual = input[2]
    Aln = InitAligner(reference)

    for hit in Aln.map(seq):

        readseq = seq
        readqual = qual

        if hit.strand == 1:
            is_reverse = False
        if hit.strand == -1:
            is_reverse = True

        readstart = hit.r_st
        readend = hit.r_en

        if is_reverse is False:

            while ReadBeforePrimer_FW(readstart, ForwardList) is True:
                readseq, readqual, readstart = slice_forward_left(
                    readstart, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readstart = hit2.r_st

            while ReadAfterPrimer_FW(readend, ReverseList) is True:
                readseq, readqual, readend = slice_forward_right(
                    readend, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readend = hit2.r_en

            while (readstart in AllList) is True:
                readseq, readqual, readstart = slice_forward_left(
                    readstart, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readstart = hit2.r_st

            while (readend in AllList) is True:
                readseq, readqual, readend = slice_forward_right(
                    readend, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readend = hit2.r_en

        if is_reverse is True:

            while ReadBeforePrimer_RV(readstart, ForwardList) is True:
                readseq, readqual, readstart = slice_reverse_left(
                    readstart, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readstart = hit2.r_st

            while ReadAfterPrimer_RV(readend, ReverseList) is True:
                readseq, readqual, readend = slice_reverse_right(
                    readend, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readend = hit2.r_en

            while (readstart in AllList) is True:
                readseq, readqual, readstart = slice_reverse_left(
                    readstart, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readstart = hit2.r_st

            while (readend in AllList) is True:
                readseq, readqual, readend = slice_reverse_right(
                    readend, readseq, readqual
                )
                for hit2 in Aln.map(readseq):
                    readend = hit2.r_en
        
        if len(readseq) == 0:
            return np.nan
        
        return readseq, readqual


if __name__ == "__main__":
    reference = flags.reference
    primerfasta = flags.primers
    fastqfile = flags.input
    output = flags.output

    ForwardList, ReverseList, AllList, BedLeft, BedRight = PrimerCoordinates(
        primerfasta, reference
    )

    ReadFrame = IndexReads(fastqfile)

    ReadFrame["ProcessedReads"] = ReadFrame.apply(
        Cut_FastQ, args=(reference, ForwardList, ReverseList, AllList), axis=1
    )

    ReadFrame.dropna(subset=['ProcessedReads'], inplace=True)
    
    ReadFrame[["ProcessedSeq", "ProcessedQual"]] = pd.DataFrame(
        ReadFrame["ProcessedReads"].tolist(), index=ReadFrame.index
    )
    ReadFrame.drop(columns=["ProcessedReads", "Sequence", "Qualities"], inplace=True)

    ## Schrijf de resultaten naar de output fastq
    ## het gebruiken van dataframe.iterrows() is te langzaam, zeker als er veel reads zijn (lees: +3 uur voor het schrijven)
    ## daarom dezelfde truc voor het inladen van data andersom gebruiken
    ## van dataframe naar een lijst met dicts die volledig in-memory is.
    ## vanaf daar loopen over de lijst met dicts en schrijven naar een fastq file
    ReadDict = ReadFrame.to_dict(orient="records")

    with open(output, "w") as fileout:
        for index in range(len(ReadDict)):
            for key in ReadDict[index]:
                if key == "Readname":
                    fileout.write("@" + ReadDict[index][key] + "\n")
                if key == "ProcessedSeq":
                    fileout.write(str(ReadDict[index][key]) + "\n" + "+" + "\n")
                if key == "ProcessedQual":
                    fileout.write(str(ReadDict[index][key]) + "\n")
