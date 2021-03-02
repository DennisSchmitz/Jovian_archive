from Bio import SeqIO
import re
import argparse
import mappy as mp
from io import StringIO
import pandas as pd
import numpy as np
import multiprocessing


arg = argparse.ArgumentParser()

arg.add_argument(
    '--input','-i',
    metavar="File",
    help="Input Fastq file",
    type=str,
    required=True
)

arg.add_argument(
    '--reference', '-ref',
    metavar="File",
    help="Input reference fasta",
    type=str,
    required=True
)

arg.add_argument(
    '--primers', '-pr',
    metavar="File",
    help="Fasta file with used primers (no ambiguity codes!)",
    type=str,
    required=True
)

arg.add_argument(
    '--output', '-o',
    metavar="File",
    help="Output FastQ File",
    type=str,
    required=True
)

arg.add_argument(
    "--threads", "-t",
    metavar="N",
    help="Number of threads that can be used in parallel",
    type=int,
    default=2,
    required=False
)

flags = arg.parse_args()

def IndexReads(fastqfile):
    #### Big thanks to this guy who suggested using a giant dict instead of df.append for speed improvements
    #### https://stackoverflow.com/a/50105723
    
    ## Laad het FastQ bestand in memory zodat we geen operaties hoeven te doen op een filehandle
    ## vermijden van filesystem-latency
    with open(fastqfile, "r") as f:
        line = f.read()
    
    ReadDict = {}
    i = 0
    
    ## Biopython wil enkel en alleen handelen op filehandles (beetje raar) dus is het nodig om via StringIO een filehandle te maken vanuit een memory-stream
    fastq_io = StringIO(line)
    for record in SeqIO.parse(fastq_io, "fastq"):

        RecordQualities = ''.join(map(lambda x: chr( x+33 ), record.letter_annotations["phred_quality"]))
        seq = record.seq
        if len(seq) == 0:
            seq = np.nan
        ReadDict[i] = {"Readname": str(record.id), "Sequence": str(seq), 'Qualities': str(RecordQualities)}
        i = i + 1

    fastq_io.close()

    ## Maak van de dict een dataframe
    ReadIndex = pd.DataFrame.from_dict(ReadDict, "index")
    
    return ReadIndex

### zoek de coordinaten van primersequenties
def search_primers(pattern, reference, id):
    for record in SeqIO.parse(reference, "fasta"):
        chrom = record.id
        
        for match in re.finditer(str(pattern), str(record.seq)):
            start_pos = match.start()
            end_pos = match.end()
            
            return chrom, start_pos, end_pos, id
        
        else:
            
            return None, None, None, None

def PrimerCoordinates(primers, ref):
    left = ['LEFT', 'PLUS', "POSITIVE"]
    right = ['RIGHT', 'MINUS', 'NEGATIVE']
    
    BedLeft = pd.DataFrame([])
    BedRight = pd.DataFrame([])
    
    for record in SeqIO.parse(primers, "fasta"):
        if any(orient in record.id for orient in left) is True:
            ref, start, end, id = search_primers(record.seq, reference, record.id)
            BedLeft = BedLeft.append(pd.DataFrame({"chrom": ref, "start": start, "stop": end, "name": id}, index=[0]), ignore_index=True)
        if any(orient in record.id for orient in right) is True:
            ref, start, end, id = search_primers(record.seq, reference, record.id)
            BedRight = BedRight.append(pd.DataFrame({"chrom": ref, "start": start, "stop": end, "name": id}, index=[0]), ignore_index=True)
    
    ForwardList = PrimerCoordinates_forward(BedLeft)
    ReverseList = PrimerCoordinates_reverse(BedRight)
    
    return ForwardList, ReverseList, BedLeft, BedRight



def PrimerCoordinates_forward(bed):
    coordlist = []
    for index, chrom in bed.iterrows():
        list = [*range(chrom.start-1, chrom.stop, 1)]
        for i in list:
            coordlist.append(i)
            
    return coordlist

def PrimerCoordinates_reverse(bed):
    coordlist = []
    for index, chrom in bed.iterrows():
        list = [*range(chrom.start+1, chrom.stop+1, 1)]
        for i in list:
            coordlist.append(i)
            
    return coordlist

def slice_forward(readstart, seq, qual):
    readstart = readstart + 1
    trimmedseq = seq[1:]
    trimmedqual = qual[1:]

    return trimmedseq, trimmedqual, readstart

def slice_reverse(readend, seq, qual):
    readend = readend - 1
    trimmedseq = seq[1:]
    trimmedqual = qual[1:]

    return trimmedseq, trimmedqual, readend

def ReadBeforePrimer_FW(readstart, coordlist):
    if readstart in coordlist:
        return False
    diff = lambda val: abs(val - readstart)
    nearest_int = min(coordlist, key=diff)

    if readstart <= nearest_int:
        return True
    else:
        return False


def ReadAfterPrimer_RV(readend, coordlist):
    if readend in coordlist:
        return False
    diff = lambda val: abs(val - readend)
    nearest_int = min(coordlist, key=diff)

    if readend >= nearest_int:
        return True
    else:
        return False
    
    
def Cut_reads(Frame):

    reference = flags.reference
    primerfasta = flags.primers
    fwlist, rvlist, BedLeft, BedRight = PrimerCoordinates(primerfasta, reference)
    
    # init the aligner on a "per block" basis to reduce overhead.
    Aln = mp.Aligner(reference, preset="sr", best_n=1)
    
    readnames = Frame['Readname'].tolist()
    sequences = Frame['Sequence'].tolist()
    qualities = Frame['Qualities'].tolist()
    
    processed_readnames = []
    processed_sequences = []
    processed_qualities = []

    for i in range(len(readnames)):
        name = readnames[i]
        seq = sequences[i]
        qual = qualities[i]

        looplimiter = 0
        for hit in Aln.map(seq):
            if looplimiter != 0:
                continue
            looplimiter +=1

            if hit.strand == 1:
                reverse = False
            if hit.strand == -1:
                reverse = True

            start = hit.r_st
            end = hit.r_en
            
            qstart = hit.q_st

            if reverse is False:
                
                if qstart != 0:
                    seq = seq[qstart:]
                    qual = qual[qstart:]
                    
                
                while ReadBeforePrimer_FW(start, fwlist) is True:
                    seq, qual, start = slice_forward(start, seq, qual)
                    hitlimiter = 0
                    for hit2 in Aln.map(seq):
                        if hitlimiter != 0:
                            continue
                        hitlimiter +=1
                        start = hit2.r_st

                stepper = 0
                while (start in fwlist) is True:
                    seq, qual, start = slice_forward(start, seq, qual)
                    if stepper == 3:
                        hitlimiter = 0
                        for hit2 in Aln.map(seq):
                            if hitlimiter != 0:
                                continue
                            hitlimiter +=1
                            start = hit2.r_st
                        stepper = 0
                        continue
                    stepper +=1

            if reverse is True:
                
                while ReadAfterPrimer_RV(end, rvlist) is True:
                    seq, qual, end = slice_reverse(end, seq, qual)
                    hitlimiter = 0
                    for hit2 in Aln.map(seq):
                        if hitlimiter != 0:
                            continue
                        hitlimiter +=1
                        end = hit2.r_en

                stepper = 0
                while (end in rvlist) is True:
                    seq, qual, end = slice_reverse(end, seq, qual)
                    if stepper == 3:
                        hitlimiter = 0
                        for hit2 in Aln.map(seq):
                            if hitlimiter != 0:
                                continue
                            hitlimiter +=1
                            end = hit2.r_en
                        stepper = 0
                        continue
                    stepper +=1

        if len(seq) == 0:
            seq = np.nan
        if len(qual) == 0:
            qual = np.nan

        processed_readnames.append(name)
        processed_sequences.append(seq)
        processed_qualities.append(qual)
        
    ProcessedReads = pd.DataFrame({
        "Readname": processed_readnames,
        "Sequence": processed_sequences,
        "Qualities": processed_qualities})
    return ProcessedReads

def parallel(df, func, workers=min(multiprocessing.cpu_count(), 128)):
    df_split = np.array_split(df, workers)
    pool = multiprocessing.Pool(workers)
    df = pd.concat(pool.map(func, df_split))
    pool.close()
    pool.join()
    return df

if __name__ == "__main__":
    reference = flags.reference
    primerfasta = flags.primers
    fastqfile = flags.input
    output = flags.output
    
    
    ReadFrame = IndexReads(fastqfile)
    
    ReadFrame.dropna(subset=['Sequence'], inplace=True)
    ReadFrame = ReadFrame.sample(frac=1).reset_index(drop=True)
    
    ProcessedReads = parallel(ReadFrame, Cut_reads, flags.threads)
    
    ProcessedReads.dropna(subset=['Sequence', 'Qualities'], inplace=True)
    
    ReadDict = ProcessedReads.to_dict(orient='records')
    
    with open(output, "w") as fileout:
        for index in range(len(ReadDict)):
            for key in ReadDict[index]:
                if key == "Readname":
                    fileout.write(
                        "@" + ReadDict[index][key] + "\n"
                    )
                if key == 'Sequence':
                    fileout.write(
                        str(ReadDict[index][key]) + "\n" + "+" + "\n"
                    )
                if key == 'Qualities':
                    fileout.write(
                        str(ReadDict[index][key]) + "\n"
                    )