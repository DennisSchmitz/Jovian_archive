from Bio import SeqIO
import re
import pysam
import argparse
import os

##* import van modin is later (na argparse sectie) zodat het aantal bruikbare threads ingesteld kan worden

arg = argparse.ArgumentParser()

arg.add_argument(
    '--input','-i',
    metavar="File",
    help="Input BAM file",
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
        list = [*range(chrom.start-1, chrom.stop+1, 1)]
        for i in list:
            coordlist.append(i)
            
    return coordlist

def PrimerCoordinates_reverse(bed):
    coordlist = []
    for index, chrom in bed.iterrows():
        list = [*range(chrom.start, chrom.stop+2, 1)]
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
    trimmedseq = seq[:-1]
    trimmedqual = qual[:-1]

    return trimmedseq, trimmedqual, readend

def IndexReads(bamfile):
    ReadDict = {}
    i = 0
    
    for read in bamfile:
        
        if read.is_unmapped is True:
            continue
        
        readname = read.query_name
        IsReverse = read.is_reverse
        ReadStart = read.reference_start + 1
        ReadEnd = read.reference_end + 1
        ReadSeq = read.query_sequence
        ReadQual = quality = ''.join(map(lambda x: chr( x+33 ), read.query_qualities))
        
        ReadDict[i] = {
            "Readname": str(readname), 
            "Sequence": str(ReadSeq), 
            "Qualities": str(ReadQual), 
            "StartPos": int(ReadStart), 
            "EndPos": int(ReadEnd), 
            "Reverse": bool(IsReverse)
            }
        
        i = i + 1
        
    ReadIndex = pd.DataFrame.from_dict(ReadDict, "index")
    
    return ReadIndex

def FlipStrand(seq, qual):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    bases = list(seq)
    bases = [complement[base] for base in bases]
    seq = "".join(bases)
    seq = seq[::-1]
    qual = qual[::-1]
    
    return seq, qual

def Cut_reads(input, bedLeft, bedRight, fwlist, rvlist):
    seq = input[1]
    qual = input[2]
    start = input[3]
    end = input[4]
    reverse = input[5]
    
    start_of_first_primer = bedLeft.start.iloc[0] - 2
    end_of_last_primer = bedRight.stop.iloc[-1] + 2

    if reverse is False:
        if start < start_of_first_primer:
            to_cut = start_of_first_primer - start
            start = start + to_cut
            seq = seq[to_cut:]
            qual = qual[to_cut:]
        
        while (start in fwlist) is True:
            seq, qual, start = slice_forward(start, seq, qual)

    if reverse is True:
        if end > end_of_last_primer:
            to_cut = end + end_of_last_primer
            end = end - to_cut
            seq = seq[:-to_cut]
            qual = qual[:-to_cut]
            
        while (end in rvlist) is True:
            seq, qual, end = slice_reverse(end, seq, qual)
        
        seq, qual = FlipStrand(seq, qual)
        
    return seq, qual

if __name__ == "__main__":
    reference = flags.reference
    primerfasta = flags.primers
    bamfile = flags.input
    output = flags.output
    
    
    ForwardList, ReverseList, BedLeft, BedRight = PrimerCoordinates(primerfasta, reference)
    
    
    bam = pysam.AlignmentFile(bamfile, "rb")
    
    ReadFrame = IndexReads(bam)
    
    ReadFrame["ProcessedReads"] = ReadFrame.apply(Cut_reads, args=(BedLeft, BedRight, ForwardList, ReverseList), axis=1)
    ReadFrame[["ProcessedSeq", "ProcessedQual"]] = pd.DataFrame(ReadFrame["ProcessedReads"].tolist(), index=ReadFrame.index)
    
    ReadFrame.drop(columns=['ProcessedReads', 'Sequence', 'Qualities', "StartPos", "EndPos", "Reverse"], inplace=True)
    
    ReadDict = ReadFrame.to_dict(orient='records')
    
    with open(output, "w") as fileout:
        for index in range(len(ReadDict)):
            for key in ReadDict[index]:
                if key == "Readname":
                    fileout.write(
                        "@" + ReadDict[index][key] + "\n"
                    )
                if key == 'ProcessedSeq':
                    fileout.write(
                        str(ReadDict[index][key]) + "\n" + "+" + "\n"
                    )
                if key == 'ProcessedQual':
                    fileout.write(
                        str(ReadDict[index][key]) + "\n"
                    )