from Bio import SeqIO
import argparse

arg = argparse.ArgumentParser()

arg.add_argument(
    '--primers',
    metavar="File",
    type=str,
    required=True
)

arg.add_argument(
    '--three',
    metavar="File",
    type=argparse.FileType("w"),
    required=True
)

arg.add_argument(
    '--five',
    metavar="File",
    type=argparse.FileType("w"),
    required=True
)

flags = arg.parse_args()

def three_prime(primerfile, three_end_out):
    for record in SeqIO.parse(primerfile, "fasta"):
        three_end_out.write(f">{record.id}" + "\n" + str(record.seq) + "X" + "\n")

def five_prime(primerfile, five_end_out):
    for record in SeqIO.parse(primerfile, "fasta"):
        five_end_out.write(f">{record.id}" + "\n" + "X" + str(record.seq) + "\n")

if __name__ == "__main__":
    filein = flags.primers
    
    three_prime(filein, flags.three)
    five_prime(filein, flags.five)