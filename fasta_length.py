## simple python script to read fasta file and output gene table of fasta headers and lengths in bp
## works for nucleotide and amino acid fasta files
## requires python 2.7 will not run with 3+

import sys
from Bio import SeqIO

FastaFile = open(sys.argv[1], 'rU')

for rec in SeqIO.parse(FastaFile, 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    print name, seqLen

FastaFile.close()
