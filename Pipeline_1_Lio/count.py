from seqUtils import *
from glob import glob
import os 
folder = glob ("/home/jpalmer/PycharmProjects/hiv-evolution-master/4_1_Edited/*.fasta")

count = 0 
for infile in folder:
    filename = os.path.basename(infile)
    outfile = open("/home/jpalmer/dart/data/unaligned/"+filename, 'w')

    with open(infile) as handle:
        fasta = parse_fasta(handle)
    

    for header in fasta:
        sliced = header.split(".")[-1]
        outfile.write('>'+sliced+'\n'+fasta[header]+'\n')


