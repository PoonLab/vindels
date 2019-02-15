from seqUtils import *
from glob import glob
import os
folder = glob ("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta_")

count = 0 
for infile in folder:
    with open(infile) as handle:
        fasta = parse_fasta(handle)
    
    os.rename(infile, infile[:-1])

    for seq in fasta:
        count += 1

print(count)

