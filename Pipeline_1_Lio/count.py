from seqUtils import *
from glob import glob
folder = glob ("/home/jpalmer/PycharmProjects/hiv-evolution-master/1_SubtypeSequences/*.fasta")

count = 0 
for infile in folder:
    with open(infile) as handle:
        fasta = parse_fasta(handle)
    

    for seq in fasta:
        count += 1

print(count)

