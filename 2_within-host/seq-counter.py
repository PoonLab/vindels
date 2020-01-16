# sequence counter 
import os 
import sys
from glob import * 
import re
from seqUtils import *

files = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/mcc/*a_recon.fasta")
count = 0
seqs = 0
for f in files:
    filename = os.path.basename(f)
    #print(filename)
    infile = open(f, "rU")
    data = parse_fasta(infile)
    seqs += len(data)
    name = filename.split("-a")[0]
    print(name)

    count +=1 
print(count)
print(seqs)
