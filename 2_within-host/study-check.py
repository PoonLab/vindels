import os 
import sys 
from seqUtils import * 
from glob import glob 
#print(glob.__file__)

infolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta")

toCheck = ["KC247412", 'KC247375']

for infile in infolder:
    with open(infile, "rU") as handle: 
        data = parse_fasta(handle)


    for header in data.keys():
        if header.split(".")[4] in toCheck:
            print(header)
        
