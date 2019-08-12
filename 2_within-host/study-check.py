import os 
import sys 
from seqUtils import * 
from glob import glob 
#print(glob.__file__)

infolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta")

toCheck = ["FJ919967", "FJ496214"]

testRange = [ "KC149"+str(x) for x in range(260,300) ]

for infile in infolder:
    with open(infile, "rU") as handle: 
        data = parse_fasta(handle)

    filename = os.path.basename(infile)
    '''if filename == "111848.fasta":
        unique = set()
        for header in data:
            unique.add(header.split(".")[4])
        
        print(unique)'''
    print(filename)
    for header in data.keys():
        print(header.split("_")[1])
        
        #if header.split(".")[4] in toCheck:
        #    print(header)
        
