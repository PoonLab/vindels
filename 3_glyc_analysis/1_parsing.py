# script to parse the sequences retrieved for the glycosylation site analysis 
# collab with Adam

import sys 
import os 
import subprocess
from seqUtils import *

filename = "/home/jpalmer/PycharmProjects/glyc-analysis/gp120-glycs.fasta"


#if not os.path.isfile(filename):
    #print("this is not a proper file")
    #exit()
if not os.path.isdir("/home/jpalmer/PycharmProjects/glyc-analysis/1_filtered/"):
    os.mkdir("/home/jpalmer/PycharmProjects/glyc-analysis/1_filtered/")


infile = open(filename, "rU")

fasta = parse_fasta(infile)

subtypes = {}
filtered = {}

for header in fasta.keys():
    fields = header.split(".")

    status = fields[4]

    if status in ["acute_infection", "AIDS", "chronic"]:
        filtered[header] = fasta[header]
        if fields[0] != "-":
            if fields[0] in subtypes.keys():
                subtypes[fields[0]][header] = fasta[header]
            else:
                subtypes[fields[0]] = {header:fasta[header]}

'''
for sub in subtypes.keys():
    print(sub)
    print(len(subtypes[sub]))
'''

#write_fasta(filtered,"/home/jpalmer/PycharmProjects/glyc-analysis/1_filtered/gp120-glycs.fasta")

print(len(fasta))
print(len(filtered))
#print(filtered.keys())