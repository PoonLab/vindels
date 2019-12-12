import glob
import sys
import os 
from seqUtils import *

path = "/home/jpalmer/PycharmProjects/hiv-withinhost/"
file1 = open(path+"4MSA/hm-screen/111848/111848-1.fasta","rU")
file2 = open(path+"4MSA/hm-screen/111848/111848-2.fasta","rU")

incsv = open(path+"3RegionSequences/variable/111848.csv","rU")
csv = incsv.readlines()
data1 = parse_fasta(file1)
data2 = parse_fasta(file2)

out1 = open(path+"3RegionSequences/variable/111848-1.csv","w+")
out2 = open(path+"3RegionSequences/variable/111848-2.csv","w+")

for line in csv:
    header = line.split(",")[0]  
    if header in data1.keys():
        out1.write(line)
    if header in data2.keys():
        out2.write(line)


out1.close()
out2.close()
