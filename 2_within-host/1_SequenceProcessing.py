from glob import glob
import sys
from seqUtils import *
import os



folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/0SequenceSets/*.fasta")

subtypes = {}
location = {}
Cstudies = []
minlen = 1000000
total = 0 
for infile in folder:

    fasta = open(infile, 'r')

    filename = os.path.basename(infile)

    filter = {}
    byStudy = set()
    #output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/1FilteredSeqs/" + filename, "w")
    data = parse_fasta(fasta)
    
    # filter : screen for sequences missing subtype, collection year, and those without sampling collection info
    for x in data:
        print(len(data[x]))
        continue

        header = x.split(".")
        total += 1 
        if len(data[x]) < minlen:
            minlen = len(data[x])
        subtypes[header[0]] = subtypes.get(header[0],0) + 1
        location[header[1]] = location.get(header[1],0) + 1
        byStudy.add(header[0])
        if header[0] != "-" and header[2] != "-" and len(data[x]) > 1400 and (header[6] != "-" or header[7] != "-" or header[8] != "-" or header[9] != "-"):
            filter[x] = data[x]

    # screen for identical sequences
    #for header in filter:
        #output.write(">" + header + "\n" + filter[header] + "\n")

    for x in byStudy:
        if x == "C":
            Cstudies.append(filename)
    #print(filter.keys())
print(minlen)
print(total)
print(subtypes)
print(location)
print(sorted(Cstudies))





