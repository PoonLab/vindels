from glob import glob
import sys
from seqUtils import *




folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/0SequenceSets/*.fasta")


for file in folder:

    fasta = open(file, 'r')

    filename = file.split("/")[-1]

    filter = {}

    output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/1FilteredSeqs/" + filename, "w")
    data = parse_fasta(fasta)

    # filter : screen for sequences missing subtype, collection year, and those without sampling collection info
    for x in data:
        header = x.split(".")

        if header[0] != "-" and header[2] != "-" and len(data[x]) > 1400 and (header[6] != "-" or header[7] != "-" or header[8] != "-" or header[9] != "-"):
            filter[x] = data[x]

    #print(filter)
    # screen for identical sequences
    for header in filter:
        output.write(">" + header + "\n" + filter[header] + "\n")







