from glob import glob
import sys
from seqUtils import *




folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/0SequenceSets/*.fasta")


for file in folder:

    fasta = open(file, 'r')

    filename = file.split("/")[-1]
    dict1 = {}
    unique = {}

    output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/1FilteredSeqs/" + filename, "w")
    data = parse_fasta(fasta)

    # screen for sequences missing subtype, collection year, and those without sampling collection info
    for x in data:
        header = x.split(".")
        print(x)
        if header[0] != "-" and header[2] != "-" and len(data[x]) > 1400 and (header[6] != "-" or header[7] != "-" or header[8] != "-" or header[9] != "-"):
            dict1[x] = data[x]


    # screen for identical sequences
    for y in dict1:
        unique[dict1[y]] = y

    print(len(data))
    print(len(dict1))
    print(len(unique))

    for seq in unique:
        output.write(">" + unique[seq] + "\n" + seq + "\n")






