from glob import glob
import sys
from seqUtils import *




folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/1SequenceSets/*.fasta")


for file in folder:

    fasta = open(file, 'r')

    filename = file.split("/")[-1]

    output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/2FilteredSeqs/" + filename, "w")
    data = parse_fasta(fasta)

    for x in data:
        header = x.split(".")
        if header[6] != "-" or header[7] != "-" or header[8] != "-" or header[9] != "-":
            output.write(">"+x+"\n"+data[x]+"\n")




