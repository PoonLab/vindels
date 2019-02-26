import sys
from seqUtils import *
from glob import glob
#processes the two separate sequence files retrieved from Genbank
#reads line by line, finds and prints accession number, finds and prints collection date if it exists


# these files contain the raw data from GenBank

folder = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/3_RegionSequences/full_length/*.fasta')


#prints the date to this
output =  open("/home/jpalmer/PycharmProjects/hiv-withinhost/accno_all.txt", 'w')
#output1 = open("accno1.txt", 'w')
#output2 = open("accno2.txt", 'w')
seq = []
for file in folder:
    #print(file)

    fasta = open(file, 'r')

    data = parse_fasta(fasta)
    for n in data:
        seq.append(n)
        output.write(n.split(".")[4] + "\n")

print(len(seq))