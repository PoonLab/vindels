import sys
from seqUtils import *
from glob import glob
#processes the two separate sequence files retrieved from Genbank
#reads line by line, finds and prints accession number, finds and prints collection date if it exists


# these files contain the raw data from GenBank

folder = glob('/home/jpalme56/PycharmProjects/hiv-evolution-master/5_1_final/*.fasta')


#prints the date to this
output =  open("accnototal.txt", 'w')
output1 = open("accno1.txt", 'w')
output2 = open("accno2.txt", 'w')
seq = []
for file in folder:
    #print(file)

    fasta = open(file, 'r')

    data = parse_fasta(fasta)
    for n in data:
        seq.append(n)
print(len(seq))
first = len(seq)/2
second = len(seq) - first
print(first)
print(second)

for i in range(first):
    output1.write(seq[i] + "\n")

seq2 = seq[first:]
for j in range(second):
    output2.write(seq2[j] + "\n")

for z in seq:
    output.write(z + "\n")