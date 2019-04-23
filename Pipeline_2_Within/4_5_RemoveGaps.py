import sys
import os
from glob import glob
from seqUtils import *



folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta")
count = 0

for infile in folder:
    filename = os.path.basename(infile)
    with open(infile, "r") as fasta:
        transposed = transpose_fasta(convert_fasta(fasta))
    
    whitelist = []
    for pos, x in enumerate(transposed):
        gaps = x.count("-")
        freq = float(gaps)/len(x)
        if freq < 0.95:
            whitelist.append(pos)

    fasta = open(infile, "r")

    with open(infile, "r") as fasta:
        data = parse_fasta(fasta)
    
    count += len(data)
    output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/4_1MSAedited/"+ filename,'w')

    #dictionary format
    for header in data.keys():
        seq = ''
        for n, char in enumerate(data[header][:]):
            if n in whitelist:
                seq += char

        output.write('>' + header + "\n")
        output.write(seq + "\n")

    output.close()
print(count)