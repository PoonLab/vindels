import sys
from glob import glob
from seqUtils import *

import os

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/*.fasta")
count = 0

for infile in folder:

    name = infile.split("/")[6].split('.fasta')[0]

    fasta = open(infile, "r")
    print(os.path.basename(infile))
    #list format
    
    transposed = transpose_fasta(convert_fasta(fasta))

    whitelist = []
    for pos, x in enumerate(transposed):

        gaps = x.count("-")

        freq = float(gaps)/len(x)

        if freq < 0.95:
            whitelist.append(pos)

    fasta.close()


    fasta = open(file, "r")

    data = parse_fasta(fasta)
    count += len(data)
    output = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/5_1_final/"+ name +"_.fasta",'w')

    #dictionary format
    
    for i in data.keys():
        accno = i[-8:]

        if "." in accno:
            accno = accno.split(".")[1]

        seq2 = ''
        for n, char in enumerate(data[i][:]):

            if n in whitelist:
                seq2 += char

        output.write('>' + accno + "\n")
        output.write(seq2 + "\n")


    fasta.close()
    output.close()
print(count)