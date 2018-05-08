import sys
from glob import glob
from seqUtils import *



folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/5MSAlignments/*.fasta")


for file in folder:

    name = file.split("/")[6].split('_MSA')[0]

    fasta = open(file, "r")
    
    #list format
    
    transposed = transpose_fasta(convert_fasta(fasta))
    
    pos = 0

    whitelist = []
    for x in transposed:
        pos += 1
    
        gaps = x.count("-")

        freq = float(gaps)/len(x)

        if freq < 0.95:
            whitelist.append(pos)

    fasta.close()


    fasta = open(file, "r")

    data = parse_fasta(fasta)
    
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/5_1_nogaps/"+ name +"_MSA_.fasta",'w')

    #dictionary format
    
    for i in data.keys():
        seq2 = ''
        for n, char in enumerate(data[i][:]):

            if n+1 in whitelist:
                seq2 += char

        output.write('>' + i + "\n")
        output.write(seq2 + "\n")


    fasta.close()
    output.close()
