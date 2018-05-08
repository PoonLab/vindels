import sys
from glob import glob
from seqUtils import *


folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/5_1_nogaps/*.fasta")

for file in folder:
    input = open(file,'r')
    name = file.split("/")[6].split("_MSA_.fasta")[0]

    print(name)
    #list format
    transposed = transpose_fasta(convert_fasta(input))

    pos = 0

    for x in transposed:
        pos += 1
        gaps = x.count("-")

        #calculates the proportion of gaps in the given position
        freq = float(gaps)/len(x)

        print(pos, freq)

        if freq < 0.2:
            stop = pos
            #print(file + "\n" + "Cut Point:" + str(stop))
            break





    #used to print the lists of gap frequencies of each MSA

    print(stop)
    input.close()


    input = open(file, 'r')
    data = parse_fasta(input)

    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/5_2_Shortened/"+ name +"_MSA-2.fasta",'w')


    #dictionary format

    for i in data.keys():
        data[i] = data[i][stop-1:]

        output.write(">" + i + "\n")
        output.write(data[i] + "\n")

    input.close()
