from glob import glob
from seqUtils import *
import re
files = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_Cherries/*.fasta")



for n in files:

    fasta = open(n, "r")

    data = parse_fasta2(fasta)


    subtype = n.split("/")[-1].split(".")[0]
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_5_Translated_Cherries/" +subtype+ ".fasta", 'w')

    for x in data:

        if(len(data[x]) < 2):
            continue

        aaSeq1 = translate_nuc(data[x][0],0)
        aaSeq2 = translate_nuc(data[x][1],0)

        output.write(">"+ x + "\n>ref\n" + aaSeq1 + "\n>query\n" + aaSeq2 + "\n")