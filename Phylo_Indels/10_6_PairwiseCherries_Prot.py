import subprocess
# TODO: use tempfile module
import os
from glob import glob
from seqUtils import *
from gotoh2 import *


files = glob("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_5_Translated_Cherries/*+.fasta")


for n in files:
    fasta = open(n, "r")

    data = parse_fasta2(fasta)

    subtype = n.split("/")[-1]
    output = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_6_Pairwise_Prot/" +subtype + "z", 'w')

    for x in data:

        if len(data[x]) < 2:
            print(x)

        seq1 = data[x][0]
        seq2 = data[x][1]


        pairwise = Aligner()
        pairwise.set_model('EmpHIV25')
        pairwise.gap_extend_penalty = 20
        pairwise.gap_open_penalty = 40
        pairwise.is_global=True


        nglyc1 = []
        nglyc2 = []



        for n in re.finditer("N[^P][ST][^P]", seq1):
            nglyc1.append(n.end()-1)
        for n in re.finditer("N[^P][ST][^P]", seq2):
            nglyc2.append(n.end()-1)

        string1 = ",".join(map(str,nglyc1))
        string2 = ",".join(map(str,nglyc2))

        result = pairwise.align(seq1, seq2)
        
        output.write(">" + x + "." + string1 + "." + string2 + '\n')
        output.write(">ref\n" + result[0] + "\n>query\n" + result[1] + '\n')


