import subprocess
# TODO: use tempfile module
import os
from glob import glob
from seqUtils import *
from gotoh2 import *



files = glob("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_Cherries/*.fasta")

for n in files:
    fasta = open(n, "r")
    data = parse_fasta2(fasta)
    subtype = n.split("/")[-1]
    output = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_1_Pairwise_Nt/" + subtype + "z", 'w')

    for x in data:

        if len(data[x]) < 2:
            print(x)

        seq1 = data[x][0]
        seq2 = data[x][1]

        pairwise = Aligner()
        pairwise.set_model('HYPHY_NUC')
        pairwise.gap_extend_penalty = 5
        pairwise.gap_open_penalty = 30
        pairwise.is_global = True
        result = pairwise.align(seq1, seq2)
        output.write(">" + x + "." + '\n')
        output.write(">ref\n" + result[0] + "\n>query\n" + result[1] + '\n')