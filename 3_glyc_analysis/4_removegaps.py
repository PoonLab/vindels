import sys
import os
from glob import glob
from seqUtils import *



infile = open("/home/jpalmer/PycharmProjects/glyc-analysis/4_msa/gaps-in/gp120.msa","rU")
output = open("/home/jpalmer/PycharmProjects/glyc-analysis/4_msa/gp120.msa",'w')
count = 0

transposed = transpose_fasta(convert_fasta(infile))
    
whitelist = []
for pos, x in enumerate(transposed):
    gaps = x.count("-")
    freq = float(gaps)/len(x)
    if freq < 0.95:
        whitelist.append(pos)
infile.close()



with open("/home/jpalmer/PycharmProjects/glyc-analysis/4_msa/gaps-in/gp120.msa","rU") as fasta:
    data = parse_fasta(fasta)

for header in data.keys():
    seq = ''
    for n, char in enumerate(data[header][:]):
        if n in whitelist:
            seq += char

    output.write('>' + header + "\n")
    output.write(seq + "\n")
