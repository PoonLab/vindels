import re
from glob import glob
import os
from gotoh2 import *
from seqUtils import *

gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line

print(ntref)
aaref = translate_nuc(ntref, 0)

inpath = "/home/jpalmer/PycharmProjects/glyc-analysis/1_filtered/gp120-glycs.fasta"
outpath = "/home/jpalmer/PycharmProjects/glyc-analysis/2_pairwise/gp120.fasta"

if not os.path.isdir(os.path.dirname(outpath)):
    os.mkdir("/home/jpalmer/PycharmProjects/glyc-analysis/2_pairwise/")

infile = open(inpath,"rU")

pairwise = {}

filename = os.path.basename(inpath)

data = parse_fasta(infile)

unequal = []

output = open(outpath, "w")

for header in data:

    # PART 1 NUCLEOTIDE BASED ALIGNMENT TO REMOVE EXTRANEOUS SEQUENCE
    nt_pair = Aligner()
    nt_pair.set_model('HYPHY_NUC')
    nt_pair.is_global = False
    nt_pair.gap_open_penalty = 30
    nt_pair.gap_extend_penalty = 10

    result = nt_pair.align(ntref, data[header])

    left, right = get_boundaries(result[0])

    ntqry = result[1][left:right].replace("-","")


    # PART 2 CODON BASED NUCLEOTIDE ALIGNMENT 
    aaqry = translate_nuc(ntqry,0)

    aa_pair = Aligner()
    aa_pair.set_model('EmpHIV25')
    aa_pair.gap_extend_penalty = 10
    aa_pair.gap_open_penalty = 30
    aa_pair.is_global = True

    result2 = aa_pair.align(aaref,aaqry)

    newref = list(ntref)
    newqry = list(ntqry)

    # reads through the amino acid alignment and adds codon gaps to the proper locations
    for i in range(len(result2[0])):
        if result2[0][i] == '-':
            newref[i * 3:i * 3] = ['-', '-', '-']

        if result2[1][i] == '-':
            newqry[i * 3:i * 3] = ['-', '-', '-']

    finalRef = "".join(newref)
    finalQry = "".join(newqry)
    
    if len(newref) != len(newqry):
        unequal.append(header)
        print(header)
        print(ntqry)
        print("")
        print(ntref)
        print("")
        print(aaref)
        print("")
        print(aaqry)
        print("--------------------")
        print("")
        continue

    output.write(">" + header + '\n')
    output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')

