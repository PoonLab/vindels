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

#print(ntref)
aaref = translate_nuc(ntref, 0)

inpath = "/home/jpalmer/PycharmProjects/glyc-analysis/7_prnseq/gp120.fasta"
outpath = "/home/jpalmer/PycharmProjects/glyc-analysis/8_aapairwise/gp120.fasta"

if not os.path.isdir(os.path.dirname(outpath)):
    os.mkdir("/home/jpalmer/PycharmProjects/glyc-analysis/8_aapairwise/")

infile = open(inpath,"rU")

pairwise = {}

filename = os.path.basename(inpath)

data = parse_fasta(infile)

unequal = []

output = open(outpath, "w")

for header in data:
    

    # PART 2 CODON BASED NUCLEOTIDE ALIGNMENT 

    aaqry = translate_nuc(data[header],0)
    aa_pair = Aligner()
    aa_pair.set_model('EmpHIV25')
    aa_pair.gap_extend_penalty = 10
    aa_pair.gap_open_penalty = 30
    aa_pair.is_global = True

    result = aa_pair.align(aaref,aaqry)

    if ("*" in result[0]) or ("*" in result[1]):
        continue


        # block below used to check the best reading frames 
        # found that they still suck overall
        # therefore, chose to ignore these sequences 
        '''
        lowest = 10^8
        for i in range(3):
            aaqry = translate_nuc(data[header],i)
            aa_pair = Aligner()
            aa_pair.set_model('EmpHIV25')
            aa_pair.gap_extend_penalty = 10
            aa_pair.gap_open_penalty = 30
            aa_pair.is_global = True

            result = aa_pair.align(aaref,aaqry)
            if result[2] < lowest:
                choice = i
                lowest = result[2]
        
        aaqry = translate_nuc(data[header],choice)
        aa_pair = Aligner()
        aa_pair.set_model('EmpHIV25')
        aa_pair.gap_extend_penalty = 10
        aa_pair.gap_open_penalty = 30
        aa_pair.is_global = True
        result = aa_pair.align(aaref,aaqry)

        print(result[0][:100])
        print(result[1][:100])
        print("")
        '''
    print(result[0])
    print(result[1])
    print("")
    
    output.write(">" + header + '\n')
    output.write(">ref\n" + result[0] + "\n>query\n" + result[1] + '\n')

output.close()
    
