import re
import os
import re
from gotoh2 import *
from seqUtils import *

gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line

def isCysteine(codon):
    if codon.upper() in ["TGT", "TGC"]:
        return True
    else:
        return False
aaref = translate_nuc(ntref, 0)

inpath = "/home/jpalmer/PycharmProjects/glyc-analysis/2_pairwise/gp120.fasta"
outpath = "/home/jpalmer/PycharmProjects/glyc-analysis/3_glycs/glycs.csv"

if not os.path.isdir(os.path.dirname(outpath)):
    os.mkdir("/home/jpalmer/PycharmProjects/glyc-analysis/3_glycs/")

infile = open(inpath,"rU")

filename = os.path.basename(inpath)
vseqdict = {}
data = parse_fasta2(infile)

output = open(outpath, "w")
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]



for header, seq  in data.items():
    ref, query = seq

    #creates an alignment index to locate the variable regions 
    index = {}
    ri = 0
    for ai, x in enumerate(ref):
        if x != "-":
            index.update({ri:ai})
            ri += 1

    #extracts the variable region sequences and loads them into VSEQ
    vseq=[]
    for n1, n2 in v_regions:
        vseq.append(query[index[n1]:index[n2]].replace("-",""))

    xCount = 0
    for idx in range(4):
        if not isCysteine(vseq[idx][0:3]):
            xCount += 1
    
    if xCount > 1:
        print(vseq)
        continue
    
    #vseq = ",".join(vseq)

    vseqdict[header] = vseq


for header in vseqdict.keys():

    for seq in vseqdict[header]:
        print(seq)
        aaseq = translate_nuc(seq,0)
        print(aaseq)
        nglycs = re.findall("N[^P][ST][^P]", aaseq)
        print(nglycs)


    #output.write(header+","+vseqdict[header]+"\n")
        






