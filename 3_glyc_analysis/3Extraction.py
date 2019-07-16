import re
from glob import glob
import os
from seqUtils import *
import numpy as np

gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line

c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]

fulldict = {}
cseqdict = {}
vseqdict = {}
vposdict = {}

infile = open("/home/jpalmer/PycharmProjects/glyc-analysis/2_pairwise/gp120.fasta","rU")
fasta = parse_fasta2(infile)

incorrect = []

for header, seq in fasta.items():  
    #extract the reference and query sequences 
    ref, query = seq

    left, right = get_boundaries(ref)
    query = query[left:right]

    #creates an alignment index to locate the variable regions 
    index = {}
    ri = 0
    for ai, x in enumerate(ref):
        if x != "-":
            index.update({ri:ai})
            ri += 1

    #------------
    #extracting the conserved region sequences for MSAs
    cseq = ""
    for c1, c2 in c_regions:
        cseq += query[index[c1]:index[c2]].replace("-","")
    cseqdict[header] = cseq
    print(ref)
    print("")
    print(query)
    print("")

    # extracting the whole gene sequence
    fulldict[header] = query.replace("-","")

    # extracting the variable loop sequences
    vseq=[]
    vpos = []
    for n1, n2 in v_regions:
        # cut out the variable region 
        vseq.append(query[index[n1]:index[n2]].replace("-",""))

        # locate the start and stop positions 
        start = str(len(query[:index[n1]].replace("-","")))
        end = str(len(query[:index[n2]].replace("-","")))
        
        vpos.append(tuple((start, end)))

    vseq = ",".join(vseq)
    #used for making the vseqdict which is used to create the library of variable region sequences 
    vseqdict[header] = vseq
    vposdict[header] = vpos




fullout = open("/home/jpalmer/PycharmProjects/glyc-analysis/3_sequences/full/gp120.csv","w")
vout = open("/home/jpalmer/PycharmProjects/glyc-analysis/3_sequences/variable/gp120.csv", "w")
cout = open("/home/jpalmer/PycharmProjects/glyc-analysis/3_sequences/conserved/gp120.fasta", "w")
fullout.write("header,seq,V1st,V1end,V2st,V2end,V3st,V3end,V4st,V4end,V5st,V5end\n")
vout.write("header,V1,V2,V3,V4,V5\n")
total = 0
patcount = 0
for header in fulldict.keys():
    #print(",".join(vposdict[header][1]))
    cout.write(">" + header + "\n" + cseqdict[header] + "\n")
    fullout.write(header + "," + fulldict[header])
    vout.write(header+ ",")
    
    for n in range(5):
        if n < 4:
            vout.write(vseqdict[header][n] + "," + ",".join(vposdict[header][n]) + ",")
            fullout.write("," + ",".join(vposdict[header][n]))
        else:
            vout.write(vseqdict[header][n] + "," + ",".join(vposdict[header][n]) + "\n")
            fullout.write("," + ",".join(vposdict[header][n]) +"\n")




