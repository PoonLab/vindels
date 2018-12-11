import re
from glob import glob
import os
from seqUtils import *


gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line

c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/2_2_pairwiseAA/*.fasta")

for file in folder:
    filename = file.split("/")[-1].split(".")[0]

    with open(file) as temp:
        fasta = parse_fasta2(temp)

    #os.mkdir("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/conserved/" + filename)

    outputv = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + filename + ".csv", "w")
    outputv.write("V1,V2,V3,V4,V5\n")

    incorrect = []
    patdict = {}

    #print(fasta.items())
    for header, seq in fasta.items():
        print(header)
        accno = header.split(".")[4]
        patid = header.split(".")[3]
        print(accno)

        ref, query = seq

        index = {}
        ri = 0
        for ai, x in enumerate(ref):
            if x != "-":
                index.update({ri:ai})
                ri += 1


        for n1, n2 in v_regions:
            seq = query[index[n1]:index[n2]]

            if n1 != 1377:
                outputv.write(seq.replace("-", "") + ",")
            else:
                outputv.write(seq.replace("-","") + "\n")



        seq = ""
        for c1, c2 in c_regions:
            seq += query[index[c1]:index[c2]]

            #if c1 != 1410:
            #    outputc.write(seq)
            #else:
            #    outputc.write(seq+"\n")
        if patid in patdict.keys():
            patdict[patid][header] = seq
        else:
            patdict[patid] = {header: seq}


    for pat in patdict:
        outputc = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/conserved/" + filename + "/" +pat+".fasta", "w")

        for header in patdict[pat]:
            outputc.write(">"+header+"\n"+patdict[pat][header]+"\n")