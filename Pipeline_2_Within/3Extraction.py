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

    print(file)

    unique = {}

    with open(file) as temp:
        fasta = parse_fasta2(temp)

    #os.mkdir("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + filename)
    #os.mkdir("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/conserved/" + filename)


    incorrect = []
    cpatdict = {}
    vpatdict = {}

    #print(fasta.items())
    for header, seq in fasta.items():
        #print(header)
        #print()
        accno = header.split(".")[4]
        patid = header.split(".")[3]
        #print(accno)

        ref, query = seq

        index = {}
        ri = 0
        for ai, x in enumerate(ref):
            if x != "-":
                index.update({ri:ai})
                ri += 1

        vseq=[]
        for n1, n2 in v_regions:
            vseq.append(query[index[n1]:index[n2]].replace("-",""))
        
        cseq = ""
        for c1, c2 in c_regions:
            cseq += query[index[c1]:index[c2]].replace("-","")


        if patid in unique.keys():
            #add sequence to the existing list
            if cseq in unique[patid].keys():
                unique[patid][cseq].append(header)

            #initialize new sequence as a list
            else:
                unique[patid][cseq] = [header]
        else:
            unique[patid] = {cseq:[header]}

            #if c1 != 1410:
            #    outputc.write(seq)
            #else:
            #    outputc.write(seq+"\n")
        '''if patid in cpatdict.keys():
            cpatdict[patid][header] = cseq
        else:
            cpatdict[patid] = {header: cseq}'''

        if patid in vpatdict.keys():
            vpatdict[patid][header] = ",".join(vseq)
        else:
            vpatdict[patid] = {header: ",".join(vseq)}
        

    for pat in unique.keys():
        outputc = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/conserved/" + filename + "/" +pat+".fasta", "w")
        outputv = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + filename + "/" + pat + ".csv", "w")
        outputv.write("header,V1,V2,V3,V4,V5\n")

        for seq in unique[pat].keys():

            header = unique[pat][seq][0]
            print(header)
            outputc.write(">"+header+"\n"+seq+"\n")

            outputv.write(header + "," + vpatdict[pat][header]+"\n")



