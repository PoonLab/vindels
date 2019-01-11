import re
from glob import glob
import os
from seqUtils import *
import numpy as np


def samecheck(list):
    return all(x == list[0] for x in list)

gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line

c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/2_2_pairwiseAA/*.fasta")

overall = {}
vpatdict = {}
studyno = 0
full = {}
for file in folder:
    study = {}
    filename = file.split("/")[-1].split(".")[0]
    studyno += 1
    #print(file)

    #will contain all patients for the given study
    #unique = {}

    with open(file) as temp:
        fasta = parse_fasta2(temp)

    incorrect = []

    for header, seq in fasta.items():
        accno = header.split(".")[4]
        patid = header.split(".")[3]

        ref, query = seq

        index = {}
        ri = 0
        for ai, x in enumerate(ref):
            if x != "-":
                index.update({ri:ai})
                ri += 1

        #loads the variable regions into VSEQ
        vseq=[]
        for n1, n2 in v_regions:
            vseq.append(query[index[n1]:index[n2]].replace("-",""))
        vseq = ",".join(vseq)

        '''cseq = ""
        for c1, c2 in c_regions:
            cseq += query[index[c1]:index[c2]].replace("-","")'''

        #load the study-wise dictionary
        if patid in study.keys():
            study[patid][header] = query.replace("-","")
        else:
            study[patid] = {header:query.replace("-","")}

        #used for generated concatenated conserved sequences for MSAs
        '''if patid in unique.keys():
            #add sequence to the existing list
            if cseq in unique[patid].keys():
                unique[patid][cseq].append(header)

            #initialize new sequence as a list
            else:
                unique[patid][cseq] = [header]
        else:
            unique[patid] = {cseq:[header]}'''


        if patid in vpatdict.keys():
            vpatdict[patid][header] = vseq

        else:
            vpatdict[patid] = {header: vseq}

    #dump the entire study-wise dictionary into the full dictionary
    #overwrite a patients entry with the new one if the new study has a larger entry
    for patid in study:
        if patid in full.keys():
            #check whether the new one is bigger
            if len(study[patid]) > len(full[patid]):
                full[patid] = study[patid]
        else:
            full[patid] = study[patid]


patcount = 0
for pat in full.keys():
    info = [[], [], [], [],[]]
    for header in full[pat].keys():
        fields = header.split(".")
        info[0].append(fields[5])
        info[1].append(fields[6])
        info[2].append(fields[7])
        info[3].append(fields[8])
        info[4].append(fields[9])
    
    #ensure that the patient has 5 time points or more 
    if int(info[4][0]) < 5:
        print(pat + " has too few timepoints")
        continue


    bool = []
    for x in range(4):
        bool.append(len(set(info[x])) > 1)
    if not any(bool):
        print(pat + " has no unique timepoints")
        continue
    
    #find the first location where bool = TRUE (idx), AKA first field containing unique timepoint values
    idx = -1
    for n, x in enumerate(bool):
        if x:
            idx = n
            break

    # arbitrarily choose the first populated time field

    outputfull = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/" + pat + ".fasta","w")
    outputv = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + pat + ".csv", "w")
    outputv.write("header,V1,V2,V3,V4,V5\n")

    patcount += 1
    for header in full[pat].keys():
        fields = header.split(".")
        #print(fields)
        if "-" in fields[idx+5]:
            continue

        # fields[0:5] . selected time scale . number of time points . which time scale was chosen
        hedit = ".".join(fields[0:5]) + "." + fields[9] + "." + str(idx) + "_" + fields[idx+5]
        #print(fields[9])
        outputfull.write(">" + hedit + "\n" + full[pat][header] + "\n")
        outputv.write(hedit + "," + vpatdict[pat][header] + "\n")
print(patcount)

'''
for pat in full.keys():
        outputfull = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/" + pat + ".fasta", "w")
        outputv = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + pat + ".csv", "w")
        outputv.write("header,V1,V2,V3,V4,V5\n")
        for header in full[pat].keys():
            outputfull.write(">" + header + "\n" + full[pat][header] + "\n")
            outputv.write(header + "," + vpatdict[pat][header]+"\n")'''










