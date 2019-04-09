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

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/2PairwiseAA/*.fasta")

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
            vseq.append(str(len(query[:index[n1]].replace("-",""))))
            vseq.append(str(len(query[:index[n2]].replace("-",""))))
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
    dates = [[], [], [], [],[]]
    for header in full[pat].keys():
        fields = header.split(".")

        #these are the four date fields to check 
        dates[0].append(fields[5])
        dates[1].append(fields[6])
        dates[2].append(fields[7])
        dates[3].append(fields[8])
        #this one is just a count of how many timepoints
        dates[4].append(fields[9])
    
    #ensure that the patient has 5 time points or more 
    if int(dates[4][0]) < 5:
        print(pat + " has too few timepoints")
        continue
    
    # unique = list of 4 elements, the number of unique dates in each field    
    unique = []
    for x in range(4):
        dateset = set(dates[x])
        unique.append(len(dateset))

    #skip the patient if they contain no unique timepoints 
    if not any(i > 1 for i in unique):
        print(pat + " has no unique timepoints in ALL date fields")
        continue
    
    #find the field with the most unique timepoints 
    #completes the bestIdx variable, with the best timepoint 
    bestIdx = -1
    for n, x in enumerate(unique):
        if x > 1 and x > bestIdx:
            bestIdx = n

    #used to detect and fix negative date values in the chosen date list (dates[bestIdx])
    hasNeg = False
    lowest = 0
    for n, j in enumerate(dates[bestIdx]):
        
        #simple filter to first get rid of any 2222 values 
        if j == "2222":
            dates[bestIdx][n] = None
            j = "-"

        #this will check every non '-' value to see if its negative
        #if its negative, record it in hasNeg and find the lowest negative value 
        if j != "-":
            dates[bestIdx][n] = int(j)
            negative = re.search('-\d*',j)
            if negative is not None:
                hasNeg = True
                if int(j) < lowest:
                    lowest = int(j)

    # for the date lists containing negative values, recenter all valid dates to 0  
    if hasNeg: 
        print("HELLO THERE")
        print(pat)
        print(len(dates[bestIdx]))
        for n, k in enumerate(dates[bestIdx]):
            try:
                k = int(k)
            except:
                continue
            
            dates[bestIdx][n] = k - lowest
         
    
    #skips the entire patient if no valuable dates are found
    if all(i in ("-", "2222", None) for i in dates[bestIdx]):
        print(pat + " had no valuable date information")
        continue

    outputfull = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/" + pat + ".fasta","w")
    outputv = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + pat + ".csv", "w")
    outputv.write("header,V1,start,stop,V2,start,stop,V3,start,stop,V4,start,stop,V5,start,stop\n")

    patcount += 1
    for n, header in enumerate(full[pat].keys()):
        fields = header.split(".")
        
        #skips the sequences that do not contain a proper date 
        if dates[bestIdx][n] in ("-", "2222", None):
            continue

        # fields[0:5] . number of time points . selected time scale . which time scale was chosen
        hedit = ".".join(fields[0:5]) + "." + fields[9] + "." + str(bestIdx) + "_" + str(dates[bestIdx][n])
        outputfull.write(">" + fields[4] +"_" + str(dates[bestIdx][n]) + "\n" + full[pat][header] + "\n")
        outputv.write(fields[4] + "_" + str(dates[bestIdx][n]) + "," + vpatdict[pat][header] + "\n")
print(patcount)





