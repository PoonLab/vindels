import sys 
import os 
import glob 
from seqUtils2 import *
import csv 

infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/ins-edit.csv", "rU")
icsv = csv.DictReader(infile,  delimiter="\t")

#print(data)

ins_overlap = {}
aaTotal = {1:0, 2:0, 3:0, 4:0, 5:0}
ngTotal = {1:0, 2:0, 3:0, 4:0, 5:0}

def extractIndels(aa_anc, aa_tip, deletions=False):
    indels = []

    #Temporary strings used to store nucleotides from insertion sequences
    if deletions:
        swap = aa_anc
        aa_anc = aa_tip
        aa_tip = swap
    
    temp = ''

    #EXTRACT ALL INDELS 
    for n, achar in enumerate(aa_anc):
        tchar = aa_tip[n]
            
        # builds the Insertion list
        #Case 1: no insertion, need to store and retrieve indel data
        if achar != "-":
            if temp:
                indels.append(tuple((temp,n-1)))
                temp = ''

        #Case 2: insertion
        elif achar == "-":
            temp = temp + tchar

    #handles the END case
    if temp:
        indels.append(tuple((temp, n - 1)))
    return indels

count = 0
interfered = {}
for line in icsv:
    
    header = line["accno"]
    anc = line["ancestor"]
    tip = line["tipseq"]

    if line["anc.glycs"] != "NA":
        aglycs = line["anc.glycs"].split(",")
    else:
        aglycs = []
    if line["tip.glycs"] != "NA":
        tglycs = line["tip.glycs"].split(",")
    else:
        tglycs = []


    aa_anc = translate_nuc(anc, 0).replace("?", "-")
    aa_tip = translate_nuc(tip, 0).replace("?", "-")

    print(header)
    print(aa_anc)
    print(aa_tip)

    if "*" in aa_anc or "*" in aa_tip:
        count += 1
        print(aa_anc)
        continue

    vloop = int(line["vloop"])
    
    # for retrieving total aa count for each vloop
    temp = aa_anc.replace("-", "")
    aaCount = len(temp)

    # for retrieving the number of PNG amino acids for each vloop
    ngCount = 4*(len(re.findall("N[^P][ST][^P]", temp)))

    aaTotal[vloop] += aaCount
    ngTotal[vloop] += ngCount
    
    # used to change all nglyc positions to zero index
    if aglycs:
        for n, x in enumerate(aglycs):
            aglycs[n] = int(x) - 1
    if tglycs:
        for n, x in enumerate(tglycs):
            tglycs[n] = int(x) - 1

    #print(aglycs)
    #print(tglycs)
    #for n, char in enumerate(anc):

    #List of tuples containing insertion sequences (0) and their end positions (1)
    insertions = []

    #Temporary strings used to store nucleotides from insertion sequences
    iTemp = ''

    aindex = alignIndex(aa_anc)

    #EXTRACT ALL INDELS 
    insertions = extractIndels(aa_anc, aa_tip, False) 

    # CHECK FOR OVERLAP USING NGLYC POSITIONS 
    for seq, pos in insertions:
        inEnd = pos
        inStart = pos - (len(seq) - 1)
    
        #print(inStart)
        #print(inEnd)
        for ngStart in tglycs:
            if ngStart == '':
                continue
            ngStart = int(ngStart)
            ngEnd = ngStart + 3
                
            #CHECK FOR OVERLAP
            if not ((inStart > ngEnd and inEnd > ngEnd) or (inStart < ngStart and inEnd < ngStart)):
                istr = str(inStart) + ":" + str(inEnd)
                gstr = str(ngStart) + ":" + str(ngEnd)
                tpl = tuple((istr, gstr))

                # CHECK FOR NGLYC SITE CHANGE  
                inSeq = aa_anc[ngStart:ngEnd+1]
                beyond = aa_anc[ngEnd+1:].replace('-', '')

                dashes = inSeq.count('-')
                #print(beyond)
                #print(dashes)
                #print(inSeq)
                if inSeq[0] == "N" and '-' in inSeq and len(beyond) > dashes :
                    i = 0
                    inSeq = list(inSeq.replace("-", ""))
                    while i < dashes and inSeq:
                        inSeq.append(beyond[i])
                        i += 1
                    inSeq = ''.join(inSeq)
                #print(inSeq)
                test = re.search("N[^P][ST][^P]", inSeq)
                #print(test)
                if not test:
                    #print("FAIL")

                    #add it to the interfered list
                    duple = False
                    if header in interfered:
                        #print("HELLO")

                        # look for the insertion in the list ; if its there already, no need to add it again (prevents redundancy)
                        for x in interfered[header]:
                            if x[0] == tpl[0]:
                                duple = True
                        if not duple:
                            interfered[header].append(tpl)
                    else:
                        interfered[header] = [tpl]

infile.close()
print(interfered)
print(len(interfered))
# OUTPUT 
ioutput = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/interfered/insertions.csv", "w")
ioutput.write("header\tpos\n")
for header in interfered:
    outlist = []
    for tup in interfered[header]:
        ins, nglyc = tup
        pos = ins.split(":")[0]
        outlist.append(pos)
    ioutput.write(header + "\t" + ",".join(map(str, outlist)) + "\n")
ioutput.close()

print(aaTotal)
print(ngTotal)
outprops = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/interfered/ngprops.csv", "w")

props = [float(ngTotal[1]) / aaTotal[1],
        float(ngTotal[2]) / aaTotal[2],
        float(ngTotal[4]) / aaTotal[4],
        float(ngTotal[5]) / aaTotal[5]]
outprops.write("type,V1,V2,V3,V4,V5\n")
outprops.write("insertions,{},{},0.0,{},{}\n" .format(props[0],props[1],props[2],props[3]))


# DELETIONS 
# ------------------------

interfered = {}

infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/del-edit.csv", "rU")
dcsv = csv.DictReader(infile,  delimiter="\t")

interfered = {}
for line in dcsv:
    
    header = line["accno"]
    anc = line["ancestor"]
    tip = line["tipseq"]

    if line["anc.glycs"] != "NA":
        aglycs = line["anc.glycs"].split(",")
    else:
        aglycs = []
    if line["tip.glycs"] != "NA":
        tglycs = line["tip.glycs"].split(",")
    else:
        tglycs = []


    aa_anc = translate_nuc(anc, 0).replace("?", "-")
    aa_tip = translate_nuc(tip, 0).replace("?", "-")

    #print(header)
    #print(aa_anc)
    #print(aa_tip)

    if "*" in aa_anc or "*" in aa_tip:
        count += 1
        print(aa_anc)
        continue

    vloop = int(line["vloop"])
    
    # for retrieving total aa count for each vloop
    temp = aa_anc.replace("-", "")
    aaCount = len(temp)

    # for retrieving the number of PNG amino acids for each vloop
    ngCount = 4*(len(re.findall("N[^P][ST][^P]", temp)))

    aaTotal[vloop] += aaCount
    ngTotal[vloop] += ngCount
    
    # used to change all nglyc positions to zero index
    if aglycs:
        for n, x in enumerate(aglycs):
            aglycs[n] = int(x) - 1
    if tglycs:
        for n, x in enumerate(tglycs):
            tglycs[n] = int(x) - 1


    #List of tuples containing insertion sequences (0) and their end positions (1)
    deletions = []

    #Temporary strings used to store nucleotides from insertion sequences
    iTemp = ''

    aindex = alignIndex(aa_tip)

    #EXTRACT ALL INDELS 
    deletions = extractIndels(aa_anc, aa_tip, True) 
    print(header)
    print(aa_anc)
    print(aa_tip)
    print(deletions)
    # CHECK FOR OVERLAP USING NGLYC POSITIONS 
    for seq, pos in deletions:
        inEnd = pos
        inStart = pos - (len(seq) - 1)
    
        print(inStart)
        print(inEnd)
        for ngStart in aglycs:
            if ngStart == '':
                continue
            ngStart = int(ngStart)
            ngEnd = ngStart + 3
                
            #CHECK FOR OVERLAP
            if not ((inStart > ngEnd and inEnd > ngEnd) or (inStart < ngStart and inEnd < ngStart)):
                istr = str(inStart) + ":" + str(inEnd)
                gstr = str(ngStart) + ":" + str(ngEnd)
                tpl = tuple((istr, gstr))

                # CHECK FOR NGLYC SITE CHANGE  
                inSeq = aa_tip[ngStart:ngEnd+1]
                beyond = aa_tip[ngEnd+1:].replace('-', '')

                dashes = inSeq.count('-')
                print(beyond)
                print(dashes)
                print(inSeq)
                if inSeq[0] == "N" and '-' in inSeq and len(beyond) > dashes :
                    i = 0
                    inSeq = [inSeq.replace("-", "")]
                    while i < dashes and inSeq:
                        inSeq.append(beyond[i])
                        i += 1
                    inSeq = ''.join(inSeq)
                print(inSeq)
                test = re.search("N[^P][ST][^P]", inSeq)
                print(test)
                if not test:
                    print("FAIL")

                    #add it to the interfered list
                    duple = False
                    if header in interfered:
                        print("HELLO")

                        # look for the insertion in the list ; if its there already, no need to add it again (prevents redundancy)
                        for x in interfered[header]:
                            if x[0] == tpl[0]:
                                duple = True
                        if not duple:
                            interfered[header].append(tpl)
                    else:
                        interfered[header] = [tpl]
infile.close()

props = [float(ngTotal[1]) / aaTotal[1],
        float(ngTotal[2]) / aaTotal[2],
        float(ngTotal[4]) / aaTotal[4],
        float(ngTotal[5]) / aaTotal[5]]
outprops.write("deletions,{},{},0.0,{},{}\n".format(props[0],props[1],props[2],props[3]))
outprops.close()

print(interfered)
print(len(interfered))      

# OUTPUT 
doutput = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/interfered/deletions.csv", "w")
doutput.write("header\tpos\n")
for header in interfered:
    outlist = []
    for tup in interfered[header]:
        dl, nglyc = tup
        pos = dl.split(":")[0]
        outlist.append(pos)
    doutput.write(header + "\t" + ",".join(map(str, outlist)) + "\n")
doutput.close()

