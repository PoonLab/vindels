import sys 
import os 
import glob 
from seqUtils import *
import csv 

infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/ins-edit.csv", "rU")
icsv = csv.DictReader(infile,  delimiter="\t")

#print(data)

for line in icsv:

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

    print(anc + "\n" + tip)

    aa_anc = translate_nuc(anc, 0).replace("?", "-")
    aa_tip = translate_nuc(tip, 0).replace("?", "-")

    print(aa_anc)
    print(aa_tip)

    if aglycs:
        for n, x in enumerate(aglycs):
            aglycs[n] = int(x) - 1
    if tglycs:
        for n, x in enumerate(tglycs):
            tglycs[n] = int(x) - 1

    print(aglycs)
    print(tglycs)
    #for n, char in enumerate(anc):



    #List of tuples containing insertion sequences (0) and their end positions (1)
    insertions = []

    #Temporary strings used to store nucleotides from insertion sequences
    iTemp = ''

    ai = 0
    aindex = {}


    #FINDS ALL INDELS
    for n, achar in enumerate(aa_anc):
        tchar = aa_tip[n]

        #builds the alignment index
        if achar != '-':
            aindex.update({ai:n})
            ai += 1
            
        # builds the Insertion list
        #Case 1: no insertion, need to store and retrieve indel data
        if achar != "-":
            if iTemp:
                insertions.append(tuple((iTemp,n-1)))
                iTemp = ''

        #Case 2: insertion
        elif achar == "-":
            iTemp = iTemp + tchar

    #handles the END case
    if iTemp:
        insertions.append(tuple((iTemp, n - 1)))
    
    if len(aglycs) != len(tglycs):

        for seq, pos in insertions:
            print(seq)
            print(pos)

            iEnd = pos
            iStart = pos - (len(seq) - 1)

            #Determine which N-glycs are disrupted
            for end in nGlycs[a]:
                if end == '':
                    continue
                end = int(end)
                start = end - 3
                

                ngStart = rqindex[a][start]
                ngEnd = rqindex[a][end]

                # minus 1 inserted so the start refers to the first amino acid

                #CHECK FOR OVERLAP
                if not ((inStart > ngEnd and inEnd > ngEnd) or (inStart < ngStart and inEnd < ngStart)):
                    istr = str(inStart) + ":" + str(inEnd)
                    gstr = str(ngStart) + ":" + str(ngEnd)

                    if header in overlap.keys():
                        overlap[header][a].append(tuple((istr, gstr)))
                    elif a == 0:
                        overlap[header] = tuple(([tuple((istr, gstr))],[]))
                    else:
                        overlap[header] = tuple(([], [tuple((istr, gstr))]))'''
        

#CSNATLNCNNAINN------GSSSNNGSCSIDDGKIKEEMKN
#CSNATLNCNNAINNGSSSNNGSSSNNGSCSIDDGKIKEEMKN
#CSFNATAELKDKTQKVHSLFYRLDLVELNEDNNSNS--NTSMYRLINC
#CSFNATAELKDKTQKVHSLFYRLDLVELNEDNNSNSNTNTSMYRLINC