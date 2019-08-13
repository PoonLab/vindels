import sys 
import os 
import glob 
from seqUtils import *
import csv 

infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/ins-edit.csv", "rU")
icsv = csv.DictReader(infile,  delimiter="\t")

#print(data)

ins_overlap = {}
aaTotal = {1:0, 2:0, 3:0, 4:0, 5:0}
ngTotal = {1:0, 2:0, 3:0, 4:0, 5:0}

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

    #print(anc + "\n" + tip)

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
    
    # for retrieving the PNG site proportion of each vloop
    temp = aa_anc.replace("-", "")
    aaCount = len(temp)
    nglycs = re.findall("N[^P][ST][^P]", temp)
    #print(nglycs)
    ngCount = 4*(len(nglycs))

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
    
    for seq, pos in insertions:
        #print(seq)

        inEnd = pos
        inStart = pos - (len(seq) - 1)
    
        #print(inStart)
        #print(inEnd)
        for ngStart in tglycs:
            if ngStart == '':
                continue
            ngStart = int(ngStart)
            ngEnd = ngStart + 3
                
            #print(ngStart)
            #print(ngEnd)

            # minus 1 inserted so the start refers to the first amino acid

            #CHECK FOR OVERLAP
            if not ((inStart > ngEnd and inEnd > ngEnd) or (inStart < ngStart and inEnd < ngStart)):
                istr = str(inStart) + ":" + str(inEnd)
                gstr = str(ngStart) + ":" + str(ngEnd)
                #print(inStart)
                #print(inEnd)
                if header in ins_overlap.keys():
                    ins_overlap[header].append(tuple((istr, gstr)))
                else:
                    ins_overlap[header] = [tuple((istr, gstr))]
infile.close()

print(ins_overlap)
i_interfered = {}

infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/ins-edit.csv", "rU")
icsv = csv.DictReader(infile,  delimiter="\t")

for line in icsv:
    
    header = line["accno"]
    anc = line["ancestor"]
    tip = line["tipseq"]

    #print(anc + "\n" + tip)

    aa_anc = translate_nuc(anc, 0).replace("?", "-")
    aa_tip = translate_nuc(tip, 0).replace("?", "-")

    print(header)
    print(aa_anc)
    print(aa_tip)
    if header in ins_overlap: 
        if ins_overlap[header]:
            print(ins_overlap[header])
            for tpl in ins_overlap[header]:
                ngStart = int(tpl[1].split(":")[0])
                ngEnd = int(tpl[1].split(":")[1]) + 1
            
                print(ngStart)
                print(ngEnd)
                inSeq = aa_tip[ngStart:ngEnd]
                beyond = aa_tip[ngEnd:].replace('-', '')

                dashes = inSeq.count('-')

                if inSeq[0] == "N" and '-' in inSeq and len(beyond) >= dashes :
                    i = 0
                    inSeq = list(inSeq.replace("-", ""))
                    while i < dashes and inSeq:
                        inSeq.append(beyond[i])
                        i += 1
                    inSeq = ''.join(inSeq)

                test = re.search("N[^P][ST][^P]", inSeq)

                if not test:
                    print(ngStart)
                    print(ngEnd)
                    #add it to the i_interfered list
                    duple = False
                    if header in i_interfered:
                        for x in i_interfered[header]:
                            if x[0] == tpl[0]:
                                duple = True
                        if not duple:
                            i_interfered[header].append(tpl)
                        else:
                            i_interfered[header] = [tpl]


print(i_interfered)


'''
ioutput = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/interfered/insertions.csv", "w")
for header in ins_overlap:
    outlist = []
    for tup in ins_overlap[header]:
        pos = tup[0].split(":")[0]
        outlist.append(pos)
    ioutput.write(header + ":" + ",".join(map(str, outlist)) + "\n")

totals = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/interfered/glycs.csv", "w")
props = [float(ngTotal[1]) / aaTotal[1],
        float(ngTotal[2]) / aaTotal[2],
        #float(ngTotal[3]) / aaTotal[3],
        float(ngTotal[4]) / aaTotal[4],
        float(ngTotal[5]) / aaTotal[5]]
totals.write("insertions,{},{},{},{}\n".format(props[0],props[1],props[2],props[3]))


dlfile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/del-edit.csv", "rU")
dcsv = csv.DictReader(dlfile,  delimiter="\t")
del_overlap = {}
count =0
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

    #print(anc + "\n" + tip)

    aa_anc = translate_nuc(anc, 0).replace("?", "-")
    aa_tip = translate_nuc(tip, 0).replace("?", "-")
    
    if "*" in aa_anc or "*" in aa_tip:
        count += 1
        continue
    
    print(header)
    print(aa_anc)
    print(aa_tip)
    vloop = int(line["vloop"])
    
    # for retrieving the PNG site proportion of each vloop
    temp = aa_anc.replace("-", "")
    aaCount = len(temp)
    nglycs = re.findall("N[^P][ST][^P]", temp)

    ngCount = 4*(len(nglycs))

    aaTotal[vloop] += aaCount
    ngTotal[vloop] += ngCount

    #print(aa_anc)
    #print(aa_tip)
    
    # used to change all nglyc positions to zero index
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
    deletions = []

    #Temporary strings used to store nucleotides from insertion sequences
    dTemp = ''

    ti = 0
    tindex = {}

    #FINDS ALL INDELS
    for n, tchar in enumerate(aa_tip):
        achar = aa_anc[n]

        #builds the alignment index for the TIP SEQUENCE 
        if tchar != '-':
            tindex.update({ti:n})
            ti += 1
            
        # builds the Insertion list
        #Case 1: no insertion, need to store and retrieve indel data
        if tchar != "-":
            if dTemp:
                deletions.append(tuple((dTemp,n-1)))
                dTemp = ''

        #Case 2: insertion
        elif tchar == "-":
            dTemp = dTemp + achar

    #handles the END case
    if dTemp:
        deletions.append(tuple((dTemp, n - 1)))

    #print(deletions)
    
    for seq, pos in deletions:
        inEnd = pos
        inStart = pos - (len(seq) - 1)
    
        #print(inStart)
        #print(inEnd)

        for ngStart in tglycs:
            if ngStart == '':
                continue
            ngStart = int(ngStart)
            ngEnd = ngStart + 3
                
            #print(ngStart)
            #print(ngEnd)

                # minus 1 inserted so the start refers to the first amino acid

            #CHECK FOR OVERLAP
            if not ((inStart > ngEnd and inEnd > ngEnd) or (inStart < ngStart and inEnd < ngStart)):
                istr = str(inStart) + ":" + str(inEnd)
                gstr = str(ngStart) + ":" + str(ngEnd)
                if header in del_overlap.keys():
                    del_overlap[header].append(tuple((istr, gstr)))
                else:
                    del_overlap[header] = [tuple((istr, gstr))]
                    
print(count)
print(del_overlap)
print(len(del_overlap))


doutput = open("/home/jpalmer/PycharmProjects/hiv-withinhost/13_nglycs/interfered/deletions.csv", "w")
for header in del_overlap:
    outlist = []
    for tup in del_overlap[header]:
        pos = tup[0].split(":")[0]
        outlist.append(pos)
    doutput.write(header + ":" + ",".join(map(str, outlist)) + "\n")


# write on line 2 of the same glycs.csv file
props = [float(ngTotal[1]) / aaTotal[1],
        float(ngTotal[2]) / aaTotal[2],
        float(ngTotal[4]) / aaTotal[4],
        float(ngTotal[5]) / aaTotal[5]]
totals.write("deletions,{},{},{},{}\n".format(props[0],props[1],props[2],props[3]))'''

#CSNATLNCNNAINN------GSSSNNGSCSIDDGKIKEEMKN
#CSNATLNCNNAINNGSSSNNGSSSNNGSCSIDDGKIKEEMKN
#CSFNATAELKDKTQKVHSLFYRLDLVELNEDNNSNS--NTSMYRLINC
#CSFNATAELKDKTQKVHSLFYRLDLVELNEDNNSNSNTNTSMYRLINC
