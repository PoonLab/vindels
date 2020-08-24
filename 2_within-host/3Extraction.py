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

nov_regions = [(0,78),(78,198),(495,603),(762,864),(987,1020)]

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/2PairwiseAA/*.fasta")
print([os.path.basename(x) for x in folder])
print(len(folder))
vseqdict = {}
studyno = 0
full = {}

for infile in folder:
    
    v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]
    patdict = {}
    filename = os.path.basename(infile)
    
    # used for handling the VN data set 
    if filename == "novitsky.fasta":
        handleNov = True
    else:
        handleNov = False

    #overwrites the v_regions indexing list with the novitsky one, if the file is Vlad's
    if handleNov:
        v_regions = nov_regions

    studyno += 1
    #print(file)

    #will contain all patients for the given study
    #unique = {}

    with open(infile) as temp:
        fasta = parse_fasta2(temp)

    incorrect = []
    
    for header, seq in fasta.items():
        if not handleNov:
            accno = header.split(".")[4]
            patid = header.split(".")[3]
        else:
            header = header.strip("\r._")
            fields = header.split("_")
            patid = fields[0]
            if patid == "0Ref":
                continue
            # NEW HEADER FORMAT 
            # Subtype . Identifier . Identifier #2 . Date  
            header = "C.-.-."+fields[0]+"."+fields[2]+".-.-_"+fields[1]
            

        #extract the reference and query sequences 
        ref, query = seq

        query = query.strip("\r")


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
            vseq.append(str(len(query[:index[n1]].replace("-",""))))
            vseq.append(str(len(query[:index[n2]].replace("-",""))))
        vseq = ",".join(vseq)


        #load the patient-wise dictionary
        if patid in patdict.keys():
            patdict[patid][header] = query.replace("-","")
        else:
            patdict[patid] = {header:query.replace("-","")}
        

        #used for making the vseqdict which is used to create the library of variable region sequences 
        if patid in vseqdict.keys():
            vseqdict[patid][header] = vseq
        else:
            vseqdict[patid] = {header: vseq}

    #dump the entire study-wise dictionary into the full dictionary
    #overwrite a patients entry with the new one if the new study has a larger entry
    for patid in patdict:
        if patid in full.keys():
            #only overwrite and replace the existing entry in FULL if the new entry is bigger
            if len(patdict[patid]) > len(full[patid]):
                #print(patid)
                #print("old: "+str(len(full[patid])))
                #print("new: "+str(len(patdict[patid])))
                full[patid] = patdict[patid]
                #print("OLD: "+str(len(full[patid])))
                #print("NEW: "+str(len(patdict[patid])))
        #if patient hasnt been loaded yet, make a new entry
        else:
            full[patid] = patdict[patid]


total = 0
patcount = {}
seqcount = 0
subcount = {}
tp = []
for pat in full.keys():
    print(pat)

    delcount = 0
    # FILTER OUT PATIENT DATA SETS 
    # =========================================
    if not pat.isalpha():

        # dates = [tp1, tp2, tp3, tp4]
        dates = [[], [], [], []]
        
        timepoints = None
        for n, header in enumerate(full[pat].keys()):
            #print(header)
            fields = header.split(".")

            #these are the four date fields to check 
            for i in range(4):
                dates[i].append(fields[i+5])

            #this one is just a count of how many timepoints
            if n == 0:
                timepoints = int(fields[9])
                print(timepoints)
        
        # unique = list of 4 elements
        # counts the number of unique dates in each date field    
        unique = []
        for x in range(4):
            dateset = set(dates[x])
            unique.append(len(dateset))
        print(unique)
        # NO UNIQUE TIMEPOINTS 
        # skip the patient if they contain no unique timepoints 
        if not any(i > 1 for i in unique):
            print(pat + " has no unique timepoints in ALL date fields")
            continue
        
        #find the field with the most unique timepoints 
        #completes the bestIdx variable with the date field containing the most unique timepoints 
        bestIdx = -1
        for n, x in enumerate(unique):
            if x > 1:
                print(n)
                if x > bestIdx:
                    bestIdx = n
        tp.append(unique[bestIdx])
        bestIdx = bestIdx + 5

        #ensure that the patient has 5 time points or more 
        if timepoints < 5:
            print(pat + " has too few timepoints")
            continue


        
        # NEGATIVES --------------------------
        #used to detect and fix negative date values in the chosen date list (dates[bestIdx])
        hasNeg = False
        lowest = 0
        chosenDates = []
        #print(type(full[pat].keys()))
        for n, header in enumerate(list(full[pat].keys())[:]):
            fields = header.split(".")
            date = fields[bestIdx]
            
            #simple filter to get rid of any 2222 values 
            if date == "2222" or date == "8888" or date == "-":
                del full[pat][header]
                del vseqdict[pat][header]
                delcount += 1


            #this will check every non '-' value to see if its negative
            #if its negative, turn on hasNeg switch and find the lowest negative value 
            else:
                chosenDates.append(date)
                negative = re.search('-\d*',date)
                date = int(date)
                if negative is not None:
                    hasNeg = True
                    if date < lowest:
                        lowest = date

        # for the date lists containing negative values, recenter all valid dates to 0  
        if hasNeg: 
            for n, header in enumerate(full[pat].keys()):
                fields = header.split(".")
                date = fields[bestIdx]
                try:
                    date = int(date)
                except:
                    print("PROBLEM WITH DATE: "+date)
                    continue

                #modify the date field and reappend each element to the dictionary
                fields[bestIdx] = str(date - lowest)
                newheader = ".".join(fields)
                full[pat][newheader] = full[pat].pop(header)
                vseqdict[pat][newheader] = vseqdict[pat].pop(header)
        
        #skips the entire patient if no valuable dates are found
        if all(i in ("-", "2222","8888", None) for i in chosenDates):
            print(pat + " had no valuable date information")
            #print(chosenDates)
            continue

        # skips the entire patient if the date range is less than 200 days 
        chosenDates = [int(x) for x in chosenDates]
        rng = max(chosenDates) - min(chosenDates)
        if rng < 200:
            print(pat +  " has too small of a date range: "+str(rng))
            continue
    else:
        unique = set()
        alldates = []
        for header in full[pat]:
            date = header.split("_")[1]
            alldates.append(date)
            unique.add(date)

        if len(unique) < 5: 
            print(pat + " has too few timepoints")    
            continue

        # skips the entire patient if the date range is less than 200 days 
        alldates = [int(x) for x in alldates]
        rng = max(alldates) - min(alldates)
        if rng < 200:
            print(pat +  " has too small of a date range: "+str(rng))
            continue

    if len(vseqdict[pat]) == 0 or len(full[pat]) == 0:
        continue    


    # SUCCESS IF REACHED THIS POINT 
    # ======================================
    
    print(pat)
    
    # outputfull = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/" + pat + ".fasta","w")
    # outputv = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/" + pat + ".csv", "w")
    # outputv.write("header,V1,start.1,stop.1,V2,start.2,stop.2,V3,start.3,stop.3,V4,start.4,stop.4,V5,start.5,stop.5\n")
    
    for n, header in enumerate(full[pat].keys()):

        # LANL sequence data sets 
        if not pat.isalpha():
            fields = header.split(".")

            date = fields[bestIdx]
            #skips the sequences that do not contain a proper date 
            if date in ("-", "2222", '8888',None):
                print("SOMETHINGS WRONG")
                print(pat)
                print(date)
            
            # FULL HEADER
            # fields[0:5] . number of time points . which time scale was chosen .  time
            newheader = ".".join(fields[0:5]) + "." + fields[9] + "." + str(bestIdx-5) + "_" + str(fields[bestIdx])
            

            subcount[fields[0]] = subcount.get(fields[0], 0) + 1
            patcount[pat] = fields[0]

            

            #outputfull.write(">" + newheader + "\n" + full[pat][header] + "\n")
            #outputv.write(newheader + "," + vseqdict[pat][header] + "\n")
            seqcount += 1
        else:
            subcount["C"] = subcount.get("C", 0) + 1
            patcount[pat] = "C"
            

            #outputfull.write(">"+ header + "\n" + full[pat][header] + "\n")
            #outputv.write(header + "," + vseqdict[pat][header] + "\n")
            seqcount += 1

    #outputfull.close()
    #outputv.close()
    
print(seqcount)
print(subcount)
print(len(patcount))
tp = np.array(tp)
print(tp[tp > 4])
