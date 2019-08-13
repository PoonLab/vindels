from glob import * 
import sys 
import os 
import csv
from seqUtils import * 



def getTerminals(anFile):
    terminals = {}
    
    with open(anFile) as handle:
        data = parse_fasta(handle)

    # pos1 and pos2 check for an ancestor that directly links to a tip 
    pos1 = re.compile("^>\([^(:,)]*:[^(:,)]*,")
    pos2 = re.compile(",[^(:,)]*:[^(:,)]*\)$")
    
    for entry in data.keys():
        #print(entry)
        #check = cherry.match(entry)
        check1 = pos1.search(entry)
        check2 = pos2.search(entry)
        #print("HELLO")
        # FORMAT : terminals[header] = [sequence, ancestor]
        if check1 != None:
            header = entry.split(",")[0].strip(">(,").split(":")[0]
            #print(check1)
            tip = data[header]
            anc = data[entry].upper()
            new_tip = []
            new_anc = []
            for t, a in zip(tip,anc):
                if t != "-" or a != "-":
                    new_tip.append(t)
                    new_anc.append(a)
            new_tip = "".join(new_tip)
            new_anc = "".join(new_anc)

            terminals[header] = [new_tip, new_anc]

        if check2 != None:
            header = entry.split(",")[-1].strip(",)").split(":")[0]
            
            tip = data[header]
            anc = data[entry].upper()
            new_tip = []
            new_anc = []
            for t, a in zip(tip,anc):
                if t != "-" or a != "-":
                    new_tip.append(t)
                    new_anc.append(a)
            new_tip = "".join(new_tip)
            new_anc = "".join(new_anc)

            terminals[header] = [new_tip, new_anc]

    # RETURNS:
    # terminals[accession number] = [sequence, ancestral sequence]
    return terminals

# Depends on vrSwitch
# vSeqFile format: header, V1, start, stop, V2, start, stop, V3, start, stop, V4, start, stop, V5, start, stop 
def getVRegions(vSeqFile):

    with open(vSeqFile) as handle:
        #used to skip the header line
        handle.readline()
        data = handle.readlines()
    
    vregions = {}
    positions = {}
    for line in data:
        fields = line.strip("\n").split(",")

        # vregions[accession number] = list(V1seq, V2seq, V3seq, V4seq, V5seq)
        vregions.update({fields[0]:[]})
        positions.update({fields[0]:[]})
        for vr in range(5):
            vregions[fields[0]].append(fields[(vr*3)+1])
            positions[fields[0]].append([fields[(vr*3)+2], fields[(vr*3)+3]])

    # RETURNS
    # vregions[accession number] = list(V1seq, V2seq, V3seq, V4seq, V5seq)
    # positions[accesion number] = list([start, stop],[start, stop],[start,stop],[start,stop],[start,stop])
    return vregions, positions
        
def vrSwitch(position, boundaries):
    for i in range(5):
        start, stop = boundaries[i]
        start = int(start)
        stop = int(stop)
        if position in range(start,stop):
            return i, position-start
    return -1, -1
    
            
# Dependencies: 1) getTerminals   2) vrSwitch   3) getVRegions  
def extractIndels(anFile, vSeqFile):
    count = 0

    # {"accno":[seq1,anseq]}
    terminals = getTerminals(anFile)

    # returns a list of 1) vregion sequences and 2) boundaries for those sequences
    vregions, boundaries = getVRegions(vSeqFile)
    iDict = {}
    dDict = {}
    vSeq = {}
    for accno in terminals.keys():
        iTemp = ''
        dTemp = ''

        #for collecting insertion and deletion sequences in each variable region 
        insertions = [[],[],[],[],[]]
        deletions = [[],[],[],[],[]]

        #for compiling the variable regions to ensure they are properly examined
        seqs = ['','','','','']
        
        # ai will count the number of nucleotides, skipping gaps 
        ai = 0
        
        # iterate through every character in the main sequence 
        for n, schar in enumerate(terminals[accno][0]):

            # also retrieve the chars in the ancestral sequence 
            achar = terminals[accno][1][n]

            #retrieves a numeric value (0,1,2,3,4) to indicate which variable region the nucleotide is in, and -1 if outside of a vloop
            vregion, pos = vrSwitch(ai, boundaries[accno]) 
            
            #this counts up the number of nucleotides found in the seq of interest
            if schar != '-':
                ai += 1

            #following code only runs if inside a variable loop
            if vregion != -1:
                
                #sanity check to ensure the code is covering the variable regions 
                seqs[vregion] += schar

                #1 : ancestral gap, sequence gap  = nothing 
                #2 : ancestral gap, sequence char = insertion 
                #3 : ancestral char, sequence gap  = deletion 
                #4 : ancestral char, sequence char  = nothing 
                
                if achar == "-":
                    #insertion in seq 1 
                    if schar != '-':
                        iTemp += schar
                        
                    #clear the dTemp 
                        if dTemp:
                            deletions[vregion].append(dTemp+"-"+str(pos))
                            dTemp = ''
                
                    #nothing -- both gaps 
                    else:
                        #clear iTemp and dTemp
                        if iTemp:
                            insertions[vregion].append(iTemp+"-"+str(pos))
                            iTemp = ''
                        if dTemp:
                            deletions[vregion].append(dTemp+"-"+str(pos))
                            dTemp = ''
                        
                elif achar != "-":
                    #deletion in seq 1
                    if schar == "-":
                        dTemp += achar
        
                        #clear iTemp
                        if iTemp:
                            insertions[vregion].append(iTemp+"-"+str(pos))
                            iTemp = ''
                        
                    #nothing -- both have character
                    else:
                        #clear iTemp and dTemp
                        if iTemp:
                            insertions[vregion].append(iTemp+"-"+str(pos))
                            iTemp = ''
                        if dTemp:
                            deletions[vregion].append(dTemp+"-"+str(pos))
                            dTemp = ''

        for n in range(5):
            seqs[n] = seqs[n].replace("-","")
        iDict[accno] = insertions
        dDict[accno] = deletions
        vSeq[accno] = seqs 
        #print(vLen)
        '''
        #SANITY CHECK 
        #ensures that the iterated sequences are the proper variable loops and that they are identical to the one found in the csv file 
        for n, vr in enumerate(vseqs):
            extracted = vr.replace("-","")
            
            csvSeq = vregions[accnos[i]][n] 
            
            start, stop = boundaries[accnos[i]][n]
            start = int(start)
            stop = int(stop)
            sliced = terminals[header][i].replace("-","")[start:stop]
            if extracted == sliced and extracted == csvSeq and sliced == csvSeq:
                print(accnos[i])
                print(sliced)
                print(extracted)
                print(csvSeq)
                count +=1'''

    return iDict, dDict, vSeq

def main():
    '''hFolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/mcc/*.fasta")
    print(hFolder[0])
    result = getTerminals(hFolder[0])
    for r in result:
        print(r)
        print(result[r][0][:250])
        print(result[r][1][:250])'''

    if not sys.argv[1].endswith("/"):
        sys.argv[1] += "/"
    hFolder = glob(sys.argv[1]+"*.fasta")
    
    vPath = '/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/'

    totalseqs = 0
    for infile in hFolder:
        print(infile)
        filename = os.path.basename(infile)

        #create names for both the csv file and the output recon file 
        csvfile = filename.split('-')[0] + ".csv"     #101827.csv
        reconfile = filename.split("_recon")[0] + ".csv"   #101827-a_15.csv

        if os.path.isfile(vPath+csvfile):
            ins_out = open("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/new-VN/ins/"+reconfile,'w')
            del_out = open("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/new-VN/del/"+reconfile,'w')
            iDict, dDict, vSeq = extractIndels(infile, vPath+csvfile)

            ins_out.write("Accno,Ins,Vloop,Vlen,Seq\n")
            del_out.write("Accno,Del,Vloop,Vlen,Seq\n")

            for accno in iDict:      
                insertions = iDict[accno] # [ [], [], [], [], [] ]
                deletions = dDict[accno] # [ [], [], [], [], [] ]
                vsequences = vSeq[accno]   # [ V1-len, V2-len, V3-len, V4-len, V5-len]
                
            # j/k count from 0 to 4 for each variable loop 
                for j, ins in enumerate(insertions):
                    insList = ":".join(ins)
                    if insList == "":
                        insList = ""
                    ins_out.write(",".join([accno,insList,str(j+1),str(len(vsequences[j])),vsequences[j]])+"\n")
                for k, dl in enumerate(deletions):
                    delList = ":".join(dl)
                    if delList == "":
                        delList = ""
                    del_out.write(",".join([accno,delList,str(k+1), str(len(vsequences[k])), vsequences[k]])+"\n")
        

    print(totalseqs)


if __name__ == '__main__': 
    main()
