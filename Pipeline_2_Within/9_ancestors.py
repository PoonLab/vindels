from glob import * 
import sys 
import os 
import csv
from seqUtils import * 

def opposite(num):
    if num == 0:
        return 1
    else:
        return 0

def getNodes(anFile):
    nodes = {}
    
    with open(anFile) as handle:
        data = parse_fasta(handle)

    cherry = re.compile("\([^()]*:[^()]*,[^()]*\)")
    
    for header in data.keys():
        #print(header)
        check = cherry.match(header)

        if check != None:
            s1, s2  = header.strip("()").split(",")
            s1 = s1.split(":")[0]
            s2 = s2.split(":")[0]
            label =  s1 + "." + s2

            nodes[label] = [data[s1],data[s2],data[header]]
    return nodes

def getVRegions(vSeqFile):

    with open(vSeqFile) as handle:
        #used to skip the header line
        handle.readline()
        data = handle.readlines()
    
    vregions = {}
    positions = {}
    for line in data:
        fields = line.strip("\n").split(",")
        vregions.update({fields[0]:[]})
        positions.update({fields[0]:[]})
        for vr in range(5):
            vregions[fields[0]].append(fields[(vr*3)+1])
            positions[fields[0]].append([fields[(vr*3)+2], fields[(vr*3)+3]])

    return vregions, positions
        
def vrSwitch(position, boundaries):

    for i in range(5):
        start, stop = boundaries[i]
        start = int(start)
        stop = int(stop)
        if position in range(start,stop):
            return i 
    return -1
    
            
            

def extractIndels(anFile, vSeqFile):
    count = 0

    # {"accno1.accno2":[seq1,seq2,anseq]}
    nodes = getNodes(anFile)
    vregions, positions = getVRegions(vSeqFile)
    iDict = {}
    dDict = {}
    for header in nodes.keys():
        accnos = header.split('.')
        iDict[header] = []
        dDict[header] = []

        #repeat once for each sequence in the cherry 
        for i in range(2):
            
            iTemp = ''
            dTemp = ''
            insertions = [[],[],[],[],[]]
            deletions = [[],[],[],[],[]]
            vseqs = ['','','','','']
            ai = 0
            for n, Achar in enumerate(nodes[header][2]):
                Schar = nodes[header][i][n]
                Ochar = nodes[header][opposite(i)][n]

                #retrieves a numeric value (0,1,2,3,4) to indicate which variable region the nucleotide is in, and -1 if outside of a vloop
                vregion = vrSwitch(ai, positions[accnos[i]])
                
                #this counts up the number of nucleotides found in the SOI (seq of interest)
                if Schar != '-':
                    ai += 1

                #following code only run if inside a variable loop
                if vregion != -1:
                    
                    #sanity check to ensure the code is covering the 
                    vseqs[vregion] += Schar
                    #1 : ancestral gap, sequence gap  = nothing 
                    #2 : ancestral gap, sequence char = insertion 
                    #3 : ancestral char, sequence gap  = deletion 
                    #4 : ancestral char, sequence char  = nothing 
                    if Achar == "-":
                        #insertion in seq 1 
                        if Schar != '-':
                            iTemp += Schar
                            #print(Achar + "" + Schar+"" + Ochar)
                        #clear the dTemp 
                            if dTemp:
                                deletions[vregion].append(dTemp+"-"+str(n-1))
                                dTemp = ''
                    
                        #nothing -- both gaps 
                        else:
                            #clear iTemp and dTemp
                            if iTemp:
                                insertions[vregion].append(iTemp+"-"+str(n-1))
                                iTemp = ''
                            if dTemp:
                                deletions[vregion].append(dTemp+"-"+str(n-1))
                                dTemp = ''
                            
                    elif Achar == "*":
                        #deletion in seq 1
                        if Schar == "-":
                            dTemp += Ochar
                            #print(Achar + "" + Schar+"" + Ochar)
                            #clear iTemp
                            if iTemp:
                                insertions[vregion].append(iTemp+"-"+str(n-1))
                                iTemp = ''
                            
                        #nothing -- both have character
                        else:
                            #clear iTemp and dTemp
                            if iTemp:
                                insertions[vregion].append(iTemp+"-"+str(n-1))
                                iTemp = ''
                            if dTemp:
                                deletions[vregion].append(dTemp+"-"+str(n-1))
                                dTemp = ''
            iDict[header].append(insertions)
            dDict[header].append(deletions)
            '''
            #SANITY CHECK 
            #ensures that the iterated sequences are the proper variable loops and that they are identical to the one found in the csv file 
            for n, vr in enumerate(vseqs):
                extracted = vr.replace("-","")
                
                csvSeq = vregions[accnos[i]][n] 
                
                start, stop = positions[accnos[i]][n]
                start = int(start)
                stop = int(stop)
                sliced = nodes[header][i].replace("-","")[start:stop]
                if extracted == sliced and extracted == csvSeq and sliced == csvSeq:
                    print(accnos[i])
                    print(sliced)
                    print(extracted)
                    print(csvSeq)
                    count +=1'''

    return iDict, dDict

def main():
    hFolder = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/*.fasta')
    vPath = '/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/'
    
    for infile in hFolder:
        filename = os.path.basename(infile)
        #create names for both the csv file and the output recon file 
        csvfile = filename.split('-')[0] + ".csv"
        reconfile = filename.split("_")[0] + ".csv"
        ins_out = open("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/insertions/"+reconfile,'w')
        del_out = open("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/deletions/"+reconfile,'w')
        iDict, dDict = extractIndels(infile, vPath+csvfile)

        ins_out.write("Accno\tDate\tIns\tVloop\n")
        del_out.write("Accno\tDate\tDel\tVloop\n")

        for key in iDict:
            header = key.split(".")
            
            for n in range(2):
                accno, date = header[n].split("_")
                
                insertions = iDict[key][n]
                deletions = dDict[key][n]
                for j, ins in enumerate(insertions):
                    insList = ",".join(ins)
                    if insList == "":
                        insList = ""
                    ins_out.write("\t".join([accno,date,insList,str(j+1)])+"\n")
                for k, dl in enumerate(deletions):
                    delList = ",".join(dl)
                    if delList == "":
                        delList = ""
                    del_out.write("\t".join([accno,date,delList,str(j+1)])+"\n")




if __name__ == '__main__': 
    main()

