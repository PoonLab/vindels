from glob import * 
from sys import *
import os 
import csv
from seqUtils import * 



def getTerminals(anFile):
    terminals = {}
    with open(anFile) as handle:
        data = parse_fasta(handle)

    # pos1 and pos2 check for an ancestor that directly links to a tip 
    pos1 = re.compile("^\([^(:,)]*:[^(:,)]*,")
    pos2 = re.compile(",[^(:,)]*:[^(:,)]*\)$")

    for entry in data.keys():
        
        check1 = pos1.search(entry)
        check2 = pos2.search(entry)

        # FORMAT : trm1[header] = [sequence, ancestor]
        if check1 != None:
  
            # header of the first sequence only 
            header = entry.split(",")[0].strip(">(,").split(":")[0]

            # retrieve the tip seq 
            tip = data[header]

            # retrieve the anc seq and make it upper case 
            anc = data[entry].upper()

            # retrieve just the relevant nts from the alignment
            '''new_tip = []
            new_anc = []
            for t, a in zip(tip,anc):
                if t != "-" or a != "-":
                    new_tip.append(t)
                    new_anc.append(a)
            new_tip = "".join(new_tip)
            new_anc = "".join(new_anc)'''
            new_tip = tip
            new_anc = anc
            terminals[header] = [new_tip, new_anc]

        if check2 != None:
 
            header = entry.split(",")[-1].strip(",)").split(":")[0]
            
            tip = data[header]
            anc = data[entry].upper()
            '''new_tip = []
            new_anc = []
            for t, a in zip(tip,anc):
                if t != "-" or a != "-":
                    new_tip.append(t)
                    new_anc.append(a)
            new_tip = "".join(new_tip)
            new_anc = "".join(new_anc)'''

            new_tip = tip
            new_anc = anc
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
    # vregions[accession number] = list(list(V1seq, V2seq, V3seq, V4seq, V5seq)...)
    # positions[accesion number] = list([[start, stop],[start, stop],[start,stop],[start,stop],[start,stop]]...)
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
def extractIndels(seqpairs, vSeqFile):
    count = 0

    # returns a list of 1) vregion sequences and 2) boundaries for those sequences
    vregions, boundaries = getVRegions(vSeqFile)

    iDict = {}
    dDict = {}
    vSeq = {}
    aSeq = {}

    # get alignment info (vregion and relative position lists) using a single iteration of the first sequence
    vregion = []
    pos = []
    #previous = []  # sanity check to ensure that all vregion and pos lists are identical 
    for accno in seqpairs.keys():
        ai = 0
        for n, char in enumerate(seqpairs[accno][0]):
            #retrieves a numeric value (0,1,2,3,4) to indicate which variable region the nucleotide is in, and -1 if outside of a vloop
            # pos counts your position WITHIN the current variable loop 
            vr, p = vrSwitch(ai, boundaries[accno])
            vregion.append(vr)
            pos.append(p)
            if char != '-':
                ai += 1
        break

    for accno in seqpairs.keys():
        iTemp = ''
        dTemp = ''
        #print(accno)

        #for collecting insertion and deletion sequences in each variable region 
        insertions = [[],[],[],[],[]]
        deletions = [[],[],[],[],[]]

        #for compiling the variable regions to ensure they are properly examined
        aseqs = ['','','','','']
        vseqs = ['','','','','']
        
        # ai will count the number of nucleotides, skipping gaps 
        
        saved = 0
        # iterate through every character in the main sequence 
        for n, schar in enumerate(seqpairs[accno][0]):
            #print schar,
            # also retrieve the chars in the ancestral sequence 
            achar = seqpairs[accno][1][n]
            #vregion, pos = vrSwitch(ai, boundaries[accno]) 
            #this counts up the number of nucleotides found in the seq of interest
            
            # sanity check to ensure the boundaries are being found properly
            #if saved != vregion[n]:
                #print(n)
            #saved = vregion[n]
            
            
            #following code only runs if inside a variable loop
            if vregion[n] != -1:
                #print(n)
                
                #sanity check to ensure the code is covering the variable regions 
                aseqs[vregion[n]] += achar
                vseqs[vregion[n]] += schar

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
                            deletions[vregion[n]].append(dTemp+"-"+str(pos[n]))
                            dTemp = ''
                
                    #nothing -- both gaps 
                    else:
                        #clear iTemp and dTemp
                        if iTemp:
                            insertions[vregion[n]].append(iTemp+"-"+str(pos[n]))
                            iTemp = ''
                        if dTemp:
                            deletions[vregion[n]].append(dTemp+"-"+str(pos[n]))
                            dTemp = ''
                        
                elif achar != "-":
                    #deletion in seq 1
                    if schar == "-":
                        dTemp += achar
        
                        #clear iTemp
                        if iTemp:
                            insertions[vregion[n]].append(iTemp+"-"+str(pos[n]))
                            iTemp = ''
                        
                    #nothing -- both have character
                    else:
                        #clear iTemp and dTemp
                        if iTemp:
                            insertions[vregion[n]].append(iTemp+"-"+str(pos[n]))
                            iTemp = ''
                        if dTemp:
                            deletions[vregion[n]].append(dTemp+"-"+str(pos[n]))
                            dTemp = ''
        newvar = ['','','','','']
        newanc = ['','','','','']
        for n in range(5):
            for a,b in zip(vseqs[n], aseqs[n]):
                if a != "-" or b != "-":
                    newvar[n] += a
                    newanc[n] += b
        iDict[accno] = insertions
        dDict[accno] = deletions
        vSeq[accno] = newvar 
        aSeq[accno] = newanc
        #print(vLen)
        
        #SANITY CHECK 
        #ensures that the iterated sequences are the proper variable loops and that they are identical to the one found in the csv file 
        '''for n, vr in enumerate(vseqs):
            extracted = vr.replace("-","")
            #extract_anc = aseqs[n].replace("-","")
            csvSeq = vregions[accno][n] 
            
            print(extracted)
            print(csvSeq)

            start, stop = boundaries[accno][n]
            start = int(start)
            stop = int(stop)
            sliced = seqpairs[accno][0].replace("-","")[start:stop]
            if extracted == sliced and extracted == csvSeq and sliced == csvSeq:
                #print(accno)
                #print(sliced)
                print(extracted)
                print(extract_anc)
                #print(csvSeq)
                count +=1'''

    return iDict, dDict, vSeq, aSeq

def main():
    '''hFolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/mcc/*.fasta")
    print(hFolder[0])
    result = getTerminals(hFolder[0])
    for r in result:
        print(r)
        print(result[r][0][:250])
        print(result[r][1][:250])'''

    #if len(sys.argv) != 3:
        #print("USAGE: python 9_ancestors.py [input Historian folder] [output folder]")
        #quit()
    #for i in range(len(sys.argv)):
        #if not sys.argv[i].endswith("/"):
            #sys.argv[i] += "/"

    sample = '''----------atgagagtgatggggatacagaggaattgtcaacgatggt\
ggatatggggcatcttaggcttttggatgatatatagtgggatggggcaa\
ttgtgggtcacagtctattatggggtacctgtatggaaggaagcaaaaac\
tactctattttgtgcatcagatgctaaagcatattcccaagaagcacata\
atatctgggctacacatgcctgcgtacctacagaccccaatccacaagaa\
atggttttagaaaatgtaacagaaaattttaacatgtggaaaaatgacat\
ggtagatcagatgcatgaggatataatcagtttatgggatcaaagcctaa\
agccatgtgtaaagttgaccccactctgtgtcactttaaactgtaccact\
gttaatgttaatgttaccc------------------agccccttaatga\
tagtgatag---ttacatgagaaactgctctttcaatatgaccacagaac\
tacaagataaaaagaaggaagagagggcattgtttcatagacttgatcta\
gtaccacttaatgag------------------aatagtagcc-------\
-----------a------------------gtatatattaacaaattgta\
atacctcaaccatagcacaggcctgtccaaagatcgcgtttgatccgatt\
cctatacattattgtgctccagctggttatgcgattctaaaatgtaataa\
taaggcattcaatggatcaggaccatgcactaatgtcagcacagtacaat\
gtacacatggaattaagccagtggtatcaactcagctactgttaaatggt\
agcctagcagaagagataataatcagatctgaaaatctgacaaacaatgt\
caaaacaataatagtacaccttaatgaatctgtagaaattttgtgtacaa\
gacccaacaataatacaagaaaaagtatgaggataggaccaggacaaaca\
ttctatgcaacaggagaaataataggagatataagacaagcacattgtaa\
cattactgaaggaaaatggaataaaactttagaaagaataaggaacaaat\
tagcagagtacttccctaataaaacaatacaatttcaaagcgcctcagga\
ggggacctagaaattacatcacatagctttaattgtggaggagaattttt\
ctattgcaacacatcaagattgtttaatcatacatct-------------\
--------a---------------------atgagacatcaaacaaaaac\
atcacactccaatgtaggataaaacaaattataaacatgtggcagggggt\
aggacgagcaatatatgcccctcccattgcaggaaacataacatgtaact\
caaatatcacaggactactattgacacgtgatggaggggagagtaatagt\
agtcatcagacagagatattcagacctgcaggaggaaatatgagggacaa\
ttggagaagtgaattatataaatataaagtggtagaaattaaaccattag\
gaatagcacccactaaggcgaaaaggagagtggtggagagagaaaaaaga'''

    #print(len(sample))
    
    # folder of historian outputs 
    hFolder = glob(sys.argv[1]+"*.fasta")
    # output path for indel csv files
    outpath = sys.argv[2]
    if not os.path.isdir(outpath+"ins/"):
        os.mkdir(outpath+"ins/")
    if not os.path.isdir(outpath+"del/"):
        os.mkdir(outpath+"del/")
    
    # path to the variable region sequences 
    vPath = '/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/'

    #totalseqs = 0
    for infile in hFolder:
        #print(infile)
        filename = os.path.basename(infile)
        #if len(sys.argv) == 4:
        #replicate = filename.split("_")[0][-1]
        #create names for both the csv file and the output recon file 
        csvfile = filename.split('-')[0] + ".csv"     #101827.csv
        reconfile = filename.split("_recon")[0] + ".csv"   #101827-a_15.csv

        if os.path.isfile(vPath+csvfile):
            ins_out = open(outpath+"ins/"+reconfile,'w')
            del_out = open(outpath+"del/"+reconfile,'w')
            print(infile)
            terminals = getTerminals(infile)
            iDict, dDict, vSeq, aSeq = extractIndels(terminals, vPath+csvfile)
            
            #print(len(iDict))
            #print(len(vSeq))
            ins_out.write("header,ins,Vloop,Vlen,Vseq,anc\n")
            del_out.write("header,del,Vloop,Vlen,Vseq,anc\n")

            for accno in iDict:
                #totalseqs+=1
                insertions = iDict[accno] # [ [], [], [], [], [] ]
                deletions = dDict[accno] # [ [], [], [], [], [] ]
                vsequences = vSeq[accno]   # [ V1-len, V2-len, V3-len, V4-len, V5-len]
                asequences = aSeq[accno]
                # j/k count from 0 to 4 for each variable loop 
                for j, ins in enumerate(insertions):
                    insList = ":".join(ins)
                    if insList == "":
                        insList = ""
                    newaccno = accno + "_"  + str(j+1) #str(replicate)
                    ins_out.write(",".join([newaccno, insList, str(j+1), str(len(vsequences[j])), vsequences[j], asequences[j]])+"\n")
                for k, dl in enumerate(deletions):
                    delList = ":".join(dl)
                    if delList == "":
                        delList = ""
                    newaccno = accno + "_" + str(k+1)
                    del_out.write(",".join([newaccno, delList, str(k+1), str(len(vsequences[k])), vsequences[k], asequences[k]])+"\n") 
    #print(totalseqs)


if __name__ == '__main__': 
    main()

