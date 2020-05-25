# script to retrieve all indels within the phylogenetic tree
# script will do the following:


# take in the root of the tree

# use a recursive function to iterate over every tree node 
    # if 
    # perform a comparison between each tree node sequence and its ancestor 

from ete2 import *
from seqUtils import * 
from glob import * 
from ancestors import *
from os.path import expanduser


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

def findChildren(node):
    #print(node)
    # finds the locations of all commas in the ancestral nodes
    commas = re.finditer(",",node)
    
    # saves all the end positions of the regex match objects
    idx = [x.end()-1 for x in commas]
    

    found = False
    i = 0

    # iterates through all comma positions until it finds the split point separating the two branches of the node 
    while found == False and i < len(idx):
        
        part1 = node[0:idx[i]]
        part2 = node[idx[i]+1:len(node)]

        left1 = len(re.findall("\(",part1))
        right1 = len(re.findall("\)", part1))

        # 1 and 2 refers to the two children of the node
        # left and right refer to the brackets that are being searched
        left2 = len(re.findall("\(",part2))
        right2 = len(re.findall("\)", part2))

        # left1 must be one higher than right1
        # right2 must be one higher than left2
        if left1 == right1+1 and left2+1 == right2:
            found = True
            #print(found)

            part1 = re.sub("^\({1}","",part1)
            part2 = re.sub(":[\d\.Ee-]+\)$","",part2)

            part1 = re.sub(":[\d\.Ee-]+$","",part1)
            return(part1, part2)
            
        i += 1
    return(None)

def extractIndels(tip, anc, vidx):
    vregion, pos = vidx

    iTemp = ''
    dTemp = ''
    #print(accno)
   
    #for collecting insertion and deletion sequences in each variable region 
    insertions = [[],[],[],[],[]]
    deletions = [[],[],[],[],[]]

    #for compiling the variable regions to ensure they are properly examined
    aseqs = ['','','','','']
    vseqs = ['','','','','']

    vSeq = {}
    aSeq = {}
   
    
    # case for retrieving the v-loops of the root node
    #if anc == "":
    #    for n, char in enumerate(tip):
    #        if vregion[n] != -1:
    #            vseqs[vregion[n]] += char
    #    vSeq[accno] = vseqs 
    #    return (vSeq, aSeq)
    # ai will count the number of nucleotides, skipping gaps 
    
    current = -1
    pidx = -1
    
    saved = 0
    # iterate through every character in the main sequence
    for n, schar in enumerate(tip):

        achar = anc[n]

        #following code only runs if inside a variable loop
        #if saved != vregion[n]:
            #print(n)
        
        # vregion is a list of numbers indicating what variable region the index belongs to 
        # [-1,-1,-1,-1,-1,-1,2,2,2,2,2,2,-1,-1,-1,-1]
        if vregion[n] != -1:
            
            # will run when entering new loop for first time 
            if vregion[n] != current:
                # start a count from zero 
                pidx = 0
                
                # update the current to the current loop number
                current = vregion[n]

            else:
                if achar != "-" or schar != "-":
                    pidx += 1

            print(pidx, end="")    

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
                        deletions[vregion[n]].append(dTemp+"-"+str(pidx))
                        dTemp = ''
            
                #nothing -- both gaps 
                else:
                    #clear iTemp and dTemp
                    if iTemp:
                        insertions[vregion[n]].append(iTemp+"-"+str(pidx))
                        iTemp = ''
                    if dTemp:
                        deletions[vregion[n]].append(dTemp+"-"+str(pidx))
                        dTemp = ''
                    
            elif achar != "-":
                #deletion in seq 1
                if schar == "-":
                    dTemp += achar
    
                    #clear iTemp
                    if iTemp:
                        insertions[vregion[n]].append(iTemp+"-"+str(pidx))
                        iTemp = ''
                    
                #nothing -- both have character
                else:
                    #clear iTemp and dTemp
                    if iTemp:
                        insertions[vregion[n]].append(iTemp+"-"+str(pidx))
                        iTemp = ''
                    if dTemp:
                        deletions[vregion[n]].append(dTemp+"-"+str(pidx))
                        dTemp = ''

    newvar = ['','','','','']
    newanc = ['','','','','']
    for n in range(5):
        for a,b in zip(vseqs[n], aseqs[n]):
            if a != "-" or b != "-":
                newvar[n] += a
                newanc[n] += b
    print(newanc)
        
    #SANITY CHECK 
    #ensures that the iterated sequences are the proper variable loops and that they are identical to the one found in the csv file 
    # for n, vr in enumerate(vseqs):
        # extract_tip = vr.replace("-","")
        # extract_anc = aseqs[n].replace("-","")
        
        # print(vr)
        # print(aseqs[n])
        # start, stop = boundaries[accno][n]
        # start = int(start)
        # stop = int(stop)

        # sliced = seqpairs[accno][0].replace("-","")[start:stop]
        #if extract_tip == extract_anc:
            #print(accno)
            #print(sliced)
            #print(extract_tip)
            #print(extract_anc)
            #print(csvSeq)'''

    return [insertions, deletions, newvar, newanc]
    
def addList(total, toAdd):
    for i in len(total):
        total[i].append(toAdd[i])


# dict[node label] = (insertions, deletions)


def isLeaf(header):
    res = re.search('^[^\(^\),\n]+$',header)
    return res != None

def treeIndelExtract(node, vidx, data):
    
    insertions = {}
    deletions = {}

    vseqs = {}
    aseqs = {}


    child1, child2 = findChildren(node)
    #print(child1)
    #print("---")
    #print(child2)
    #print(child1)
    #print(child2)
    #print(data[node])

    if isLeaf(child1):
        output = extractIndels(data[child1], data[node], vidx)
        sys.exit()
        insertions.update({child1:output[0]})
        deletions.update({child1:output[1]})
        vseqs.update({child1:output[2]})
        aseqs.update({child1:output[3]})
    else:       
        below = treeIndelExtract(child1, vidx, data)
        insertions.update(below[0])
        deletions.update(below[1])
        vseqs.update(below[2])
        aseqs.update(below[3])

        # 
        output = extractIndels(data[child1], data[node], vidx)
        insertions.update({child1:output[0]})
        deletions.update({child1:output[1]})
        vseqs.update({child1:output[2]})
        aseqs.update({child1:output[3]})

    if isLeaf(child2):
        output = extractIndels(data[child2], data[node], vidx)
        sys.exit()
        insertions.update({child2:output[0]})
        deletions.update({child2:output[1]})
        vseqs.update({child2:output[2]})
        aseqs.update({child2:output[3]})
    else:
        below = treeIndelExtract(child2, vidx, data)
        insertions.update(below[0])
        deletions.update(below[1])
        vseqs.update(below[2])
        aseqs.update(below[3])

        output = extractIndels(data[child2], data[node], vidx)
        insertions.update({child2:output[0]})
        deletions.update({child2:output[1]})
        vseqs.update({child2:output[2]})
        aseqs.update({child2:output[3]})
    
    #rootData = extractIndels(data[node], "", node, vregion)
    #print(rootData)
    return [insertions, deletions, vseqs, aseqs]
    
    #------------------
    

'''#print(folder)
if len(sys.argv) != 4:
    print("USAGE: python 9_1_ancestors.py [input Historian folder] [tree folder] [output folder]")
    quit() 
for i in range(len(sys.argv)):
    if not sys.argv[i].endswith("/"):
        sys.argv[i] += "/"
'''

sys.argv = ['',"8Historian/mcc/",'7_5_MCC/final/','9Indels/mcc/wholetree/']

folder = glob(sys.argv[1]+"*.fasta") #/rep/*.fasta
home = expanduser("~")
vpath = home+'/PycharmProjects/hiv-withinhost/3RegionSequences/variable/'
tpath = sys.argv[2]              
opath = sys.argv[3]

if not os.path.isdir(opath+"ins"):
    os.mkdir(opath+'ins')
if not os.path.isdir(opath+"del"):
    os.mkdir(opath+'del')


for f in folder:
    infile = open(f, "rU")
    filename = os.path.basename(f)
    print(filename)

    # retrieved the converted fasta 
    cdata = convert_fasta(infile)
    infile.close()

    # read the tree using the Tree using ETE3
    trefile = filename.split("_recon")[0] + ".tree"
    tree = Tree(tpath+trefile)
    
    # get the root label and sequence from the last entry of cdata
    root, rootseq = cdata[len(cdata)-1]
    

    # get the boundaries from the 'variable regions' csvfile
    # returns a list of 1) vregion sequences and 2) boundaries for those sequences
    csvfile = filename.split('-')[0] + ".csv"
    vregions, boundaries = getVRegions(vpath+csvfile)
        
    # previous = []  # sanity check to ensure that all vregion and pos lists are identical 
    
    # CALIBRATION -- use an arbitrary sequence to calibrate the vregion list and the pos list
    vregions = []
    pos = []
    header, seq = cdata[0]
    
    ai = 0
    for n, char in enumerate(seq):
        #retrieves a numeric value (0,1,2,3,4) to indicate which variable region the nucleotide is in, and -1 if outside of a vloop
        #pos counts your position WITHIN the current variable loop 
        vr, p = vrSwitch(ai, boundaries[header])
        vregions.append(vr)
        pos.append(p)
        if char != '-':
            ai += 1
    # previous = vidx

    
    infile = open(f, "rU")
    pdata = parse_fasta(infile)

    result = treeIndelExtract(root, (vregions, pos), pdata)
    
    #print(vregions)
    break 
    
    tsvout = filename.split("_recon")[0] + ".tsv"   #.tsv
    ioutput = open(opath+'ins/'+tsvout, 'w+')
    doutput = open(opath+'del/'+tsvout,'w+') 
    header = "header\tindel\tvloop\tvlen\ttip\tanc\n"
    ioutput.write(header)
    doutput.write(header)
    #fastaout = open(opath+csvout,"w+")
    #sys.exit()
    
    for accno in result[0]:
        i = result[0][accno]
        d = result[1][accno]
        v = result[2][accno]
        a = result[3][accno]

        for x, ins in enumerate(i):
            inslist = ":".join(ins)
            if inslist == "":
                inslist = ""
            newaccno = accno + "_" + str(x+1)
            ioutput.write("\t".join([newaccno, inslist, str(x+1), str(len(v[x])), v[x], a[x]]) + "\n")

            dellist = ":".join(d[x])
            if dellist == "":
                dellist = ""
            newaccno = accno + "_" + str(x+1)
            doutput.write("\t".join([newaccno, dellist, str(x+1), str(len(v[x])), v[x], a[x]]) + "\n")

    ioutput.close()
    doutput.close()
    






    
