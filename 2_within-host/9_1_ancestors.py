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
t = Tree("(((((((C.-.-.QJ.16.-.-_198:13.49219,((C.-.-.QJ.07.-.-_198:1.45601,C.-.-.QJ.13.-.-_198:1.45601):14.30421,C.-.-.QJ.02.-.-_198:15.76022):270.71697):20.55638,C.-.-.QJ.15.-.-_198:307.03358):30.77864,C.-.-.QJ.08.-.-_198:337.81221):22.21124,((C.-.-.QJ.10.-.-_344:120.04968,(C.-.-.QJ.11.-.-_344:65.53720,C.-.-.QJ.05R.-.-_377:98.53720):54.51249):29.65029,(((((C.-.-.QJ.07R.-.-_471:43.50810,C.-.-.QJ.14R.-.-_471:43.50810):10.44232,C.-.-.QJ.12R.-.-_471:53.95041):12.64329,C.-.-.QJ.16R.-.-_471:66.59370):54.59621,(C.-.-.QJ.01R.-.-_471:31.51815,C.-.-.QJ.08R.-.-_471:31.51815):89.67176):128.24699,(C.-.-.QJ.04.-.-_344:70.61993,C.-.-.QJ.13.-.-_344:70.61993):51.81697):27.26308):83.33848):434.35455,((C.-.-.QJ.06.-.-_198:91.95698,C.-.-.QJ.14.-.-_198:91.95698):329.99340,((C.-.-.QJ.11.-.-_071:9.35085,C.-.-.QJ.03.-.-_071:9.35085):284.19195,(((C.-.-.QJ.03.-.-_198:408.18850,((C.-.-.QJ.18R.-.-_008:32.94331,C.-.-.QJ.07.-.-_071:95.94331):131.80356,((C.-.-.QJ.14.-.-_071:128.00682,(C.-.-.QJ.04.-.-_071:125.53527,((C.-.-.QJ.11.-.-_198:237.71934,((C.-.-.QJ.19R.-.-_008:30.87580,C.-.-.QJ.05.-.-_071:93.87580):9.09507,C.-.-.QJ.15.-.-_071:102.97087):7.74847):1.35688,(C.-.-.QJ.02.-.-_344:377.75836,((((C.-.-.QJ.01.-.-_198:130.88335,(C.-.-.QJ.06.-.-_071:0.76698,C.-.-.QJ.01.-.-_071:0.76698):3.11637):40.76912,C.-.-.QJ.05.-.-_198:171.65247):27.28449,(C.-.-.QJ.04.-.-_198:96.93034,C.-.-.QJ.12.-.-_198:96.93034):102.00663):21.81617,(C.-.-.QJ.09.-.-_344:343.53552,(((C.-.-.QJ.04R.-.-_471:214.41347,((C.-.-.QJ.11R.-.-_471:72.25259,((C.-.-.QJ.02R.-.-_471:5.29524,C.-.-.QJ.17R.-.-_471:5.29524):35.56875,C.-.-.QJ.06R.-.-_471:40.86399):31.38860):26.89292,C.-.-.QJ.05R.-.-_471:99.14551):115.26796):100.64535,(C.-.-.QJ.09R.-.-_471:172.07219,((C.-.-.QJ.03.-.-_344:17.00900,(C.-.-.QJ.04R.-.-_377:31.51312,C.-.-.QJ.03R.-.-_471:125.51312):18.49588):9.66072,C.-.-.QJ.01.-.-_344:26.66971):18.40248):142.98663):93.05436,(C.-.-.QJ.15R.-.-_471:367.37135,((C.-.-.QJ.02R.-.-_377:85.63056,((C.-.-.QJ.03R.-.-_377:22.91475,C.-.-.QJ.06R.-.-_377:22.91475):26.07031,((C.-.-.QJ.09R.-.-_377:15.43959,C.-.-.QJ.07R.-.-_377:15.43959):21.98943,C.-.-.QJ.10R.-.-_377:37.42902):11.55603):36.64551):65.22438,C.-.-.QJ.08.-.-_344:117.85494):122.51640):40.74184):62.42234):23.21761):11.00522):7.31786):13.45905):2.47155):12.61969,C.-.-.QJ.17.-.-_071:140.62652):87.12035):53.44163):3.30199,((C.-.-.QJ.09.-.-_071:78.32393,(C.-.-.QJ.08.-.-_071:39.11806,(C.-.-.QJ.19.-.-_071:36.24508,(C.-.-.QJ.20.-.-_071:24.03750,C.-.-.QJ.10.-.-_071:24.03750):12.20758):2.87298):39.20587):69.01089,C.-.-.QJ.13.-.-_071:147.33482):137.15567):6.34994,C.-.-.QJ.02.-.-_071:290.84042):2.70237):1.40759):99.44262):228.64957,C.-.-.QJ.16.-.-_071:623.04257):94.07056,C.-.-.QJ.12.-.-_071:717.11313):0.00000;" )

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
    commas = re.finditer(",",node)
    idx = [x.end()-1 for x in commas]
    
    found = False
    i = 0
    while found == False and i < len(idx):
        
        part1 = node[0:idx[i]]
        part2 = node[idx[i]+1:len(node)]

        left1 = len(re.findall("\(",part1))
        right1 = len(re.findall("\)", part1))

        #print(left1)
        #print(right1)
        #print(part1)
        #print(part2)
        left2 = len(re.findall("\(",part2))
        right2 = len(re.findall("\)", part2))

        if left1 == right1+1 and left2+1 == right2:
            found = True
            #print(found)

            part1 = re.sub("^\({1}","",part1)
            part2 = re.sub(":[\d\.Ee-]+\)$","",part2)

            part1 = re.sub(":[\d\.Ee-]+$","",part1)
            return(part1, part2)
            
        i += 1
    return(None)

#print(getIndels(t))

def extractIndels(tip, anc, accno, vregions):
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
    for n, schar in enumerate(tip):

        achar = anc[n]

        #following code only runs if inside a variable loop
        #if saved != vregion[n]:
            #print(n)
        saved = vregion[n]

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
    #vSeq[accno] = newvar 
    #aSeq[accno] = newanc
    #print(vLen)
        
    #SANITY CHECK 
    #ensures that the iterated sequences are the proper variable loops and that they are identical to the one found in the csv file 
    '''for n, vr in enumerate(vseqs):
        extract_tip = vr.replace("-","")
        extract_anc = aseqs[n].replace("-","")
        
        # start, stop = boundaries[accno][n]
        # start = int(start)
        # stop = int(stop)

        # sliced = seqpairs[accno][0].replace("-","")[start:stop]
        if extract_tip == extract_anc:
            #print(accno)
            #print(sliced)
            print(extract_tip)
            print(extract_anc)
            #print(csvSeq)'''

    return ({accno: insertions},{accno:deletions})

def addList(total, toAdd):
    for i in len(total):
        total[i].append(toAdd[i])


# dict[node label] = (insertions, deletions)


def isLeaf(header):
    res = re.search('^[^\(^\),\n]+$',header)
    return res != None

def treeIndelExtract(node, vregion, data):
    
    insertions = {}
    deletions = {}

    child1, child2 = findChildren(node)
    #print(child1)
    #print("---")
    #print(child2)
    #print(child1)
    #print(child2)
    #print(data[node])

    if isLeaf(child1):
        output = extractIndels(data[child1], data[node], child1, vregion)
        insertions.update(output[0])
        deletions.update(output[1])
    else:       
        below = treeIndelExtract(child1, vregion, data)
        insertions.update(below[0])
        deletions.update(below[1])

        output = extractIndels(data[child1], data[node], child1, vregion)
        insertions.update(output[0])
        deletions.update(output[1])

    if isLeaf(child2):
        output = extractIndels(data[child2], data[node], child2, vregion)
        insertions.update(output[0])
        deletions.update(output[1])
    else:
        below = treeIndelExtract(child2, vregion, data)
        insertions.update(below[0])
        deletions.update(below[1])

        output = extractIndels(data[child2], data[node], child2, vregion)
        insertions.update(output[0])
        deletions.update(output[1])
    
    return (insertions, deletions)
    
    #------------------
    

#print(folder)

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/rep/*.fasta")
vpath = '/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/'
tpath = '/home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/final/'
opath = '/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/rep/wholetree/'

if not os.path.isdir(opath+"ins"):
    os.mkdir(opath+'ins')
if not os.path.isdir (opath+"del"):
    os.mkdir(opath+'del')

for f in folder:
    infile = open(f, "rU")
    filename = os.path.basename(f)
    print(filename)
    cdata = convert_fasta(infile)

    trefile = filename.split("_recon")[0] + ".tree"
    tree = Tree(tpath+trefile)
    csvfile = filename.split('-')[0] + ".csv"
    #root = t.get_tree_root()
    root = cdata[len(cdata)-1][0]
    #print(root)
    infile.close()

    #  LOADING VREGIONS 
    # ------------------
    # returns a list of 1) vregion sequences and 2) boundaries for those sequences
    vregions, boundaries = getVRegions(vpath+csvfile)
        
    # get alignment info (vregion and relative position lists) using a single iteration of the first sequence
    vregion = []
    pos = []
    #previous = []  # sanity check to ensure that all vregion and pos lists are identical 

    header, seq = cdata[0]
    ai = 0
    for n, char in enumerate(seq):
        #retrieves a numeric value (0,1,2,3,4) to indicate which variable region the nucleotide is in, and -1 if outside of a vloop
        # pos counts your position WITHIN the current variable loop 
        vr, p = vrSwitch(ai, boundaries[header])
        vregion.append(vr)
        pos.append(p)
        if char != '-':
            ai += 1
    previous = vregion
    #print(root)
    
    infile = open(f, "rU")
    pdata = parse_fasta(infile)

    result = treeIndelExtract(root, vregion, pdata)
    #print(result[0])
    tsvout = filename.split("_recon")[0] + ".tsv"
    ioutput = open(opath+'ins/'+tsvout, 'w+')
    doutput = open(opath+'del/'+tsvout,'w+') 
    header = "header\tV1\tV2\tV3\tV4\tV5\n"
    ioutput.write(header)
    doutput.write(header)
    for n, node in enumerate(result[0]):
        data = [":".join(vloop) for vloop in result[0][node]]
        #print(",".join([node, data, "\n"]))
        ioutput.write(node+"\t".join(data)+"\n")

    for n, node in enumerate(result[1]):
        data = [":".join(vloop) for vloop in result[1][node]]
        #print(",".join([node, data, "\n"]))
        doutput.write(node+"\t".join(data)+"\n")
    ioutput.close()
    doutput.close()

    #tips = [ x for x in data.keys() if re.search('^[^\(\):\n]+$',x) != None]

        






    