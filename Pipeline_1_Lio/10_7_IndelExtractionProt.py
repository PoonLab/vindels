from glob import glob
from seqUtils import *
from collections import defaultdict

files = glob("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_6_Pairwise_Prot/*.fastaz")


for x in files:
    alignment = open(x, 'r')

    data = parse_fasta2(alignment)

    filename = x.split("/")[-1].split(".")[0]

    #output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_7_NGlycs/" + filename +"_overlap.fasta", "w")
    output1 = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_7_NGlycs_total/" + filename +"_t.fasta", "w")
    output2 = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_8_NGlycs_i/" + filename + "_i.fasta", 'w')
    #Records EVERY N glyc sites in all pairs of cherry sequences
    ngCount = {}

    #Records all N-glyc sites that overlap with insertion sequences
    interfered = {}
    overlap = {}
    total = {}


    for i in data:

        header = ".".join(i.split(".")[0:3])


        #List of tuples containing insertion sequences (0) and their end positions (1)
        rInsert = []
        qInsert= []

        #Temporary strings used to store nucleotides from insertion sequences
        rTemp = ''
        qTemp = ''

        ri = 0
        qi = 0
        rindex = {}
        qindex = {}


        for n, rchar in enumerate(data[i][0]):
            qchar = data[i][1][n]

            #builds the index dictionaries
            if rchar != '-':
                rindex.update({ri:n})
                ri += 1
            if qchar != '-':
                qindex.update({qi:n})
                qi += 1

            # builds the Insert dictionaries: sequences (keys), end positions (values)
            #Case 1: no insertion, need to store and retrieve indel data
            if rchar != "-" and qchar != "-":
                if rTemp:
                    rInsert.append(tuple((rTemp,n-1)))
                    rTemp = ''
                if qTemp:
                    qInsert.append(tuple((qTemp,n-1)))
                    qTemp = ''

            #Case 2: query sequence has an insertion
            elif rchar == "-" and qchar != "-":
                qTemp = qTemp + qchar

                if rTemp:
                    rInsert.append(tuple((rTemp, n - 1)))
                    rTemp = ''
            #Case 3: reference sequence has an insertion
            elif rchar != "-" and qchar == "-":
                rTemp = rTemp + rchar
                if qTemp:
                    qInsert.append(tuple((qTemp,n-1)))
                    qTemp = ''
        if rTemp:
            rInsert.append(tuple((rTemp, n - 1)))
            rTemp = ''
        if qTemp:
            qInsert.append(tuple((qTemp, n - 1)))
            qTemp = ''

        '''
        #INSERTION OUTPUT: Used for generating the amino acid insertion sequences found in folder 10_2 
        #E.g. >HEADER \n >ref \n NGATSR,HGTYTH,KLGEST,
        if rInsert or qInsert:
            output.write(">" + i + "\n")
            output.write(">ref\n")
            if rInsert:
                for key in rInsert:
                    output.write(key[0] + ",")
            output.write("\n")
            output.write(">query\n")
            if qInsert:
                for key in qInsert:
                    output.write(key[0] + ",")
            output.write("\n")
        '''

        #--- N-linked glycosylation sites ---
        #--- 1: Reference sequences (first sequence in each pair)

        rGlycs = i.split(".")[3].split(",")
        qGlycs = i.split(".")[4].split(",")

        # For each nGlyc site, check every insertion in the same sequence for overlap; if there is overlap, add to the overlap dictionary
        for end in rGlycs:
            if end == '':
                continue
            end = int(end)

            # Uses minus 3 so that the start refers to the FIRST amino acid of the nGlyc
            start = end - 3

            #Translates glyc positions into aligned positions
            ngStart = rindex[start]
            ngEnd = rindex[end]

            #Add glyc site to the total dictionary
            if header in total.keys():
                total[header][0].append(ngEnd)
            else:
                total[header] = tuple(([ngEnd],[]))



            #Determine which N-glycs are disrupted
            for j in rInsert:
                inEnd = j[1]
                inStart = inEnd - (len(j[0]) - 1)   # minus 1 inserted so the start refers to the first amino acid




                #CHECK FOR OVERLAP
                if not ((ngStart > inEnd and ngEnd > inEnd) or (ngStart < inStart and ngEnd < inStart)):

                    if header in overlap.keys():
                        overlap[header][0].append(tuple((ngStart,ngEnd)))
                    else:
                        overlap[header] = tuple(([tuple((ngStart,ngEnd))],[]))

        #--- 2: Query sequences (second sequence in each pair)
        for end in qGlycs:
            if end == '':
                continue

            end = int(end)
            start = end - 3

            #Uses the original positions to determine the aligned position of the N-glycs

            ngStart = qindex[start]
            ngEnd = qindex[end]


            #for loading the total dictionary
            if header in total.keys():
                total[header][1].append(ngEnd)
            else:
                total[header] = tuple(([],[ngEnd]))


            #Determine which N-glycs are disrupted
            for j in qInsert:
                inEnd = j[1]
                inStart = inEnd - (len(j[0])-1)

                #CHECK FOR OVERLAP
                if not ((ngStart > inEnd and ngEnd > inEnd) or (ngStart < inStart and ngEnd < inStart)):
                    if header in overlap.keys():
                        overlap[header][1].append(tuple((ngStart,ngEnd)))
                    else:
                        overlap[header] = tuple(([],[tuple((ngStart,ngEnd))]))

    #FINAL OUTPUT : N-glycosylation sites
    #Loading the interfered dictionary




    for t in data:
        header = ".".join(t.split(".")[0:3])

        # only perform if this sequence pair has overlapping N glyc sites
        if header in overlap:
            rList = t.split(".")[3].split(",")
            qList = t.split(".")[4].split(",")
            rGlycs = []
            qGlycs = []

            ri = 0
            qi = 0
            rindex = {}
            qindex = {}

            #creates the alignment index again
            for n, rchar in enumerate(data[t][0]):
                qchar = data[t][1][n]

                if rchar != '-':
                    rindex.update({ri:n})
                    ri += 1
                if qchar != '-':
                    qindex.update({qi:n})
                    qi += 1


            if overlap[header][0]:
                for tpl in overlap[header][0]:
                    start = tpl[0]
                    end = tpl[1] + 1

                    inSeq = data[t][1][start:end]
                    beyond = data[t][1][end:].replace('-', '')

                    dashes = len(inSeq) - len(inSeq.replace('-', ''))

                    if inSeq[0] == "N" and '-' in inSeq and len(beyond) >= dashes :
                        i = 0
                        inSeq = list(inSeq.replace("-", ""))
                        while i < dashes and inSeq:
                            inSeq.append(beyond[i])

                            i += 1

                        inSeq = ''.join(inSeq)

                    test = re.search("N[^P][ST][^P]", inSeq)

                    if not test:
                        if header in interfered:
                            interfered[header][0].append(tpl)
                        else:
                            interfered[header] = tuple(([tpl],[]))

            if overlap[header][1]:
                for tpl in overlap[header][1]:
                    start = tpl[0]
                    end = tpl[1] + 1

                    inSeq = data[t][0][start:end]
                    beyond = data[t][0][end:].replace('-', '')

                    dashes = len(inSeq) - len(inSeq.replace('-', ''))

                    if inSeq[0] == "N" and '-' in inSeq and len(beyond) >= dashes:
                        i = 0
                        inSeq = list(inSeq.replace("-", ""))
                        while i < dashes and inSeq:
                            inSeq.append(beyond[i])

                            i += 1

                        inSeq = ''.join(inSeq)

                    test = re.search("N[^P][ST][^P]", inSeq)

                    if not test:
                        if header in interfered:
                            interfered[header][1].append(tpl)
                        else:
                            interfered[header] = tuple(([],[tpl]))


    print(x)
    print(total)
    print(overlap)
    print(interfered)

    for m in total:
        output1.write(str(m) + ":" + ",".join(map(str, total[m][0])) + ":" + ",".join(map(str, total[m][1])) + "\n")


    '''for n in overlap:
        out1 = []
        out2 = []
        for x in overlap[n][0]:
            out1.append(list(x)[1])
        for y in overlap[n][1]:
            out2.append(list(y)[1])
        #print(str(n) + "." + ",".join(map(str,out1)) + "." + ",".join(map(str,out2)) + "\n")'''

    print("INTERFERED------------")
    for n in interfered:
        out1 = []
        out2 = []
        for x in interfered[n][0]:
            out1.append(list(x)[1])
        for y in interfered[n][1]:
            out2.append(list(y)[1])
        output2.write(str(n) + ":" + ",".join(map(str, out1)) + ":" + ",".join(map(str, out2)) + "\n")



'''for w in overlap:
    output.write(str(w) + "." + ",".join(map(str,overlap[w][0])) + "." + ",".join(map(str,overlap[w][1])) + "\n")
for v in total:
    output2.write(str(v) + "." + ",".join(map(str,total[v][0])) + "." + ",".join(map(str,total[v][1])) + "\n")'''




'''#sets of overlapped glyc sites
    overlapSet1 = set(overlap[header][0])
    overlapSet2 = set(overlap[header][1])

    #load the rGlycs and qGlycs
    if not(len(rList) == 1 and rList[0] == ''):
        for o in rList:
            o = int(o)
            rGlycs.append(tuple((rindex[o-3], rindex[o])))
    if not(len(qList) == 1 and qList[0] == ''):
        for p in qList:
            p = int(p)
            qGlycs.append(tuple((qindex[p-3], qindex[p])))

    print(rGlycs)
    print(qGlycs)


    ngSet1 = set(rGlycs)
    ngSet2 = set(qGlycs)

    #convert both lists to integers
    if len(rGlycs) == 1 and rGlycs[0] == '':
        ngSet1 = set()
    else:
        ngSet1 = set(map(int, rGlycs))

    if len(qGlycs) == 1 and qGlycs[0] == '':
        ngSet2 = set()
    else:
        ngSet2 = set(map(int, qGlycs))'''






