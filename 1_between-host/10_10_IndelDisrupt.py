from glob import glob
from seqUtils import *
from collections import defaultdict

files = glob("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_6_Pairwise_Prot/*.fastaz")


for x in files:
    print(x)
    alignment = open(x, 'r')

    data = parse_fasta2(alignment)

    filename = x.split("/")[-1].split(".")[0]

    #output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_7_NGlycs/" + filename +"_overlap.fasta", "w")
    #totOut = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_9_Indels_total/" + filename +"_t.fasta", "w")
    #intOut = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/10_10_Indels_i/" + filename + "_i.fasta", 'w')

    #Records EVERY N glyc sites in all pairs of cherry sequences
    ngCount = {}

    #Records all N-glyc sites that overlap with insertion sequences
    interfered = {}
    overlap = {}
    total = {}


    for i in data:

        header = ".".join(i.split(".")[0:3])

        #print(header)
        #print(data[i][0])
        #print(data[i][1])

        #List of tuples containing insertion sequences (0) and their end positions (1)
        rqInsert = tuple(([],[]))

        #Temporary strings used to store nucleotides from insertion sequences
        rTemp = ''
        qTemp = ''

        ri = 0
        qi = 0
        rqindex = tuple(({},{}))


        #FINDS ALL INDELS
        for n, rchar in enumerate(data[i][0]):
            qchar = data[i][1][n]

            #builds the index dictionaries
            if rchar != '-':
                rqindex[0].update({ri:n})
                ri += 1
            if qchar != '-':
                rqindex[1].update({qi:n})
                qi += 1


            # builds the Insert dictionaries: sequences (keys), end positions (values)
            #Case 1: no insertion, need to store and retrieve indel data
            if rchar != "-" and qchar != "-":
                if rTemp:
                    rqInsert[0].append(tuple((rTemp,n-1)))
                    rTemp = ''
                if qTemp:
                    rqInsert[1].append(tuple((qTemp,n-1)))
                    qTemp = ''

            #Case 2: query sequence has an insertion
            elif rchar == "-" and qchar != "-":
                qTemp = qTemp + qchar

                if rTemp:
                    rqInsert[0].append(tuple((rTemp, n - 1)))
                    rTemp = ''
            #Case 3: reference sequence has an insertion
            elif rchar != "-" and qchar == "-":
                rTemp = rTemp + rchar
                if qTemp:
                    rqInsert[1].append(tuple((qTemp,n-1)))
                    qTemp = ''

        #handles the END case
        if rTemp:
            rqInsert[0].append(tuple((rTemp, n - 1)))
            rTemp = ''
        if qTemp:
            rqInsert[1].append(tuple((qTemp, n - 1)))
            qTemp = ''
        #print(rindex)
        #print(qindex)

        #--- N-linked glycosylation sites ---
        #--- 1: Reference sequences (first sequence in each pair)

        nGlycs = tuple((i.split(".")[3].split(","), i.split(".")[4].split(",")))


        # For each nGlyc site, check every insertion in the same sequence for overlap; if there is overlap, add to the overlap dictionary
        for a in range(2):
            for info in rqInsert[a]:
                inEnd = info[1]
                inStart = inEnd - (len(info[0]) - 1)
                # Uses minus 3 so that the start refers to the FIRST amino acid of the nGlyc


                #Translates glyc positions into aligned positions


                #Add indel site to the total dictionary
                if header in total.keys():
                    total[header][a].append(inEnd)
                elif a == 0:
                    total[header] = tuple(([inEnd],[]))
                else:
                    total[header] = tuple(([], [inEnd]))


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
                            overlap[header] = tuple(([], [tuple((istr, gstr))]))


    #FINAL OUTPUT : N-glycosylation sites
    #Loading the interfered dictionary


    #print(total)
    #print(overlap)


    #verification step

    for t in data:
        header = ".".join(t.split(".")[0:3])

        # only perform if this sequence pair has overlapping N glyc sites
        if header in overlap:
            rList = t.split(".")[3].split(",")
            qList = t.split(".")[4].split(",")
            rqGlycs = tuple(([],[]))

            ri = 0
            qi = 0
            rqindex = tuple(({},{}))

            #creates the alignment index again
            '''for n, rchar in enumerate(data[t][0]):
                qchar = data[t][1][n]

                if rchar != '-':
                    rqindex[0].update({ri:n})
                    ri += 1
                if qchar != '-':
                    rqindex[1].update({qi:n})
                    qi += 1'''

            for m in range(2):
                if m == 0:
                    n = 1
                else:
                    n = 0

                if overlap[header][m]:
                    for tpl in overlap[header][m]:
                        ngStart = int(tpl[1].split(":")[0])
                        ngEnd = int(tpl[1].split(":")[1]) + 1


                        inSeq = data[t][n][ngStart:ngEnd]
                        beyond = data[t][n][ngEnd:].replace('-', '')

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
                            #add it to the interfered list
                            duple = False
                            if header in interfered:
                                for x in interfered[header][m]:
                                    if x[0] == tpl[0]:
                                        duple = True
                                if not duple:
                                    interfered[header][m].append(tpl)
                            elif m == 0:
                                interfered[header] = tuple(([tpl], []))
                            else:
                                interfered[header] = tuple(([],[tpl]))



    #print(x)
    #print(total)
    #print(overlap)
    print(interfered)

    #totOut.write("header:seq1:seq2\n")
    #for j in total:
        #totOut.write(str(j) + ":" + ",".join(map(str, total[j][0])) + ":" + ",".join(map(str, total[j][1])) + "\n")
        


    #intOut.write("header:seq1:seq2\n")
    #print("INTERFERED------------")
    for n in interfered:
        out = tuple(([],[]))

        for o in range(2):
            for tup in interfered[n][o]:
                y = tup[0].split(":")[1]
                out[o].append(y)

        #intOut.write(str(n) + ":" + ",".join(map(str, out[0])) + ":" + ",".join(map(str, out[1])) + "\n")










