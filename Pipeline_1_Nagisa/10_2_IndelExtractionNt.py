from glob import glob
from seqUtils import *

files = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_1_Pairwise_Nt/*++.fastaz")


#Contains counts of nucleotides for all sequences


#List of the accession numbers in the order their counts were added to the noIndel dictionary
seqOrder = []
ntOutput = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/noindel.csv",'w')
ntOutput.write("filename,A.1,A.2,A.3,A.4,A.5,C.1,C.2,C.3,C.4,C.5,G.1,G.2,G.3,G.4,G.5,T.1,T.2,T.3,T.4,T.5\n")
ntOutput2 = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/indel.csv",'w')
ntOutput2.write("filename,A.1,A.2,A.3,A.4,A.5,C.1,C.2,C.3,C.4,C.5,G.1,G.2,G.3,G.4,G.5,T.1,T.2,T.3,T.4,T.5\n")
for x in sorted(files):
    alignment = open(x, 'r')

    data = parse_fasta2(alignment)

    noIndel = {}
    indel = {}

    subtype = x.split("/")[-1].split("++")[0]
    filename = subtype + "_x.fasta"
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_2_NtIndels/" + filename, "w")


    for i in data:

        #Head1 refers to the reference (first) sequence in the pair
        #Head2 refers to the query (second) sequence in the pair
        head1 = i.split(".")[0] + "." + i.split(".")[2]
        head2 = i.split(".")[1] + "." + i.split(".")[2]
        vregion = i.split(".")[2]
        #This skips the template at the top of each file; also skips any problematic pairwise alignments
        if 'accno1.accno2."Vr"' in i or len(data[i]) == 1:
            print(i)
            continue

        #List of tuples containing insertion sequences (0) and their end positions (1)
        rInsert = []
        qInsert= []


        #Temporary strings used to store nucleotides from insertion sequences
        rTemp = ''
        qTemp = ''

        #Temporary counters for the reference (1) and query sequences (2)
        counter = {"A": 0, "C": 0, "G": 0, "T": 0}
        counter2 = {"A": 0, "C": 0, "G": 0, "T": 0}

        #builds the rInsert and qInsert dictionaries containing insertion sequences (keys) and end positions (values)
        for n, rchar in enumerate(data[i][0]):
            qchar = data[i][1][n]

            if rchar != "-" and qchar != "-":
                if rTemp:
                    rInsert.append(tuple((rTemp,n-1)))
                    rTemp = ''
                if qTemp:
                    qInsert.append(tuple((qTemp,n-1)))
                    qTemp = ''
                
                #NON-INDEL NT Counts-----------
                if (rchar in "ATCG") and (qchar in "ATCG"):
                    counter[rchar] = counter.get(rchar,0) + 1
                    counter[qchar] = counter.get(qchar,0) + 1
                #--------------------

            # Case 2: query sequence has an insertion
            elif rchar == "-" and qchar != "-":
                qTemp = qTemp + qchar

                if rTemp:
                    rInsert.append(tuple((rTemp, n - 1)))
                    rTemp = ''
                if qchar in "ATCG":
                    counter2[qchar] = counter2.get(qchar, 0) + 1
            # Case 3: reference sequence has an insertion
            elif rchar != "-" and qchar == "-":
                rTemp = rTemp + rchar
                if qTemp:
                    qInsert.append(tuple((qTemp, n - 1)))
                    qTemp = ''

                if rchar in "ATCG":
                    counter2[rchar] = counter2.get(rchar, 0) + 1
        if rTemp:
            rInsert.append(tuple((rTemp, n - 1)))
            rTemp = ''
        if qTemp:
            qInsert.append(tuple((qTemp, n - 1)))
            qTemp = ''


        #Nucleotide Counts
        #Transfer the contents of counters 1 & 2 to the noIndel dictionary


       #INSERTION OUTPUT: Used for generating the insertion sequences found in folder 10_2
       #E.g. >HEADER \n >ref \n ACGTCGA,ACCAGA,AGT,
        #print(rInsert)
        #print(qInsert)
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


        #LOAD THE NOINDEL DICTIONARY
        for j in counter.keys():
            label = j+'.'+vregion

            noIndel[label] = noIndel.get(label,0) + counter[j]

        for k in counter2.keys():
            label = k +"."+vregion
            indel[label] = indel.get(label,0) + counter2[k]

        #seqOrder.append(head2)
        
        if subtype =="F1":
            noIndel["A.3"] = 0
            noIndel["C.3"] = 0
            noIndel["G.3"] = 0
            noIndel["T.3"] = 0
            indel["A.3"] = 0
            indel["C.3"] = 0
            indel["G.3"] = 0
            indel["T.3"] = 0
        
        #print(counter)
        #print(counter2)
        if vregion == "1":
            print(counter)
            print(vregion)
            print("----")

        #print(key)
    ntOutput.write(subtype + ',')
    ntOutput2.write(subtype + ',')

    for key in sorted(noIndel.keys()):
        #print(subtype)
        #print(key)
        ntOutput.write(str(noIndel[key]))
        if key != "T.5":
            ntOutput.write(",")
    ntOutput.write("\n")

    for key in sorted(indel.keys()):
        #print(subtype)
        #print(key)
        ntOutput2.write(str(indel[key]))
        if key != "T.5":
            ntOutput2.write(",")
    ntOutput2.write("\n")
    print(noIndel)
    print(indel)

















