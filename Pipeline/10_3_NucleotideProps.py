from glob import glob
from seqUtils import *
from collections import defaultdict

def parse_fasta3(handle):
    # Modified parse fasta to return a dictionary of lists containing the reference [0] and the query [1]
    res = {}
    nt = []
    head = ''
    for i in handle:
        if i == "\n":
            continue
        elif i[0] == '>' or i[0] == '#':
            if ">ref" in i:
                continue
            elif ">query" in i:                    #Append the reference sequence
                res.update({head:[[],[]]})
                res[head][0] = nt
                nt = []                            # reset containers
            elif len(nt) > 0 and ">query" not in i:             #Append the query sequence // Occurs when reaching the header of a new sequence
                res[head][1] = nt
                nt = []
                #Parse the new header
                head = i.strip('\n')[1:]

            else:
                head = i.strip('\n')[1:]         #Occurs on the first header
        else:
            nt = i.strip(',\n').upper().split(",")

    res[head][1] = nt
    return res

#Script to generate the nucleotide proportions in indels
indelout = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/indels.csv",'w')
indelout.write("filename,A.1,A.2,A.3,A.4,A.5,C.1,C.2,C.3,C.4,C.5,G.1,G.2,G.3,G.4,G.5,T.1,T.2,T.3,T.4,T.5\n")


folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_2_NtIndels/*.fasta")
dict = {}
for file in sorted(folder):
    input = open(file,"r")

    data = parse_fasta3(input)

    indelProp = {}
    print(file)
    subtype = file.split("/")[-1].split("x")[0].strip("_")

    for i in data:
        print(i)
        print(data[i])
        vregions = i.split(".")[2]
        if data[i][0]:
            for x in data[i][0]:
                indelProp["A."+vregions] = indelProp.get("A."+vregions,0) + int(x.count("A"))
                indelProp["C."+vregions] = indelProp.get("C."+vregions,0) + int(x.count("C"))
                indelProp["G."+vregions] = indelProp.get("G."+vregions,0) + int(x.count("G"))
                indelProp["T."+vregions] = indelProp.get("T."+vregions,0) + int(x.count("T"))
                print(x)
        if data[i][1]:
            for x in data[i][1]:
                indelProp["A." + vregions] = indelProp.get("A."+vregions,0) +  int(x.count("A"))
                indelProp["C." + vregions] = indelProp.get("C."+vregions,0) + int(x.count("C"))
                indelProp["G." + vregions] = indelProp.get("G."+vregions,0) + int(x.count("G"))
                indelProp["T." + vregions] = indelProp.get("T."+vregions,0) + int(x.count("T"))
                print(x)
        #print(indelProp)
        dict[vregions] = dict.get(vregions,0) + 1

    if subtype == "F1":
        indelProp["A.3"] = 0
        indelProp["C.3"] = 0
        indelProp["G.3"] = 0
        indelProp["T.3"] = 0



    indelout.write(subtype + ',')

    for key in sorted(indelProp.keys()):
        print(key)
        indelout.write(str(indelProp[key]))
        if key != "T.5":
            indelout.write(",")
    indelout.write("\n")

'''
folder2 = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/9_0_nonindel/*.csv")
nonout = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/noindel2.csv",'w')
nonout.write("filename,A.1,A.2,A.3,A.4,A.5,C.1,C.2,C.3,C.4,C.5,G.1,G.2,G.3,G.4,G.5,T.1,T.2,T.3,T.4,T.5\n")

for name in sorted(folder2):
    input = open(name, "r")

    counter = {}

    subtype = name.split("/")[-1].split("+")[0]

    for line in input:
        data = line.strip("\n").split(",")
        seq1 = data[2].strip('"')
        seq2 = data[4].strip('"')
        vregion = data[5]
        if seq1 == "seq1":
            continue

        counter["A." + vregion] = counter.get("A."+vregion, 0) + int(seq1.count("A")) + int(seq2.count("A"))
        counter["C." + vregion] = counter.get("C." + vregion, 0) + int(seq1.count("C")) + int(seq2.count("C"))
        counter["G." + vregion] = counter.get("G." + vregion, 0) + int(seq1.count("G")) + int(seq2.count("G"))
        counter["T." + vregion] = counter.get("T." + vregion, 0) + int(seq1.count("T")) + int(seq2.count("T"))


    nonout.write(subtype + ',')

    for key in sorted(counter.keys()):
        print(key)
        nonout.write(str(counter[key]))
        if key != "T.5":
            nonout.write(",")
    nonout.write("\n")
print(dict)
'''