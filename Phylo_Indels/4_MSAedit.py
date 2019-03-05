from seqUtils import *
from glob import glob

# MANUAL EDITING
#listed sequences hinder the accurate formation of alignments
#delete the specific sequences provided in the given lists to improve overall alignment quality

cfolder = glob('/home/jpalme56/PycharmProjects/hiv-evolution-master/4_Conserved/*.fasta')
vfolder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/3_RegionSequences/VRegions-final/*.csv")
blacklist =  {"01_AE":['KP411841'],
              "02_AG":['KP411843'],
              "C": ['KU319547','KP411838','MF373131','KU319550','KU319539','MF373138'],
              "B":['KP411824','KP411825','JQ403020','KT427845','DQ339453','KT427832']}

for infile in cfolder:
    fasta = open(infile,"r")
    
    subtype = infile.split("/")[-1].split("_CR")[0]
    data = parse_fasta(fasta)
    if subtype in blacklist:

        for i in data.keys():
            for x in blacklist[subtype]:
                if x in i:
                    print(i)
                    del data[i]

    output_file = open('/home/jpalme56/PycharmProjects/hiv-evolution-master/4_1_Edited/' + subtype + "_CR.fasta", 'w')
    for n in data.keys():
        output_file.write(">"+ n +"\n"+ data[n]+"\n")

    output_file.close()


for infile in vfolder:
    input = open(infile,"r")

    subtype = infile.split("/")[-1].split("_VR")[0]
    output_file = open('/home/jpalme56/PycharmProjects/hiv-evolution-master/3_RegionSequences/VRegions-final/' + subtype + "_VR.csv", 'w')

    print(infile)
    if subtype in blacklist:
        for line in input:
            accno = line.split(",")[0]
            exclude = False
            for x in blacklist[subtype]:
                if x == accno:
                    print(x)
                    print(accno)
                    exclude = True
            if not exclude:
                output_file.write(line)

    else:
        for line in input:
            output_file.write(line)


    output_file.close()



