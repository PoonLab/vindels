from seqUtils import *
from glob import glob

#Converts csv cherries into fasta (ref/query) format
files = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_Cherries/*+.csv")

count = []
for x in files:
    csv = open(x, 'r')

    subtype = x.split("/")[-1].split(".csv")[0]
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_Cherries/" + subtype + "+.fasta", 'w')


    for line in csv:

        data = line.split(",")

        for n, i in enumerate(data):
            data[n] = data[n].strip('"\n')


        #remove all gaps from the sequences
        data[2] = data[2].replace('-', '')
        data[4] = data[4].replace('-', '')


        if "accno1" in data[1]:
            continue

        if len(data[2]) == 0 or len(data[4]) == 0:
            continue






        output.write(">" + data[1] + "." + data[3] + "." + data[5] + "\n>ref\n"+data[2]+"\n>query\n"+data[4]+"\n")
print(count)