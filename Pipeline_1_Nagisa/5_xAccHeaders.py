import os
from glob import glob

from seqUtils import *
count=0

#Code for changing all the headers of the MSAlignment files to only accession numbers
files = {}
paths = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/5_2_Shortened/*.fasta")
names = os.listdir("/home/jpalme56/PycharmProjects/hiv-evolution-master/5_2_Shortened/")


for x in paths:

    subtype = x.split("/")[6].split("MSA")[0]
    print(subtype)
    fasta = open(x,'r')

    x = x.strip("/")

    data = parse_fasta(fasta)

    #output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/5_3_AccHeaders/" + subtype + "MSA.fasta" , 'w')

    for i in data.keys():
        accno = i[-8:]

        if "." in accno:
            accno = accno.split(".")[1]


        #output.write(">" + accno + "\n" + data[i] + "\n")




