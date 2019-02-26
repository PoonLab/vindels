import os
from glob import glob
from seqUtils import *


#makes a dictionary that matches accession numbers to their subtype

output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/dictionary.txt", 'w')



for file in glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/5MSAs/*.msa"):


    fasta = open(file, "r")

    data = parse_fasta(fasta)

    for i in data:
        accno = i[-8:]
        if "." in accno:
            accno = accno.split(".")[1]

        date = i.split(".")[2]
        subtype = i.split(".")[0].strip(">")


        output.write(accno + ","+date+ ","+ subtype)
        output.write('\n')





