import os
from glob import glob
from seqUtils import *


#uses the subtype dictionary.txt file to sort accession numbers + dates into the relevant output file


count = 0

input = open("genbank_dates.txt" , 'r')
dictionary = open("dictionary.txt", 'r')


dates = {}
subtypes = {}


for line in input:
    info = line.strip("\n\r").split(".")

    dates[info[0]] = info[1]

# GENERATED ABOVE:
# dates: keys= ACC#s , items= dates



for n in dictionary:
    data = n.strip("\n\r").split(",")

    if data[2] in subtypes.keys():
        subtypes[data[2]].append(data[0])

    else:
        subtypes[data[2]] = []
        subtypes[data[2]].append(data[0])
#subtypes: keys = Subtype , items= accession #s



for x in subtypes:
    
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Dates/" + x + "_dates.txt", 'w')
    

    for i in subtypes[x]:    #where i is the accession number
        output.write(i + "," + dates[i])
        output.write('\n')

    output.close()


