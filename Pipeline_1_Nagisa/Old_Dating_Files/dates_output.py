import os
from glob import glob
from seqUtils import *


#uses the subtype dictionary.txt file to sort accession numbers + dates into the relevant output file


count = 0

input = open("all_dates2.txt" , 'r')
subdict = open("dictionary.txt", 'r')


dates = {}
subtypes = {}

for line in input:
    info = line.strip("\n\r").split(".")

    dates[info[0]] = info[1]

#dates: keys= ACC#s , items= dates


for n in subdict:
    data = n.strip("\n\r").split(".")

    if data[1] in subtypes.keys():
        subtypes[data[1]].append(data[0])

    else:
        subtypes[data[1]] = []
        subtypes[data[1]].append(data[0])




for x in subtypes:
    
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Dates/" + x + "_dates.txt", 'w')
    

    for i in subtypes[x]:    #where i is the accession number
        output.write(i + "," + dates[i])
        output.write('\n')

    output.close()


