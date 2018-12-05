import re
from glob import glob
from seqUtils import *

folder = glob("/home/jpalmer/PycharmProjects/hiv-evolution-master/3RegionSequences/VRegions_edit/*.csv")

output = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/ngTotal.csv", 'w')
output.write("subtype,vloop,count,accno\n")
for file in folder:

    input = open(file,"r")
    subtype = file.split("/")[-1].split("_VR")[0]

    counter = {1:0,2:0,3:0,4:0,5:0}

    #output.write(subtype + ",")

    length = 0
    for line in input:
        data = line.split(",")\

        for i in range(5):
            aaseq = translate_nuc(data[i+1],0)

            count = len(re.findall("N[^P][ST][^P]", aaseq))

            output.write("{},{},{},{}\n".format(subtype,str(i+1),str(count),data[0]))




    #for n in range(5):
     #   counter[n+1] = float(counter[n+1]) / length

    #print(counter)
    #output.write("{},{},{},{},{}\n".format(counter[1],counter[2],counter[3],counter[4],counter[5]))



