import re
from glob import glob
from seqUtils import *

folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/10_5_Translated_Cherries/*.fasta")

output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/ngCounts.csv", 'w')
output.write("subtype,V1,V2,V3,V4,V5\n")
for file in sorted(folder):

    fasta = parse_fasta2(open(file,"r"))

    subtype = file.split("/")[-1].split("+")[0]

    aaTotal = {1:0, 2:0, 3:0, 4:0, 5:0}
    ngTotal = {1:0, 2:0, 3:0, 4:0, 5:0}


    for n in fasta:
        vloop = int(n.split(".")[2])

        aaCount = len(fasta[n][0]) + len(fasta[n][1])
        nGlycs1 = re.findall("N[^P][ST][^P]", fasta[n][0])
        nGlycs2 = re.findall("N[^P][ST][^P]", fasta[n][1])

        print(nGlycs1)
        print(nGlycs2)

        ngCount = 4*(len(nGlycs1) + len(nGlycs2))
        print(ngCount)

        '''print(n)
        print(fasta[n][0])
        print(fasta[n][1])
        print(aaCount)
        print(ngCount)'''

        aaTotal[vloop] += aaCount
        ngTotal[vloop] += ngCount
    #print(subtype)
    print(aaTotal)
    print(ngTotal)



    if subtype != 'F1':
        props = [float(ngTotal[1]) / aaTotal[1],
                 float(ngTotal[2]) / aaTotal[2],
                 float(ngTotal[3]) / aaTotal[3],
                 float(ngTotal[4]) / aaTotal[4],
                 float(ngTotal[5]) / aaTotal[5]]
        print(props)
        output.write("{},{},{},{},{},{}\n" .format(subtype, props[0],props[1],props[2],props[3],props[4]))
    else:
        props = [float(ngTotal[1]) / aaTotal[1],
                 float(ngTotal[2]) / aaTotal[2],
                 float(ngTotal[4]) / aaTotal[4],
                 float(ngTotal[5]) / aaTotal[5]]
        print(props)
        output.write("{},{},{},0.0,{},{}\n" .format(subtype, props[0],props[1],props[2],props[3]))
