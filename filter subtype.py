from seqUtils import *

fasta_file = open("hiv-db.fasta", 'r')


data = parse_fasta(fasta_file)  #Returns headers and sequences as key/value pairs

subtypes = {}

count = 0

#Sort by subtype
for i in data:
    info = i.split(".")
    if info[0] not in subtypes.keys():
        subtypes[info[0]] = {}
        subtypes[info[0]][i] = data[i]
    else:
        subtypes[info[0]][i] = data[i]


#Removing sequences less than 1400 nt
for i in subtypes.keys():
    for x in subtypes[i].keys():
        if len(subtypes[i][x]) < 1400:
            del subtypes[i][x]

#Delete the subtype if there is no data remaining
    if len(subtypes[i].keys()) == 0:
        del subtypes[i]



'''
#Testing for sequence lengths
for i in subtypes:
    for x in subtypes[i].keys():
        print(len(subtypes[i][x]))
'''


#For output
for i in subtypes:
    output_file = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Subtype Sequences/"+ i, 'w')

    for x in subtypes[i]:
        output_file.write(x)
        output_file.write("\n")
        output_file.write(subtypes[i][x])
        output_file.write("\n")


print(count)



