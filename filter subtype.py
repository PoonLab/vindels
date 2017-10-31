from seqUtils import *

fasta_file = open("hiv-db.fasta", 'r')

output_file = open("Subtype_Sorted.txt", 'w')


data = parse_fasta(fasta_file)
subtypes = {}

for i in data:
    info = i.split(".")
    if info[0] not in subtypes.keys():
        subtypes[info[0]] = {}
        subtypes[info[0]][i] = data[i]
    else:
        subtypes[info[0]][i] = data[i]

for i in subtypes:
    output_file.write(i)
    output_file.write("\n")

    for x in subtypes[i]:
        output_file.write(x)
        output_file.write("\n")
        output_file.write(subtypes[i][x])
        output_file.write("\n")








