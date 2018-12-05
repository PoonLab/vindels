from seqUtils import *

fasta_file = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/hiv-db.fasta", 'r')


data = parse_fasta(fasta_file)  #Returns headers and sequences as key/value pairs

subtypes = {}

count = 0
seqcount = 0
#Sort by subtype
for i in data:
    info = i.split(".")
    seqcount += 1
    if info[0] not in subtypes.keys():
        subtypes[info[0]] = {}
        subtypes[info[0]][i] = data[i]
    else:
        subtypes[info[0]][i] = data[i]

print(seqcount)

subcount = 0
for i in subtypes.keys():
    subcount += 1
    for x in subtypes[i].keys():
        header = x.split(".")

        if header[0] == "-" or header[2] == "-" or len(subtypes[i][x]) < 1400:   #Removing sequences without a date sequences & less than 1400 nt
            del subtypes[i][x]

#Delete the subtype if there is no data remaining
    if len(subtypes[i].keys()) == 0:
        del subtypes[i]
print(subcount)

'''
#Testing for sequence lengths
for i in subtypes:
    for x in subtypes[i].keys():
        print(len(subtypes[i][x]))
'''
filtered = 0

#For output
for i in subtypes:

    filtered += 1

    output_file = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/1SubtypeSequences/" + i + "++.fasta", 'w')
    for x in subtypes[i]:
        count += 1

        output_file.write(">"+ x)
        output_file.write("\n")
        output_file.write(subtypes[i][x])
        output_file.write("\n")
        
    output_file.close()

print(filtered)
print(count)



