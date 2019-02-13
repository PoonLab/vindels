from seqUtils import *
import sys
fasta_file = open(sys.argv[1], 'r')


data = parse_fasta(fasta_file)  #Returns headers and sequences as key/value pairs

subdict = {}
filtered = {}

count = 0
seqcount = 0
#Sort by subtype
for i in data:
    header = i.split(".")
    seqcount += 1
    if header[0] != "-" and header[2] != "-" and len(data[i][x]) > 1400:
        if header[0] not in subdict.keys():
            subdict[header[0]] = {}
            subdict[header[0]][i] = data[i]
        else:
            subdict[header[0]][i] = data[i]

print(seqcount)

subcount = 0
for i in subdict.keys():
    subcount += 1
    for x in subdict[i].keys():
        header = x.split(".")

        if header[0] != "-" and header[2] != "-" and len(subdict[i][x]) > 1400:   #Removing sequences without a date sequences & less than 1400 nt
            if i not in filtered.keys():
                filtered[i] = {x:subdict[i][x]}
            else:
                filtered[i][x] = subdict[i][x]

print(subcount)

'''
#Testing for sequence lengths
for i in subtypes:
    for x in subtypes[i].keys():
        print(len(subtypes[i][x]))
'''
#filtered = 0

#For output
for i in filtered:

    output_file = open(sys.argv[1] + "/" + i + ".fasta", 'w')
    for x in filtered[i]:
        count += 1

        output_file.write(">"+ x)
        output_file.write("\n")
        output_file.write(filtered[i][x])
        output_file.write("\n")
        
    output_file.close()

print(filtered)
print(count)



