from seqUtils import *
import sys
# USAGE : python [input file name] [output directory]
# This script assumes that each sequence contains a header which is formatted in the follow way
# [subtype].[country].[sampling year].[name] ...etc.   according to the standard LANL file output

if len(sys.argv) < 3:
    print 'python [input fasta file] [output directory]'
    print "This file will a) filter out sequences that do not contain a collection year and sample\
    date and are less than 1400 nt, and b) sort them based on subtype"
    



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
    if header[0] != "-" and header[2] != "-" and len(data[i]) > 1400:
        if header[0] not in subdict.keys():
            subdict[header[0]] = {}
            subdict[header[0]][i] = data[i]
        else:
            subdict[header[0]][i] = data[i]

print(seqcount)

subcount = 0

#output each subtype into a different file
for subtype in subdict:
    subcount+=1
    if len(subdict[subtype]) == 0:
        print(subtype + " is empty")
        continue

    output_file = open(sys.argv[2] +"/" +subtype + ".fasta", 'w')
    for seq in subdict[subtype]:
        count += 1

        output_file.write(">"+ seq)
        output_file.write("\n")
        output_file.write(subdict[subtype][seq])
        output_file.write("\n")
        
    output_file.close()

print(subcount)
print(count)



