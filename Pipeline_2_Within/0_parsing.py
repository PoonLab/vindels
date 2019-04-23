from seqUtils import * 
import sys 
import glob

file = open("/home/jpalmer/Downloads/21AR_env.fas","r")

data = parse_fasta(file)
count = 0
set1 = set()
set2 = set()
set3 = set()
for header in data:
    #print(header)
    fields = header.split("_")
    set1.add(fields[0])
    set2.add(fields[1])
    set3.add(fields[2])
    seq = data[header].replace("\r","")
    seq = seq.replace("-","")
    print(seq)

    count+=1
print(set1)
print(set2)
print(set3)

print(count)