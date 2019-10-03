import os 
from seqUtils import *
import subprocess
import random

'''
path = "/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/111848/subsample/"

with open(path+"111848.fasta", "rU") as handle:
    data = parse_fasta(handle)


rlist = random.sample(range(0,1031),600)
names = [x for n, x in enumerate(data.keys()) if n in rlist]

output = open(path+"subset.fasta","w+")
for header in names:
    output.write(">"+header+"\n"+data[header].replace("-","")+"\n")
output.close()

# realign this output file 

# split it up into two files

# read files in here 

with open(path+"subset-1.fasta", "rU") as handle:
    msa1 = convert_fasta(handle)
with open(path+"subset-2.fasta", "rU") as handle:
    msa2 = convert_fasta(handle)
# determine consensus sequence of each 
con1 = consensus(msa1)
con2 = consensus(msa2)
# add consensus sequences and all remaining sequences to master file 
with open(path+"subset-final.fasta", "w+") as handle:
    handle.write(">REF1\n"+con1+"\n>REF2\n"+con2+"\n")

write_fasta_c(msa1,path+"subset-final.fasta","a+")
write_fasta_c(msa2,path+"subset-final.fasta","a+")
infile = open(path+"subset-final.fasta","rU")
for line in infile:
    if len(line) != 1613:
        print(line)
infile.close()
'''
path = "/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/111848/"

with open(path+"111848-1.fasta","rU") as handle:
    msa1 = convert_fasta(handle)

with open(path+"111848-2.fasta","rU") as handle:
    msa2 = convert_fasta(handle)
    
cnsus1 = consensus(msa1)
cnsus2 = consensus(msa2)

print(cnsus1)
print(cnsus2)

with open(path+"111848-final.fasta", "w") as handle:
    handle.write(">REF1\n"+cnsus1+"\n>REF2\n"+cnsus2+"\n")

write_fasta_c(msa1,path+"111848-final.fasta","a+")
write_fasta_c(msa2,path+"111848-final.fasta","a+")

'''
# TESTING 
# ----------------------------------------

with open(path+"111848-1.fasta", "rU") as handle:
    data1 = parse_fasta(handle)
with open(path+"111848-2.fasta","rU") as handle:
    data2 = parse_fasta(handle)

cf1 = [[a,b] for a, b in zip(data1.keys(), data1.values())]
cf2 = [[a,b] for a, b in zip(data2.keys(), data2.values())]


data1.update(data2)
rlist = random.sample(range(0,1031),600)
names = [x for n, x in enumerate(data1.keys()) if n in rlist]
sub1 = [x for x in cf1 if x[0] in names ]
sub2 = [x for x in cf2 if x[0] in names ]

#print(len(names))

con1 = consensus(sub1)
con2 = consensus(sub2)

test = open(path+"800.fasta","w+")

#print(rlist)
test.write(">REF1\n"+con1.upper()+"\n>REF2\n"+con2.upper()+"\n")
for n, header in enumerate(data1):
    if n in rlist:
        test.write(">"+header+"\n"+data1[header].replace("-","")+"\n")


out = open(path+"800-final.msa","w+")
out.write(">REF1\n"+cnsus1.upper()+"\n>REF2\n"+cnsus2.upper()+"\n")
out.close()
with open(path+"800.msa", "rU") as msa:
    data = parse_fasta(msa)
    write_fasta(data, path+"800-final.msa", "a+")
out.close()



infile = open(path+"111848-final.fasta","rU")
for line in infile:
    print(len(line))'''

