import os 
from seqUtils import *
import subprocess
import random

path = "/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/111848/"
with open(path+"111848-1.fasta","rU") as handle:
    msa1 = convert_fasta(handle)

with open(path+"111848-2.fasta","rU") as handle:
    msa2 = convert_fasta(handle)
    
cnsus1 = consensus(msa1).lower()
cnsus2 = consensus(msa2).lower()

print(cnsus1)
print(cnsus2)

with open(path+"111848-final.fasta", "w") as handle:
    handle.write(">REF1\n"+cnsus1+"\n>REF2\n"+cnsus2+"\n")

output = open(path+"111848-final.fasta", "w+") 
call = subprocess.Popen(["cat", path+"111848-1.fasta"],stdout=output)
call2 = subprocess.Popen(["cat", path+"111848-2.fasta"],stdout=output)
output.close()

# TESTING 
# ----------------------------------------

with open(path+"111848-1.fasta", "rU") as handle:
    data1 = parse_fasta(handle)
with open(path+"111848-2.fasta","rU") as handle:
    data2 = parse_fasta(handle)

cf1 = [[a,b] for a, b in zip(data1.keys(), data1.values())]
cf2 = [[a,b] for a, b in zip(data2.keys(), data2.values())]


data1.update(data2)
rlist = random.sample(range(0,1031),798)
names = [x for n, x in enumerate(data1.keys()) if n in rlist]
sub1 = [x for x in cf1 if x[0] in names ]
sub2 = [x for x in cf2 if x[0] in names ]

print(len(names))

con1 = consensus(sub1)
con2 = consensus(sub2)

test = open(path+"800.fasta","w+")

#print(rlist)
test.write(">REF1\n"+con1.upper()+"\n>REF2\n"+con2.upper()+"\n")
for n, header in enumerate(data):
    if n in rlist:
        test.write(">"+header+"\n"+data[header].replace("-","")+"\n")


out = open(path+"800-final.msa","w+")
out.write(">REF1\n"+cnsus1.upper()+"\n>REF2\n"+cnsus2.upper()+"\n")
out.close()
with open(path+"800.msa", "rU") as msa:
    data = parse_fasta(msa)
    write_fasta(data, path+"800-final.msa", "a+")
out.close()





