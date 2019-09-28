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

with open(path+"/111848-final.fasta", "w") as handle:
    handle.write(">REF1\n"+cnsus1+"\n>REF2\n"+cnsus2+"\n")

output = open(path+"/111848-final.fasta", "a+") 
call = subprocess.Popen(["cat", path+"111848-1.fasta"],stdout=output)
call2 = subprocess.Popen(["cat", path+"111848-2.fasta"],stdout=output)


rlist = random.sample(range(0,800),800)
print(rlist)
with open(path+"/800.fasta", "w") as handle:
    for n, header in enumerate()



