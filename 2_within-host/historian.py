import os 
import re
import sys 
from glob import * 
from subprocess import check_output, run, PIPE, Popen, STDOUT
import numpy as np


hdir = sys.argv[1]

sizes = os.popen("ls -l " + hdir+'/failed/' + " | awk '{print $5}'").read().split("\n")[1:-1]
sizes = [int(x) for x in sizes]
names = os.popen("ls -l " + hdir+'/failed/' + " | awk '{print $9}'").read().split("\n")[1:-1]


sizes = np.array(sizes)
print(sizes)
print(sizes[np.where(sizes==0)])

names = np.array(names)
outname = names[np.where(sizes==0)]


tree = [x.split("_recon")[0] + ".tree.sample" for x in outname]
msa =  [x.split("-")[0] + ".fasta" for x in outname]


for i in range(len(tree)):
    sys.stdout.write("Beginning analysis on file: " + outname[i]+"\n")
    
    outfile = open(hdir+'output/'+outname[i],'w+')
    
    cmd = ['/home/jpalmer/historian/bin/historian', '-guide',hdir+"msa/"+msa[i], '-tree', hdir+'trees/'+tree[i],'-ancseq','-output','fasta']
    
    proc = run(cmd, check=True, stdout=outfile).stdout
    outfile.close()

