import os 
import re
import sys 
from glob import * 
from subprocess import check_call, PIPE, Popen, STDOUT
import numpy as np
import time

if len(sys.argv) != 3:
    print("USAGE python historian.py [tree folder] [output folder]")
    sys.exit()

for i in range(len(sys.argv)):
    if not sys.argv[i].endswith('/'):
        sys.argv[i] += '/'


tpath = sys.argv[1]
trees = [os.path.basename(x) for x in glob(tpath+"*")]

opath = sys.argv[2]
mpath = '/home/jpalmer/work/4MSA/final/'

finished = [os.path.basename(x) for x in glob(opath+"*")]
trees = [x for x in trees if x.split(".")[0]+"_recon.fasta" not in finished]

msa =  [mpath+x.split("-")[0] + ".fasta" for x in trees]
out =  [opath+x.split(".tree")[0] + "_recon.fasta" for x in trees]

trees = [tpath+x for x in trees]
start = time.time()

print(len(trees))
