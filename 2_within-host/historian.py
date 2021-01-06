import os 
import re
import sys 
from glob import * 
from subprocess import check_call, PIPE, Popen, STDOUT
import numpy as np
import time

if not sys.argv[1].endswith('/'):
    sys.argv[1] += '/'


tpath = sys.argv[1]+'trees/'
opath = sys.argv[1]+'output/'
mpath = '/home/jpalmer/work/4MSA/final/'

trees = glob(tpath+'*')


msa =  [mpath+os.path.basename(x.split("-")[0] + ".fasta") for x in trees]
out =  [opath+os.path.basename(x.split(".tree")[0] + "_recon.fasta") for x in trees]

start = time.time()

for i in range(len(trees)):
    print("Beginning analysis on file: " + trees[i]+"\n")
    #sys.stdout.write(msa[i]+'\n')
    #sys.stdout.write(out[i]+'\n')

    outfile = open(out[i],'w+') 
    status = check_call(['/home/jpalmer/historian/bin/historian', '-guide', msa[i], '-tree', trees[i],'-ancseq','-mcmc','-output','fasta'], stdout=outfile) 
    outfile.close()
    current = time.time()
    print('Elapsed: {}'.format(current-start))

