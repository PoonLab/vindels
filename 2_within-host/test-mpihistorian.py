#!/usr/bin/env python
import os, sys
from subprocess import check_call, Popen, PIPE
import datetime
from glob import *
import re 


print(len(sys.argv))
if len(sys.argv) < 2:
    print 'python mpihistorian.py [formatted folder]'


if not sys.argv[1].endswith('/'):
    sys.argv[1] += '/'

msa_path = sys.argv[1]+"msa/"
tree_path = sys.argv[1]+"trees/"
out_path = sys.argv[1]+"output/"

files = os.listdir(tree_path)
files = [f for f in files if f.endswith('.tree')]

#os.system('export LD_LIBRARY_PATH=/usr/lib/jvm:$LD_LIBRARY_PATH')

for i in range(len(files)):
    if i < 4:
        continue
    # for editing one particular file in the data set
    if "original" in files[i]:
        files[i] = re.sub("-original","",files[i])

    # for creating the MSA name
    msa_name = re.sub("-[ab].*$",".fasta", files[i])
    print(msa_name)
    print(files[i])

    # for creating the output file name
    out_name = files[i].split(".")[0]+"_recon.fasta"

    # starting message 
    sys.stdout.write("%s Process starting task %d of %d on file %s\n"%(datetime.datetime.now(), i, len(files),files[i]))
    print(out_name)

            
    outfile = open(out_path+out_name, 'w')
    #p = Popen(" ".join(['historian', '-guide', msa_path+msa_name, '-tree', tree_path+files[i],'-ancseq','-output','fasta','>',out_path+out_name]),shell=True)
    status = check_call(['/home/jpalmer/historian/bin/historian', '-v','-guide', msa_path+msa_name, '-tree', tree_path+files[i],'-ancseq','-output','fasta'], stdout=outfile)
    outfile.close()

    # completion terminal message
    sys.stdout.write("%s Process completed task %d of %d:\n%s\n"%(datetime.datetime.now(), i, len(files), files[i]))



