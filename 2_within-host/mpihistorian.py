#!/usr/bin/env python
import os, sys
from subprocess import check_call, Popen, PIPE
import datetime
from glob import *

from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()


#if len(sys.argv) < 2:
#    print 'python mpihistorian.py'


#if not sys.argv[1].endswith('/'):
#    sys.argv[1] += '/'

msa_path = "/home/jpalmer/historian-mcc/msa/"
tree_path = "/home/jpalmer/historian-mcc/trees/"
out_path = "/home/jpalmer/historian-mcc/output/"

files = os.listdir(tree_path)
files = [f for f in files if f.endswith('.tree')]

os.system('export LD_LIBRARY_PATH=/usr/lib/jvm:$LD_LIBRARY_PATH')

for i in range(len(files)):
    if i%nprocs == my_rank and my_rank < len(files):
        msa_name = files[i].split("-")[0]+".fasta"
        out_name = files[i].split(".")[0]+"_recon.fasta"
        # starting terminal message 
        #sys.stdout.write("%s Process %d of %d starting task %d of %d\n"% (datetime.datetime.now(), my_rank, nprocs, i, len(files)))
        #print(msa_path+msa_name)
        print(out_path+out_name)
        #print(tree_path+files[i])
        #p = Popen(" ".join(['historian', '-vv', '-guide', msa_path+msa_name, '-tree', tree_path+files[i],'-ancseq','-output','fasta','>', out_path+out_name]), shell=True)
        outfile = open(out_path+out_name, 'w')
        status = check_call(['historian', '-guide', msa_path+msa_name, '-tree', tree_path+files[i],'-ancseq','-output','fasta'], stdout=outfile)
        outfile.close()

        # completion terminal message
        #sys.stdout.write("%s Process %d of %d completed task %d of %d:\n%s\n"%(datetime.datetime.now(), my_rank, nprocs, i, len(files), result[-1]))

MPI.COMM_WORLD.Barrier()
MPI.Finalize()


