#!/usr/bin/env python
import os, sys
from subprocess import check_call, Popen, PIPE
import datetime

from mpi4py import MPI
my_rank = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()


if len(sys.argv) < 2:
    print 'python mpibeast.py [path to XML]'

path_to_xml = sys.argv[1]
if not path_to_xml.endswith('/'):
    path_to_xml += '/'

files = os.listdir(path_to_xml)
files = [f for f in files if f.endswith('.xml')]

os.system('export LD_LIBRARY_PATH=/usr/lib/jvm:$LD_LIBRARY_PATH')

for i in range(len(files)):
    if i%nprocs == my_rank and my_rank < len(files):
        sys.stdout.write("%s Process %d of %d starting task %d of %d\n"%(datetime.datetime.now(), my_rank, nprocs, i, len(files)))
        #p = check_call("export LD_LIBRARY_PATH=/usr/lib/jvm; /usr/bin/java -jar -Djava.library.path=/usr/local/lib /usr/local/lib/beast.jar %s%s" % (path_to_xml, files[i]), shell=True)
        p = Popen(['java', '-jar', '-Djava.library.path=/usr/local/lib', '/usr/local/lib/beast.jar', path_to_xml+files[i]], stdout=PIPE)
        #p = Popen(['java', '-jar', '/home/thuy/programs/BEAST/BEASTv1.8.2/lib/beast.jar', path_to_xml+files[i]], stdout=PIPE)
        result = p.communicate()[0].split('\n')
        sys.stdout.write("%s Process %d of %d completed task %d of %d:\n%s\n"%(datetime.datetime.now(), my_rank, nprocs, i, len(files), result[-1]))

MPI.COMM_WORLD.Barrier()
MPI.Finalize()


