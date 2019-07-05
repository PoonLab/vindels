import os, sys
import subprocess  
import datetime
from glob import *

msa_path = "/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/"
tree_path = "/home/jpalmer/PycharmProjects/hiv-withinhost/7_5_MCC/rescaled/"
out_path = "/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/test/"

files = glob(tree_path+"*.tree")

for n, i in enumerate(files):
    name = os.path.basename(i)
    msa_name = name.split("-")[0]+".fasta"
    out_name = name.split(".")[0]+"_recon.fasta"
    # starting terminal message 
    #sys.stdout.write("%s Process %d of %d starting task %d of %d\n"% (datetime.datetime.now(), my_rank, nprocs, i, len(files)))
    print(os.path.isfile(msa_path+msa_name))
    print(os.path.isfile(out_path+out_name))
    print(os.path.isfile(tree_path+name))

    
    outfile = open(out_path+out_name, 'w')
    #p = Popen(['historian', '-guide', msa_path+msa_name, '-tree', tree_path+name,'-ancseq','-output','fasta'], stdout=outfile)
    status = subprocess.check_call(['/home/jpalmer/historian/bin/historian', '-v','-guide', msa_path+msa_name, '-tree', tree_path+name,'-ancseq','-output','fasta'], stdout=outfile)
        #'historian', '-guide', msa_path+msa_name, '-tree', tree_path+name,'-ancseq','-output','fasta'],shell=False, stderr=subprocess.STDOUT)
    outfile.close()
    
    if n > 4:
        break


