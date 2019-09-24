# script to perform a check for hypermutation on all patient MSAs

from glob import glob 
import sys
import os 
import subprocess
from seqUtils import * 
import pandas
import re

msafolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta")
for infile in msafolder:
    
    print(os.path.basename(infile))
    with open(infile, "rU") as handle:
        msa = parse_fasta(handle)
    call = subprocess.check_output(["python","/home/jpalmer/Poplars/poplars/hypermut.py", "--consensus", infile ])

    total, hm = re.split("\n\n", call)

    '''total = re.sub(r"[^\n\S]+(\S)",r"\t\1",total)
    #print(total)
    total = total.split("\n")[2:]
    total = [x.split("\t")[0] for x in total]
    print(len(total))'''

    hm = hm.split("\n")[1:-1]
    hm = [x.split(" ")[0] for x in hm]
    print(len(hm))
    print(hm)

    #output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/"+os.path.basename(infile), "w")
    '''for header in msa:
        if header not in hm:
            output.write(">"+header+"\n"+msa[header]+"\n")'''

