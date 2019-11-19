from glob import glob  
from beauti import parse
import os
import re
import sys

from seqUtils import *

p = re.compile('.+_([0-9]+)\..+')

if len(sys.argv) != 3:
    print("USAGE: python 5_Beauti.py [template file name] [unique run id]")
    sys.exit()

template_id = sys.argv[1]
run_id = sys.argv[2]

if os.path.isdir('/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/' + run_id + "/"):
    print("ERROR: Run name already exists ")
    sys.exit()
else:
    os.mkdir('/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/' + run_id + "/")
    
if not os.path.isfile('/home/jpalmer/vindels/2_within-host/' + template_id):
    print("ERROR: Template file not found" + '/home/jpalmer/vindels/2_within-host/' + template_id)
    sys.exit()


template_file = '/home/jpalmer/vindels/2_within-host/'+ template_id

#parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/*.fasta')
subtypes = {}
for f in files:
    with open(f, "r") as handle: 
        data = parse_fasta(handle)
    
    #for header in data:
        #fields = header.split(".")
        #subtypes[fields[0]] = subtypes.get(fields[0],0) + 1

    for rep in ['a','b']:
        name = str(os.path.basename(f)).split(".")[0]  + '-' + rep
        name = re.sub("-original","",name)
        print(name)
        xmlpath = '/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/' + run_id + "/"
        out = xmlpath + name + ".xml"
        #os.mkdir(xmlpath+name+'/')
        #stem = xmlpath+name+"/"

        #override stem to test output location
        stem = '/home/jpalmer/' + run_id + '/output/' + name
        #print(stem)
        print(out)
        #parse(template_file, f, stem, out, 'days', 1)

