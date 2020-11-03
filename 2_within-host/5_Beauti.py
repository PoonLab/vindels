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

work = "/home/jpalmer/PycharmProjects/hiv-withinhost/"
template_path = '/home/jpalmer/vindels/2_within-host/templates/'
out_path = work+'5BEAST/'

template_id = sys.argv[1]
run_id = sys.argv[2]


template_file = template_path+template_id
# Check if template file exists 
if not os.path.isfile(template_file):
    print("ERROR: Template file not found " + template_file)
    sys.exit()

# Check for existing directory 
if os.path.isdir(out_path + run_id + "/"):
    print("ERROR: Run name already exists ")
    sys.exit()
else:
    os.mkdir(out_path + run_id + "/")
    


#parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob(work+'4MSA/final/*.fasta')
subtypes = {}
for f in files:
    with open(f, "r") as handle: 
        data = parse_fasta(handle)
    
    #for header in data:
        #fields = header.split(".")
        #subtypes[fields[0]] = subtypes.get(fields[0],0) + 1

    for rep in ['a','b']:
        name = str(os.path.basename(f)).split(".")[0]  + '-' + rep
        #name = re.sub("-original","",name)
        print(name)
        xmlpath = out_path + run_id + "/"
        out = xmlpath + name + ".xml"
        #os.mkdir(xmlpath+name+'/')
        #stem = xmlpath+name+"/"

        #override stem to test output location
        stem = '/home/jpalmer/' + run_id + '/output/' + name
        #print(stem)
        print(out)
        parse(template_file, f, stem, out, 'days', 1)

