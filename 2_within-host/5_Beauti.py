from glob import glob  
from beauti import parse
import os
import re
from seqUtils import *

p = re.compile('.+_([0-9]+)\..+')

template_file = '/home/jpalmer/vindels/2_within-host/template-constant.xml'

#parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/signal/*.fasta')
subtypes = {}
for f in files:
    with open(f, "r") as handle: 
        data = parse_fasta(handle)
    
    #for header in data:
        #fields = header.split(".")
        #subtypes[fields[0]] = subtypes.get(fields[0],0) + 1

    for rep in ['a','b']:
        name = str(os.path.basename(f)).split(".")[0]  + '-' + rep
        print(name)
        xmlpath = '/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/constant-1e4/'
        out = xmlpath + name + ".xml"
        #os.mkdir(xmlpath+name+'/')
        #stem = xmlpath+name+"/"

        #override stem to test output location
        stem = '/home/jpalmer/11BEAST-constantlog/output/' + name
        print(stem)
        print(out)
        parse(template_file, f, stem, out, 'days', 1)

