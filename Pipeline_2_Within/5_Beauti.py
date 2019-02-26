from glob import glob
from beauti import parse
import os
import re

p = re.compile('.+_([0-9]+)\..+')

template_file = '/home/jpalmer/vindels/Pipeline_2_Within/template.xml'

#parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta')

for f in files:

    for rep in ['a','b']:
        name = str(os.path.basename(f)).split(".")[0] + '-' + rep
        print(name)
        xmlpath = '/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/'
        out = xmlpath + name + ".xml"
        #os.mkdir(xmlpath+name+'/')
        #stem = xmlpath+name+"/"

        #override stem to test output location
        stem = '/home/jpalmer/6BEASTout/' + name
        print(stem)
        print(out)
        parse(template_file, f, stem, out, 'days', 1)

