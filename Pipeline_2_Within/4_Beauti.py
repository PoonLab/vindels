from glob import glob
from beauti import parse
import os
import re

p = re.compile('.+_([0-9]+)\..+')

template_file = 'template.xml'

#parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta')

for f in files:
    if "101827" in f:
        name = str(os.path.basename(f)).split(".fasta")[0]
        newpath = '/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/'

        #os.mkdir(newpath+name+'/')
        stem = newpath+name+"/"

        ofn = stem + name + ".xml"
        parse(template_file, f, stem, ofn, 'days', 1)
        break
