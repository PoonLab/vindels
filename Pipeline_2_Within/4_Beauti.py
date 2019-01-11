from glob import glob
from beauti import parse
import os
import re

p = re.compile('.+_([0-9]+)\..+')
template_file = '../1500/template.xml'

# parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob('../1500/*.sub.fa')

for f in files:
    fn = os.path.basename(f)
    rep = re.sub(p, '\\1', fn)
    stem = 'step' + rep
    ofn = '../1500/{}.xml'.format(stem)
    parse(template_file, f, stem, ofn, '1')
