import sys
from glob import glob
from seqUtils import *
import os 

def removeGaps(fasta, outdir):
    name = os.path.basename(fasta)

    alignment = open(fasta, "r")
    
    transposed = transpose_fasta(convert_fasta(alignment))

    whitelist = []
    for pos, x in enumerate(transposed):

        gaps = x.count("-")

        freq = float(gaps)/len(x)

        if freq < 0.95:
            whitelist.append(pos)

    alignment.close()
    
    alignment = open(fasta, "r")

    data = parse_fasta(alignment)

    output = open(outdir + name,'w')
    
    for header in data.keys():
        newseq = ''
        for n, char in enumerate(data[header][:]):

            if n in whitelist:
                newseq += char

        output.write('>' + header + "\n")
        output.write(newseq + "\n")


    alignment.close()
    output.close()

if not sys.argv[1].endswith("/"):
    sys.argv[1] += "/"
if not sys.argv[2].endswith("/"):
    sys.argv[2] += "/"
infolder = glob(sys.argv[1] + "*.fasta")
outfolder = sys.argv[2]

for infile in infolder:
    removeGaps(infile, outfolder)

