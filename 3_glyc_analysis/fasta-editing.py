from seqUtils import * 

infile = open("/home/jpalmer/PycharmProjects/glyc-analysis/novitsky.fasta", "rU")

fasta = parse_fasta(infile)

for header, seq in fasta.items():
    seq = seq.replace("\n", "")

outfile = open("/home/jpalmer/PycharmProjects/glyc-analysis/novitsky-edited.fasta", "w")
for header, seq in fasta.items():
    if len(header.split("_")[0]) > 2:
        print(header)
    outfile.write(">"+header+"\n"+seq.replace("-","")+"\n")
outfile.close()

