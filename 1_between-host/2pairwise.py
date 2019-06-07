from seqUtils import *
import subprocess
# TODO: use tempfile module
import os



gp120 = open("hxb2_gp120_sequence.txt", 'r')

refseq = ''
for line in gp120:
    refseq += line

print(refseq)

#ref = 'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCGGGGGGGGGGGGGTTTCCCCCCCCCCTTTTAAA'


def pairwise(ref, query):
    with open('temp.seq', 'w') as handle:
        handle.write('>ref\n{}\n>query\n{}\n'.format(ref, query))
    
    p = subprocess.check_output(['mafft', '--quiet', 'temp.seq'])
    return p


for filename in os.listdir("/home/jpalme56/PycharmProjects/hiv-evolution-master/1_SubtypeSequences"):

    fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/1_SubtypeSequences/"+filename, 'r')


    #fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/test.fasta", 'r')
    subtype = filename.split(".")[0]
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/2PairwiseAlignments/"+ subtype + "_pairwise.fasta", 'w')

    data = parse_fasta(fasta_in)

    for n in data:
        align = pairwise(refseq, data[n])
        output.write(">" + n + '\n')
        output.write(align)



