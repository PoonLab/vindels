from seqUtils import *
import subprocess
# TODO: use tempfile module
import os
import re
from gotoh2 import *
from operator import itemgetter
gp120 = open("hxb2_gp120_sequence.txt", 'r')

refseq = ''
for line in gp120:
    refseq += line.strip("\n")

aaRef = translate_nuc(refseq,0)


def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')

    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0, len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])

    if right:
        res[1] = len(str) - len(right[0])

    return res

for filename in os.listdir("/home/jpalme56/PycharmProjects/hiv-evolution-master/1SubtypeSequences"):
    fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/1SubtypeSequences/" + filename, 'r')


    #fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/test.fasta", 'r')
    subtype = filename.split(".")[0]


    data = parse_fasta(fasta_in)


    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/2_1_AAPairwise/" + subtype + "++.fasta", 'w')



    counter = []
    incorrect = []
    for n in data:
        #print(n)

        nucSeq = {}


        nt_pair = Aligner()
        nt_pair.set_model('HYPHY_NUC')
        nt_pair.is_global = False
        nt_pair.gap_open_penalty = 30
        nt_pair.gap_extend_penalty = 10

        result = nt_pair.align(data[n], refseq)

        left, right = get_boundaries(result[1])

        ntQry = result[0][left:right].replace('-','')


        #print("&&&&")
        '''for offset in range(3):
            temp = translate_nuc(ntQry,offset)
            nucSeq[offset] = temp
            counter.append(temp.count("*"))'''

        aaQry = translate_nuc(ntQry,0)
        if "*" in aaQry:
            incorrect.append(n)
            continue

        aa_pair = Aligner()
        aa_pair.set_model('EmpHIV25')
        aa_pair.gap_extend_penalty = 10
        aa_pair.gap_open_penalty = 30
        aa_pair.is_global = True


        result2 = aa_pair.align(aaQry,aaRef)

        newQry = list(ntQry)
        newRef = list(refseq)
        #print(n)
        #print(result2[0])
        #print(result2[1])


        for i in range(len(result2[1])):
            if result2[0][i] == '-':
                newQry[i*3:i*3] = ['-', '-', '-']

            if result2[1][i] == '-':
                newRef[i*3:i*3] = ['-', '-', '-']


        if len(newRef) != len(newQry):
            incorrect.append(n)
            continue
        #print("".join(newQry))
        #print("".join(newRef))



        output.write(">" + n + '\n')
        output.write(">ref\n" + ''.join(newRef) + "\n>query\n" + ''.join(newQry) + '\n')
    print(filename)
    print(len(incorrect))