from seqUtils import *
from gotoh2 import *
from glob import glob
import os

#nucleotide version of the reference sequence
gp120 = open("hxb2_gp120_sequence.txt", 'r')

refseq = ''
for line in gp120:
    refseq += line.strip("\n")

#amino acid version of the reference sequence
aaRef = translate_nuc(refseq,0)

#used for cutting off flanking regions
def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')

    # returns indices of without gap prefix and suffix
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

    #used to retrieve subtype name from my file name e.g. 01_AE.fasta >> 01_AE
    subtype = filename.split(".")[0]

    data = parse_fasta(fasta_in)

    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/2_2_Test/" + subtype + "++.fasta", 'w')

    stops = []
    unequal = []


    for header in data:

        #nt alignment to remove any unwanted gene regions
        nt_pair = Aligner()
        nt_pair.set_model('HYPHY_NUC')
        nt_pair.is_global = False
        nt_pair.gap_open_penalty = 30
        nt_pair.gap_extend_penalty = 10

        result = nt_pair.align(refseq,data[header])

        #finds the boundaries of the gene of interest
        left, right = get_boundaries(result[0])

        #cuts the query sequence at the boundaries and changes it back to regular nt sequence
        ntQry = result[1][left:right].replace('-','')

        #translate ntQry to amino acids
        aaQry = translate_nuc(ntQry,0)

        #skips any non functional sequences (you can check the genbank ID with the print line)
        if "*" in aaQry:
            stops.append(header)
            #print(header)
            continue


        #amino acid alignment
        aa_pair = Aligner()
        aa_pair.set_model('EmpHIV25')
        aa_pair.gap_extend_penalty = 10
        aa_pair.gap_open_penalty = 30
        aa_pair.is_global = True

        result2 = aa_pair.align(aaRef,aaQry)

        # prints the amino acid alignment
        print(header)
        print(result2[0])
        print(result2[1])
        print("\n")

        # convert the nucleotide sequences to lists so they can be edited (strings cant be mutated)
        newRef = list(refseq)
        newQry = list(ntQry)


        #reads through the amino acid alignment and adds codon gaps to the proper locations
        for i in range(len(result2[0])):
            if result2[0][i] == '-':
                newRef[i*3:i*3] = ['-', '-', '-']

            if result2[1][i] == '-':
                newQry[i*3:i*3] = ['-', '-', '-']


        if len(newRef) != len(newQry):
            unequal.append(header)
            print(ntQry)
            print(refseq)

            continue

        #aligned ref and query sequences
        finalRef = "".join(newRef)
        finalQry = "".join(newQry)

        print(finalRef)
        print(finalQry)
        print("\n")

        #output
        output.write(">" + header + '\n')
        output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')

    #reports the number of faulty stop codon sequences
    print(len(stops))
    #reports the number of unequal sequences
    print(len(unequal))
