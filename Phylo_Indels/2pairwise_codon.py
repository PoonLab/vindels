from seqUtils import *
from gotoh2 import *
from glob import glob
import os

# USAGE : python [input folder directory] [output directory]
# This script assumes that each sequence contains a header which is formatted in the follow way
# [subtype].[country].[sampling year].[name] ...etc.   according to the standard LANL file output

if len(sys.argv) < 3:
    print 'python [input dir] [output dir]'
    print "This file will perform pairwise alignments between each sequence found in the input dir\
        and HXB2 and output the result to a new dir"
    exit()

#nucleotide version of the reference sequence
gp120 = open("hxb2_gp120_sequence.txt", 'r')

#load the gp120 reference 
ntRef = ''
for line in gp120:
    ntRef += line.strip("\n")

#amino acid version of the reference sequence
aaRef = translate_nuc(ntRef,0)

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

folder = glob(sys.argv[1] + "/*.fasta")

for infile in folder:
    openfile = open(infile, 'r')
    filename = os.path.basename(infile)
    data = parse_fasta(openfile)

    output = open(sys.argv[2] + "/" + filename, 'w')

    stops = []
    unequal = []

    for header in data:

        #nt alignment soley to remove any unwanted gene regions
        nt_pair = Aligner()
        nt_pair.set_model('HYPHY_NUC')
        nt_pair.is_global = False
        nt_pair.gap_open_penalty = 30
        nt_pair.gap_extend_penalty = 10

        result = nt_pair.align(ntRef,data[header])

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

        # testing the amino acid alignment
        print(header)
        print(ntRef)
        print(ntQry)
        print(result2[0])
        print(result2[1])
        print("\n")

        # convert the nucleotide sequences to lists so they can be edited 
        newRef = list(ntRef)
        newQry = list(ntQry)

        #reads through the amino acid alignment and adds codon gaps to the proper locations
        for i in range(len(result2[0])):
            if result2[0][i] == '-':
                newRef[i*3:i*3] = ['-', '-', '-']

            if result2[1][i] == '-':
                newQry[i*3:i*3] = ['-', '-', '-']

        # unequal lengths of final alignment typically indicated a faulty sequence
        if len(newRef) != len(newQry):
            unequal.append(header)
            #print(ntQry)
            #print(ntRef)

            continue

        #aligned ref and query sequences
        finalRef = "".join(newRef)
        finalQry = "".join(newQry)

        print(finalRef)
        print(finalQry)
        print("\n")

        output.write(">" + header + '\n')
        output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')

    #reports the number of faulty stop codon sequences
    print(len(stops))
    #reports the number of unequal sequences
    print(len(unequal))
