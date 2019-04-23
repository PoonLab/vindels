import re
from glob import glob
import os
from gotoh2 import *
from seqUtils import *



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


gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line

vlad_ntref = ntref[390:]

aaref = translate_nuc(ntref, 0)
vlad_aaref = translate_nuc(vlad_ntref, 0)


folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/1FilteredSeqs/*.fasta")

pairwise = {}

for infile in folder:

    filename = os.path.basename(infile)
    
    
    with open(infile, "r") as handle:
        data = parse_fasta(handle)
    unequal = []
    
    '''if filename != "novitsky.fasta":#print(filename)
        for header in data:

            accno = header.split(".")[4]
            patid = header.split(".")[3]

            nt_pair = Aligner()
            nt_pair.set_model('HYPHY_NUC')
            nt_pair.is_global = False
            nt_pair.gap_open_penalty = 30
            nt_pair.gap_extend_penalty = 10

            result = nt_pair.align(ntref, data[header])

            left, right = get_boundaries(result[0])

            ntqry = result[1][left:right].replace("-","")

            aaqry = translate_nuc(ntqry,0)


            aa_pair = Aligner()
            aa_pair.set_model('EmpHIV25')
            aa_pair.gap_extend_penalty = 10
            aa_pair.gap_open_penalty = 30
            aa_pair.is_global = True

            result2 = aa_pair.align(aaref,aaqry)

            newref = list(ntref)
            newqry = list(ntqry)

            # reads through the amino acid alignment and adds codon gaps to the proper locations
            for i in range(len(result2[0])):
                if result2[0][i] == '-':
                    newref[i * 3:i * 3] = ['-', '-', '-']

                if result2[1][i] == '-':
                    newqry[i * 3:i * 3] = ['-', '-', '-']

            finalRef = "".join(newref)
            finalQry = "".join(newqry)

            if len(newref) != len(newqry):
                unequal.append(header)
                print(header)
                print(aaref)
                print(aaqry)
                continue

            output.write(">" + header + '\n')
            output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')'''
    if filename == "novitsky.fasta":
        print(filename)
        output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/2PairwiseAA/novitsky_.fasta", 'w')
        for header in data:
            
            seq = data[header].replace("\r","").replace("-","")        
            nt_pair = Aligner()
            nt_pair.set_model('HYPHY_NUC')
            nt_pair.is_global = False
            nt_pair.gap_open_penalty = 30
            nt_pair.gap_extend_penalty = 10

            result = nt_pair.align(vlad_ntref, seq)
            #print(result[0])
            #print("")
            #print(result[1])
            #print("")

            #retrieves the gp120 boundaries from the reference sequence 
            left, right = get_boundaries(result[0])

            #retrieves the query sequence within the boundaries dictated by the reference sequence 
            ntqry = result[1][left:right].replace("-","")
            #print(ntqry)
            #print("")
            aaqry = translate_nuc(ntqry,0)


            aa_pair = Aligner()
            aa_pair.set_model('EmpHIV25')
            aa_pair.gap_extend_penalty = 10
            aa_pair.gap_open_penalty = 50
            aa_pair.is_global = True

            result2 = aa_pair.align(vlad_aaref,aaqry)

            #print(result2[0])
            #print("")
            #print(result2[1])
            #print("")

            newref = list(vlad_ntref)
            newqry = list(ntqry)

            # reads through the amino acid alignment and adds codon gaps to the proper locations
            for i in range(len(result2[0])):
                if result2[0][i] == '-':
                    newref[i * 3:i * 3] = ['-', '-', '-']

                if result2[1][i] == '-':
                    newqry[i * 3:i * 3] = ['-', '-', '-']
            
            finalRef = "".join(newref)
            finalQry = "".join(newqry)

            #print("")
            #print(finalRef)
            #print("")
            #print(finalQry)
            #print("")
            if len(newref) != len(newqry):
                unequal.append(header)
                print(header)
                print(aaref)
                print(aaqry)
                continue

            output.write(">" + header + '\n')
            output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')
        output.close()









