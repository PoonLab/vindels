import re
from glob import glob
import os
from gotoh2 import *
from seqUtils import *
import time

start_time = time.time()
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


'''def pairwiseAlign(ntqry, ntref=""):
    if ntref == "":
        ntref = ATGAGAGTGAAGGAGAAATATCAGCACTTGTGGAGATGGGGGTGGAGATGGGGCACCATGCTCCTTGGGATGTTGATGAT\
CTGTAGTGCTACAGAAAAATTGTGGGTCACAGTCTATTATGGGGTACCTGTGTGGAAGGAAGCAACCACCACTCTATTTT\
GTGCATCAGATGCTAAAGCATATGATACAGAGGTACATAATGTTTGGGCCACACATGCCTGTGTACCCACAGACCCCAAC\
CCACAAGAAGTAGTATTGGTAAATGTGACAGAAAATTTTAACATGTGGAAAAATGACATGGTAGAACAGATGCATGAGGA\
TATAATCAGTTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATT\
TGAAGAATGATACTAATACCAATAGTAGTAGCGGGAGAATGATAATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAAT\
ATCAGCACAAGCATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAATACCAATAGATAATGA\
TACTACCAGCTATAAGTTGACAAGTTGTAACACCTCAGTCATTACACAGGCCTGTCCAAAGGTATCCTTTGAGCCAATTC\
CCATACATTATTGTGCCCCGGCTGGTTTTGCGATTCTAAAATGTAATAATAAGACGTTCAATGGAACAGGACCATGTACA\
AATGTCAGCACAGTACAATGTACACATGGAATTAGGCCAGTAGTATCAACTCAACTGCTGTTAAATGGCAGTCTAGCAGA\
AGAAGAGGTAGTAATTAGATCTGTCAATTTCACGGACAATGCTAAAACCATAATAGTACAGCTGAACACATCTGTAGAAA\
TTAATTGTACAAGACCCAACAACAATACAAGAAAAAGAATCCGTATCCAGAGAGGACCAGGGAGAGCATTTGTTACAATA\
GGAAAAATAGGAAATATGAGACAAGCACATTGTAACATTAGTAGAGCAAAATGGAATAACACTTTAAAACAGATAGCTAG\
CAAATTAAGAGAACAATTTGGAAATAATAAAACAATAATCTTTAAGCAATCCTCAGGAGGGGACCCAGAAATTGTAACGC\
ACAGTTTTAATTGTGGAGGGGAATTTTTCTACTGTAATTCAACACAACTGTTTAATAGTACTTGGTTTAATAGTACTTGG\
AGTACTGAAGGGTCAAATAACACTGAAGGAAGTGACACAATCACCCTCCCATGCAGAATAAAACAAATTATAAACATGTG\
GCAGAAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAGTGGACAAATTAGATGTTCATCAAATATTACAGGGCTGCTAT\
TAACAAGAGATGGTGGTAATAGCAACAATGAGTCCGAGATCTTCAGACCTGGAGGAGGAGATATGAGGGACAATTGGAGA\
AGTGAATTATATAAATATAAAGTAGTAAAAATTGAACCATTAGGAGTAGCACCCACCAAGGCAAAGAGAAGAGTGGTGCA\
GAGAGAAAAAAGA # this is gp120

    ntqry = re.sub("-","",ntqry)

    nt_pair = Aligner()
    nt_pair.set_model('HYPHY_NUC')
    nt_pair.is_global = False
    nt_pair.gap_open_penalty = 30
    nt_pair.gap_extend_penalty = 10

    result = nt_pair.align(ntref, ntqry)

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

    if len(finalRef) != len(finalQry):
        unequal.append(header)
        print(header)
        print(aaref)
        print(aaqry)
        print(finalRef)
        print(finalQry)
        continue'''


d = {}
print(d)
e = {"stuff":[1,2,3], "other":['a','b','c']}
d.update(e)
print(d)



pairwiseAlign("")
gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
gp120 = ""
with open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r') as handle:
    for line in handle:
        line = line.strip("\n")
        gp120 += line

vlad_ntref = gp120[390:]
vlad_aaref = translate_nuc(vlad_ntref, 0)

print(vlad_aaref)

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/1FilteredSeqs/*.fasta")

pairwise = {}
total = 0 
vtotal = 0
subtypes = {}
for infile in folder:

    filename = os.path.basename(infile)
    
    ntref = gp120
    aaref = translate_nuc(gp120, 0)
    
    with open(infile, "r") as handle:
        data = parse_fasta(handle)
    unequal = []

    #output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/2PairwiseAA/"+filename, "w")
    
    # appropriately handle the VN data set 

    for header in data:
        
        if "_" in header:
            st = "C"
            vtotal += 1
        else:
            st = header.split(".")[0]
            total += 1
        subtypes[st] = subtypes.get(st,0) + 1
        if filename == "novitsky.fasta":
            ntqry = data[header].replace("\r","").replace("-","")      
            aaref = vlad_aaref
            ntref = vlad_ntref
        else:
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

        if len(finalRef) != len(finalQry):
            unequal.append(header)
            print(header)
            print(aaref)
            print(aaqry)
            print(finalRef)
            print(finalQry)
            continue

        #output.write(">" + header + '\n')
        #output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')
    #output.close()
    
    

    '''aaqry = translate_nuc(ntqry,0)


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

            print(result2[0])
            print("")
            print(result2[1])
            print("")
            print(finalRef)
            print("")
            print(finalQry)
            print("")
            if len(newref) != len(newqry):
                unequal.append(header)
                print("ERROR")
                #print(header)
                #print(aaref)
                #print(aaqry)
                continue

            #output.write(">" + header + '\n')
            #output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')
        #output.close()'''


print("--- %s Seconds ---" % (time.time() - start_time))
print(total)
print(vtotal)
print(total+vtotal)
print(subtypes)
