from seqUtils import *
import subprocess
# TODO: use tempfile module
from glob import glob
from gotoh2 import *

refseq = """MRVKEKYQHLWRWGWRWGTMLLGMLMICSATEKLWVTVYYGVPVWKEATTTLFCASDAKAYDTEVHNVWATHACVPTDPN
PQEVVLVNVTENFNMWKNDMVEQMHEDIISLWDQSLKPCVKLTPLCVSLKCTDLKNDTNTNSSSGRMIMEKGEIKNCSFN
ISTSIRGKVQKEYAFFYKLDIIPIDNDTTSYKLTSCNTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCT
NVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTI
GKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFYCNSTQLFNSTWFNSTW
STEGSNNTEGSDTITLPCRIKQIINMWQKVGKAMYAPPISGQIRCSSNITGLLLTRDGGNSNNESEIFRPGGGDMRDNWR
SELYKYKVVKIEPLGVAPTKAKRRVVQREKR""".replace("\n",'')

refv5 = refseq[456:473]

folder = glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/3_RegionSequences/VRegions_mod/*.csv")

for file in folder:
    print(file)
    csv = open(file,'r')

    name = file.split("/")[-1]

    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/3_RegionSequences/VRegions_modTEST/"+name, 'w')

    for line in csv:
        line = line.split(",")
        v5 = line[5].strip("\n")

        if v5 != "":
            aav5 = translate_nuc(v5,0)

            pairwise = Aligner()
            pairwise.set_model('EmpHIV25')
            pairwise.gap_extend_penalty = 5
            pairwise.gap_open_penalty = 40
            pairwise.is_global = True



            result = pairwise.align(refv5, aav5)

            print(line[0])
            #print(v5)
            print(aav5)

            print(result[0])
            print(result[1])

            newv5 = list(v5)
            for n in range(len(result[1])):

                if result[1][n] == '-':
                    newv5[n * 3:n * 3] = ['-', '-', '-']

            newv5 = "".join(newv5)
            print(newv5)
            ri = 0
            index = {}
            for ai, char in enumerate(result[0]):

                if char != '-':
                    index.update({ri:ai})
                    ri += 1
            #print(result[1][index[2] + 1:index[14]].replace("-", ''))

            final = newv5[(index[2]+1)*3:((index[14])*3)].replace("-",'')   #translate the position, add 1 to make it 1-indexed so that you can multiply by 3, multiply by 3, apply to a 0-index slice, result is the first nt of v5

            print(final)

            print('')

            output.write(",".join(line[0:5]) + ","+ final + "\n")
        else:
            print(line[0])