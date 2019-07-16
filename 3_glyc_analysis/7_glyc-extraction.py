from seqUtils import *
import glob
import os 
import sys 
import re 
import csv

infile = open("/home/jpalmer/PycharmProjects/glyc-analysis/8_aapairwise/gp120.fasta", 'rU')

aapairs = parse_fasta2(infile)

positions = {"acute":{},"chronic":{}}


for header, seq in aapairs.items():

    #if header == "C.HM215294.2003.35293.acute_infection.-.2547":
            
    ref, query = seq 

    # CREATE A CUSTOM ALIGNMENT INDEX FOR THE REFERENCE SEQEUENCE ridx[reference position] = alignment position
    ridx = {}
    ai = 0
    for n, char in enumerate(ref):
        if char != "-":
            ridx.update({n:ai})
            ai += 1

    rnogaps = ref.replace("-","")
    qnogaps = query.replace("-","")

    qidx = alignidx(query)
    
    #print(qidx)
    #print(query)
    nglycs = re.finditer("N[^P][ST][^P]", qnogaps)
    #print(header)
    #print(query)
    status = header.split(".")[4]
    if status == "acute_infection":
        print(header.split(".")[4])
        stage = "acute"
    else:
        stage = "chronic"
    
    for ng in nglycs:
        #print(query[qidx[ng.start()]:qidx[ng.end()-1]+1])

        ng_qry = qidx[ng.start()]
        
        #print(ref[ng_qry-4:ng_qry+8])
        #print(query[ng_qry-4:ng_qry+8])
        #print("")
        if ng_qry in ridx:
            #print(rnogaps[ridx[ng_qry]:ridx[ng_qry]+4])
            positions[stage][ridx[ng_qry]] = positions[stage].get(ridx[ng_qry], 0) + 1


'''with open('/home/jpalmer/PycharmProjects/glyc-analysis/9_glycs/glycs.csv', 'w') as outfile:
    writer = csv.dictWriter(outfile, fieldnames=["Position", "Count"])
    writer.writeheader()
    writer.writerows(positions)'''


outfile = open('/home/jpalmer/PycharmProjects/glyc-analysis/9_glycs/acute.csv', 'w')
outfile.write("position,count\n")
for pos in positions["acute"]:
    outfile.write(str(pos)+","+str(positions["acute"][pos])+"\n")
outfile.close()

outfile = open('/home/jpalmer/PycharmProjects/glyc-analysis/9_glycs/chronic.csv', 'w')
outfile.write("position,count\n")
for pos in positions["chronic"]:
    outfile.write(str(pos)+","+str(positions["chronic"][pos])+"\n")
outfile.close()



        

#MRVTGIRRNYQHL-WGWGTMLIWLLMSCSAAENLWVTVYYGVPVWKEATTTLFCASDAKGYDTEVHNVWATHACVPTDPNPQEIVLGNVTENFNMWKNNMVEQMHEDVISLWDQSLKPCVKLTPLCVTLNCTNWKNNTNTNSTTTSSPESPLEKGEIKNCSFNVTSGIRDKVQKEYALFYRLDIVPIDNENDSFRLISCNTSVITQACPKITFEPIPIHY?APAGFAILK?NNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEDVVIRSSNFSDNAKIIIVHLNETVNITCIRPNNNTRKSINI--GPGRAFYTTGAITGDIRQAHCNLSTSQWNNTLKQIVAKLRGQFGN-RTIVFNHSSGGDPEIVMHSFNCGGEFFYCNTTALFNSTW-NTNGTWKGT--TGENDTITLQCRIKQIINMWQEVGKAMYAPPIRGQIRCTSNITGILLTRDGGKSNNNSNGSEIFRPGGGDMRDNWRSELYKYKVVKIEPLGIAPSKARRRVVQREKR