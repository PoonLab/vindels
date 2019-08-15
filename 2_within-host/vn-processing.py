from seqUtils import * 

infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/novitsky.fasta", "rU")

fasta = parse_fasta(infile)

patDict = {}



for header,seq in fasta.items():
    
    # fix the sequences 
    seq = seq.replace("\n", "")
    seq = seq.replace("-","")
    #print(seq)

    # extract the fields from the header 
    pat, date, id = header.split("_")[0:3]
    
    '''
    # for glyc-analysis 
    if len(pat) == 1: 
        status = "acute_infection"
    else:
        status = "chronic"
    '''

    # form a new header 
    newheader = "C.-.-."+pat+"."+id+".-.-_"+date
    #print(pat)
    #print(newheader)
    #outfile.write(">"+ newheader + "\n" + seq + "\n" )

    if id not in patDict.keys():
        patDict[id] = {}
    patDict[id][newheader] = seq
total = 0
for pat in patDict.keys():
    unique = set()
    total += len(patDict[pat])
    for header in patDict[pat]:
        date = header.split("_")[1]
        unique.add(date)
    
    if len(unique) < 5: 
        del patDict[pat]
    else:
        print(pat)
print(len(patDict))
print(total)

for pat in patDict.keys():
    outfile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/test/"+pat+".fasta", "w")
    #print(len(patDict[pat]))
    for header,seq  in patDict[pat].items():
        #print(header)
        outfile.write(">"+ header + "\n" + seq + "\n" )
    outfile.close()





    