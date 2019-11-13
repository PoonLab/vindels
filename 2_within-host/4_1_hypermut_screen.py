# script to perform a check for hypermutation on all patient MSAs

from glob import glob 
import sys
import os 
import subprocess
from seqUtils import * 
import pandas as pd
import re

msa_path = "/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/prelim/"
msafolder = glob(msa_path+"*.fasta")
full_path = "/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/"

print(msafolder)
with open(msa_path+"111848/blacklist-848.txt","rU") as handle:
    blacklist = handle.readlines()

blacklist = [x.strip("\n") for x in blacklist]

for msafile in msafolder:
    
    filename = os.path.basename(msafile)
    print(filename)
    # opens the sequences to be edited
    with open(msafile, "rU") as handle:
        msa = parse_fasta(handle)

    # extract all dates 
    dates = []
    for header in msa:
        date = header.split("_")[1]
        dates.append(int(date))
    lowest = min(dates)
    
    # creates a data frame with MSA headers, MSA seqs, and extracted dates
    df = pd.DataFrame({'headers':list(msa.keys()), 'date': dates, 'sequences':list(msa.values())})
    print(len(df))
    tp1 = df
    #tp1 = tp1.where(tp1['date'] == lowest)
    #tp1 = tp1.dropna()
    
    print(len(tp1))
    #print(tp1)

    # take the headers/seqs from TIME POINT 1 and load them into a list of [header,seq] lists
    fasta1 = [[a,b] for a,b in zip(tp1['headers'], tp1['sequences'])]
    #print(fasta1)
    cnsus = consensus(fasta1).upper()
    
    # creates a temp MSA file with the consensus sequence at the top
    tpath = "/home/jpalmer/vindels/2_within-host/hm-temp.fasta"
    #if filename in ["111848-1.fasta", "111848-2.fasta"]:
        #tpath = "/home/jpalmer/vindels/2_within-host/"+filename
    tfile = open(tpath,"w+")
    tfile.write(">REF\n"+cnsus+"\n")
    for header in msa:
        tfile.write(">"+header+"\n"+msa[header]+"\n")
    tfile.close()

    # performs the hypermut screen on the prelim MSA 
    call = subprocess.check_output(["python","/home/jpalmer/Poplars/poplars/hypermut.py", tpath ])
    
    # extracts the hypermutated sequences
    total, hm = re.split("\n\n", call)

    # parse the total list of sequences 
    total = total.split("\n")[1:]
    total = [x.split(" ")[0] for x in total]
    
    # parse the list of hypermutated sequences 
    hm = hm.split("\n")[1:-1]
    hm = [x.split(" ")[0] for x in hm]

    # open the full length sequences for parsing  
    fullfile = open(full_path+filename, "rU")
    fulldata = parse_fasta(fullfile)

    if filename in ["111848-1.fasta", "111848-2.fasta"]:
        print(len(total) - len(hm))
        for name in blacklist:
            try:
                total.remove(name)
            except:
                print".",
        print(len(total) - len(hm))

    # write only the non-hypermutated sequences to a screened folder
    output = open(full_path+"hm-screen-full/"+filename, "w+")
    count = 0
    '''if filename=="OS.fasta":
        print(total)
        print(len(total))
        print(hm)
        print(fulldata.keys())
        print(len(fulldata.keys()))'''
    for header in total:
        if header not in hm:
            count += 1
            output.write(">"+header+"\n"+fulldata[header]+"\n")
    print(count)
