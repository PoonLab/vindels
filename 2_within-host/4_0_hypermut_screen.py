# script to perform a check for hypermutation on all patient MSAs

from glob import glob 
import sys
import os 
import subprocess
from seqUtils import * 
import pandas as pd
import re

msafolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/prelim/*.fasta")
full_path = "/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/"
for msafile in msafolder:
    
    filename = os.path.basename(msafile)
    print(filename)
    # opens the sequences to be edited
    with open(msafile, "rU") as handle:
        msa = parse_fasta(handle)

    # remove preceding zeros in all dates 
    dates = []
    for header in msa:
        main, date = header.split("_")
        #print(date)
        dates.append(int(date))
    lowest = min(dates)
    
    df = pd.DataFrame({'headers':list(msa.keys()), 'date': dates, 'sequences':list(msa.values())})
    tp1 = df.where(df['date'] == lowest)
    tp1 = tp1.dropna()
    #print(tp1)
    fasta1 = [[a,b] for a,b in zip(tp1['headers'], tp1['sequences'])]
    #print([x[0] for x in fasta1])
    cnsus = consensus(fasta1).lower()
    
    tpath = "/home/jpalmer/vindels/2_within-host/hm-temp.fasta"
    tfile = open(tpath,"w")
    tfile.write(">REF\n"+cnsus+"\n")
    for header in msa:
        tfile.write(">"+header+"\n"+msa[header]+"\n")
    tfile.close()

    #print(filename)
    print(len(msa))

    # performs the hypermut screen on the prelim MSA 
    call = subprocess.check_output(["python","/home/jpalmer/Poplars/poplars/hypermut.py", tpath ])
    
    # extracts the hypermutated sequences
    total, hm = re.split("\n\n", call)
    #print(total)
    hm = hm.split("\n")[1:-1]
    hm = [x.split(" ")[0] for x in hm]

    total = total.split("\n")[1:]
    total = [x.split(" ")[0] for x in total]
    print(len(total))
    print(len(hm))
    #print(hm)
    #print(len(hm))

    fullfile = open(full_path+filename, "rU")
    fulldata = parse_fasta(fullfile)
    # writes only the non-hypermutated sequences to a screened folder
    output = open(full_path+"hm-screen/"+filename, "w")
    count = 0
    for header in total:
        if header not in hm:
            count += 1
            output.write(">"+header+"\n"+fulldata[header]+"\n")
    print(count)
