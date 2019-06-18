from glob import * 
import sys 
import os 
import csv
from seqUtils import * 

vfolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/abbrev/*.csv")

for infile in vfolder:

    filename = os.path.basename(infile)
    accdict = {}
    dictname = "/home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/prelim_multi/dict/"+filename.split(".")[0]+"-a.dictionary"
    if os.path.isfile(dictname):
        dictfile = open(dictname,"rU")

        # load the accession no dictionary
        for line in dictfile:
            fields = line.strip("\n\r").split(",")
            accdict[fields[0]] = fields[1]



        vfile = open(infile, "rU")
        #reads the header
        vfile.readline()

        outfile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/variable/"+filename,"w")
        outfile.write("header,V1,start,stop,V2,start,stop,V3,start,stop,V4,start,stop,V5,start,stop\n")
        
        
        for line in vfile:
            fields = line.split(",")
            accno = fields[0].split("_")[0]
            
            fields[0] = accdict[accno]

            outfile.write(",".join(fields))
