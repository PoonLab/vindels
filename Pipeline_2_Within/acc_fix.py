from glob import * 
import sys 
import os 
import csv
from seqUtils import * 

ifolder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/insertions/*.csv")
dpath = "/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/deletions/"

for infile in ifolder:

    filename = os.path.basename(infile)
    accdict = {}
    dictfile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/dictionaries/"+filename.split(".")[0]+".dictionary","rU")

    for line in dictfile:
        fields = line.strip("\n\r").split(",")
        accdict[fields[0]] = fields[1]


    icsv = open(infile, "rU")
    icsv.readline()

    insout = open("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/ins_fix/"+filename,"w")
    insout.write("Accno,Ins,Vloop,Vlen,Seq\n")
    for line in icsv:
        fields = line.split(",")
        accno = fields[0].split("_")[0]
        
        fields[0] = accdict[accno]

        insout.write(",".join(fields))


    dcsv = open(dpath+filename,"rU")
    dcsv.readline()
    delout = open("/home/jpalmer/PycharmProjects/hiv-withinhost/9Indels/del_fix/"+filename,"w")
    delout.write("Accno,Del,Vloop,Vlen,Seq\n")

    for line in dcsv:
        fields = line.split(",")
        accno = fields[0].split("_")[0]
        
        fields[0] = accdict[accno]

        delout.write(",".join(fields))