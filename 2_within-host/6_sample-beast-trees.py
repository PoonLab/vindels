from Bio import Phylo
import cStringIO
from glob import glob
import re
import random
import os 
import sys
import shutil

def sample_beast(infile, outdir, numsample=5):
    
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if not outdir.endswith("/"):
        outdir += "/"

    name = os.path.basename(infile).split(".")[0]
    input = open(infile,'rU')

    # load a list of state numbers 
    states = []
    for line in input:
        data = line.split()
        if len(data) == 6:
            states.append(data[1].lstrip("STATE_"))
    input.close()
    total = len(states) - 1
    #print(total)
    start = int(total*0.1) + 1
    #print(start)
    rsample = []

    randpool = list(range(start, total))
    
    for x in range(numsample):
        #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
        random.shuffle(randpool)
        x = randpool.pop()
        rsample.append(str(x*int(states[1])))
    sample_count = 0

    input2 = open(infile,'rU')

    seqDict = {}


    #d = open(outdir+"dict/"+name+".dictionary", "w")

    for line in input2:

        #locates the translate table 
        search = re.search("\d*[^\S\t\n\r]'.*'", line)
        data = line.split()

        # this will load the translate dictionary (when search has been found)
        if search != None:
            line = line.strip("\t\n,").split()
            header = line[1].strip("'")
            #date = line[1].strip("'").split("_")[1]
            #print(line[0])
            #print(date)
            
            #NEW SEQUENCE HEADER FORMAT
            seqDict[line[0]] = header
        
        # this finds and processes each tree state in the rsample
        elif len(data) == 6 and data[1].lstrip("STATE_") in rsample:
            sample_count += 1
            state = data[1].lstrip("STATE_")
            if len(state) >1:
                state = state[:-5]
            rawtree = data[-1]

            #locates and fixes all the rate comments found within each tree state 
            fixed = re.sub("\[&rate[^\]]*\]", "", rawtree)

            handle = cStringIO.StringIO(fixed)
            #print(type(handle))
            #read the tree as a newick phylo object
            tree = Phylo.read(handle, 'newick')
            #convert all the tip names using the seqDict dictionary

            for tip in tree.get_terminals():
                tip.name = seqDict[tip.name]
                
            Phylo.write(tree, outdir+name+"_"+str(sample_count)+"_"+str(state)+".tree.sample", 'newick')
            #print(tree)


    # used for making dictionaries to convert sequences back to full header (bc i messed up)
    #for num in seqDict.keys():
    #    d.write(seqDict[num].split(".")[4]+","+seqDict[num]+"\n")
    #print(seqDict)


'''if sys.argv[1][-1] != "/":
    sys.argv[1] += "/"
if sys.argv[2][-1] != "/":
    sys.argv[2] += "/"
'''
def main():
    

    if len(sys.argv) != 3:
        print("USAGE: python 6_sample-beast-trees.py [input trees folder] [output trees folder]")
        quit()
    for i in range(len(sys.argv)):
        if not sys.argv[i].endswith("/"):
            sys.argv[i] += "/"

    infolder = glob(sys.argv[1]+"*.time.trees")
    outfolder = sys.argv[2]
    for infile in infolder:
        #filename = os.path.basename(infile).split(".")[0]

        '''
        patfolder = outfolder + filename + "/"
        if not os.path.isdir(patfolder):
            os.makedirs(patfolder)
        else:
            shutil.rmtree(patfolder,ignore_errors=False)
            os.makedirs(patfolder)
        '''
        sample_beast(infile,outfolder, 20)



if __name__ == '__main__': 
    main()





'''
nwk = line.split()[-1]
handle = StringIO(nwk)
tree = Phylo.read(handle, 'newick')
for tip in tree.get_terminals():
    tip.name = dictionary[tip.name]

Phylo.write(tree, infile="")'''