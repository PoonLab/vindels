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
            states.append(int(data[1].lstrip("STATE_")))
    input.close()
    total = len(states) * 2 - 1
    #print(total)
    start = int(total*0.1) + 2
    #print(start)
    rsample = []

    #print(states)
    rsample = random.sample(states, numsample)
    rsample = [str(x) for x in rsample]
    '''for x in range(numsample):
        #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
        random.shuffle(randpool)
        x = randpool.pop()
        rsample.append(str(x*int(states[1])))
    '''
    sample_count = 0
    
    #print(rsample)
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
    


def main():
    

    if len(sys.argv) != 3: #and len(sys.argv) != 4:
        print("USAGE: python 6_sample-beast-trees.py [input trees folder] [output trees folder] [optional comma-delimiter exclusion]")
        quit()
    for i in range(len(sys.argv)):
        if not sys.argv[i].endswith("/"):
            sys.argv[i] += "/"

    #if len(sys.argv) == 4:
    #    exclude = sys.argv[4].split(",")
        
    infolder = glob(sys.argv[1]+"*.time.trees")
    outfolder = sys.argv[2]
    
    for infile in infolder:
        #filename = os.path.basename(infile).split(".")[0]
        #skip = False
        #for item in exclude:
        #    if re.search(item, infile) != None:
        #        skip = True
        print(os.path.basename(infile))
        sample_beast(infile,outfolder, 200)



if __name__ == '__main__': 
    main()

