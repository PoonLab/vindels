from Bio import Phylo
import cStringIO
from glob import glob
import re
import random
import os 
import sys
import shutil


def getSeqDict(infile):
    #locates the translate table 
    
    seqDict = {}
    lines = open(infile, "rU")

    for line in lines:
        search = re.search("\d*[^\S\t\n\r]'.*'", line)
        data = line.split()

        # this will load the translate dictionary (when search has been found)
        if search != None:
            line = line.strip("\t\n,").split()
            header = line[1].strip("'")
            seqDict[line[0]] = header
    return seqDict

def replaceHeaders(infile, tree, outdir):
    seqDict = getSeqDict(infile)

    fixed = re.sub("\[&rate[^\]]*\]", "", tree)

    handle = cStringIO.StringIO(fixed)
    #print(type(handle))
            #read the tree as a newick phylo object
    tree = Phylo.read(handle, 'newick')
            #convert all the tip names using the seqDict dictionary

    for tip in tree.get_terminals():
        tip.name = seqDict[tip.name]
    Phylo.write(tree, outdir+os.path.basename(infile))



def sample_beast(infile, outdir, numsample=5):
    
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if not outdir.endswith("/"):
        outdir += "/"

    name = os.path.basename(infile).split(".")[0]
    tfile = open(infile,'rU')

    # load a list of state numbers 
    states = []
    for line in tfile:
        data = line.split()
        if len(data) == 6:
            states.append(int(data[1].lstrip("STATE_")))
    tfile.close()
    total = len(states) - 1
    #print(total)
    start = int(total*0.1) + 1
    #print(start)
    
    # Check : upper and lower bounds of sampling 
    #print(states[start:(total+1)])
    #print(states[start])
    rsample = random.sample(states[start:(total+1)], numsample)
    print(len(rsample))
    #rsample = [x for x in rsample]
    
    # Sanity check; should be 0 
    #print([x<20000000 for x in rsample])
    '''for x in range(numsample):
        #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
        random.shuffle(randpool)
        x = randpool.pop()
        rsample.append(str(x*int(states[1])))
    '''
    sample_count = 0
    #print(rsample)
    seqDict = getSeqDict(infile)

    tfile2 = open(infile,'rU')

    #d = open(outdir+"dict/"+name+".dictionary", "w")
    
    for line in tfile2:
        data = line.split()
        # this finds and processes each tree state in the rsample
        if len(data) == 6 and int(data[1].lstrip("STATE_")) in rsample:
            sample_count += 1
            state = data[1].lstrip("STATE_")
            #print(sample_count)
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
        print("USAGE: python 6_sample-beast-trees.py [input trees folder] [output trees folder]")
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

