from Bio import Phylo
import cStringIO
from glob import glob
import re
import random
import os 


def sample_beast(infile, outdir, numsample=10 ):
    tcount = 0
    states = []
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if outdir[-1] != "/":
        outdir += "/"

    name = os.path.basename(infile).split(".")[0]
    input = open(infile,'rU')

    for line in input:
        data = line.split()
        if len(data) == 6:
            states.append(data[1].lstrip("STATE_"))
    input.close()

    total = len(states) - 1
    start = int(total*0.1) + 1

    rsample = []
    for x in range(numsample):
        #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
        rsample.append(str(random.randint(start,total)*int(states[1])))


    input2 = open(infile,'rU')

    seqDict = {}
    

    for line in input2:

        #locates the translate table 
        search = re.search("\d*[^\S\t\n\r]'.*'", line)
        data = line.split()

        # this will load the translate dictionary (when search has been found)
        if search != None:
            line = line.strip("\t\n,").split()
            header = line[1].strip("'").split(".")
            date = line[1].strip("'").split("_")[1]
            print(line[0])
            print(date)
            seqDict[line[0]] = header[4] + "_" + date
        
        # this finds and processes each tree state in the rsample
        elif len(data) == 6 and data[1].lstrip("STATE_") in rsample:
            state = data[1].lstrip("STATE_")
            if len(state) >1:
                state = state[:-5]
            rawtree = data[-1]

            #locates and fixes all the rate comments found within each tree state 
            fixed = re.sub("\[&rate[^\]]*\]", "", rawtree)

            handle = cStringIO.StringIO(fixed)
            print(type(handle))
            #read the tree as a newick phylo object
            tree = Phylo.read(handle, 'newick')
            #convert all the tip names using the seqDict dictionary

            for tip in tree.get_terminals():
                tip.name = seqDict[tip.name]
                
            Phylo.write(tree, outdir+name+"_"+state+".tree.sample", 'newick')
            print(tree)
    print(seqDict)


#folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout/*.trees")
#for infile in folder:
folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout2/*.time.trees")

for infile in folder:
    sample_beast(infile,"/home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/prelim2/", 1)




'''states =[]
for line in ininfile:
    data = line.split()
    if len(data) == 6:
            states.append(int(data[1].lstrip("STATE_")))
ininfile.close()

total = len(states) - 1
start = int(total*0.1) + 1
print(start)
rsample = []

for x in range(10):
    #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
    rsample.append(str(random.randint(start,total)*int(states[1])))

print(rsample)'''
'''for x in range(10):
    #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
    rsample.append(str(random.randint(101,1000)*100000))
print(rsample)'''



        #write the tree using Phylo






'''
nwk = line.split()[-1]
handle = StringIO(nwk)
tree = Phylo.read(handle, 'newick')
for tip in tree.get_terminals():
    tip.name = dictionary[tip.name]

Phylo.write(tree, infile="")'''


# open infile stream to BEAST tree log

# iterate through stream on a by-line basis - pass through TWICE

# first pass through, count the number of trees (record states = MCMC step number)

# determine random sample of trees given user-defined sample size and burn-in

# second pass through

# locate translate block (associates integer indices with sequence labels)

# store index-label combinations in a dictionary

# locate tree block (each Newick tree string is prefixed with "tree").