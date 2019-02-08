from Bio import Phylo
from io import StringIO
from glob import glob
import re
import random


#folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout/*.trees")
#for file in folder:
infile = open("/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout/28465.time.trees",'rU')

rsample = []
for x in range(10):
    #chose 101 just in case (left 101 states for the burn in and sampled from last 900 states out of 1001)
    rsample.append(str(random.randint(101,1000)*100000))
print(rsample)

seqDict = {}
states = []
for line in infile:

    search = re.search("\d*[^\S\t\n\r]'.*'", line)
    data = line.split()
    if search != None:
        
        fields = line.strip("\t\n,").split()
        fields[1] = fields[1].strip("'")
        #print(fields)
        seqDict[fields[0]] = fields[1]
    elif len(data) == 6 and data[1].lstrip("STATE_") in rsample:
        state = data[1].lstrip("STATE_")
        handle = data[-1]

        #remove all comments from the tree lines 
        remove = re.sub("\[&rate[^\]]*\]", "", handle)
        #print(remove)
        break

        #read the tree as a newick phylo object
        tree = Phylo.read(handle, 'newick')
        #convert all the tip names using the seqDict dictionary
        for tip in tree.get_terminals():
            print(tip)

        #write the tree using Phylo


        
        

    
#print(seqDict)
        





'''
nwk = line.split()[-1]
handle = StringIO(nwk)
tree = Phylo.read(handle, 'newick')
for tip in tree.get_terminals():
    tip.name = dictionary[tip.name]

Phylo.write(tree, file="")'''


# open file stream to BEAST tree log

# iterate through stream on a by-line basis - pass through TWICE

# first pass through, count the number of trees (record states = MCMC step number)

# determine random sample of trees given user-defined sample size and burn-in

# second pass through

# locate translate block (associates integer indices with sequence labels)

# store index-label combinations in a dictionary

# locate tree block (each Newick tree string is prefixed with "tree").