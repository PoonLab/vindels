from Bio import Phylo
import cStringIO
from glob import glob
import re
import random
import os 


def sample_beast(file, numsample=10):
    tcount = 0
    states = []
    

    path = os.path.dirname(file) + "/"
    name = os.path.basename(file).split(".")[0]
    input = open(file,'rU')

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


    input2 = open(file,'rU')

    seqDict = {}

    for line in input2:
        search = re.search("\d*[^\S\t\n\r]'.*'", line)
        data = line.split()

        if search != None:
            fields = line.strip("\t\n,").split()
            fields[1] = fields[1].strip("'")
            seqDict[fields[0]] = fields[1]
        elif len(data) == 6 and data[1].lstrip("STATE_") in rsample:
            state = data[1].lstrip("STATE_")
            rawtree = data[-1]

            #fixed all comments from the tree lines 
            fixed = re.sub("\[&rate[^\]]*\]", "", rawtree)
            x = fixed.decode('utf-8')

            handle = cStringIO.StringIO(fixed)
            print(type(handle))
            #read the tree as a newick phylo object
            tree = Phylo.read(handle, 'newick')
            #convert all the tip names using the seqDict dictionary

            for tip in tree.get_terminals():
                tip.name = seqDict[tip.name].split(".")[4]
            Phylo.write(tree, path+name+"_"+state+".tree.sample", 'newick')
            


#folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout/*.trees")
#for file in folder:
folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout/*time.trees")

for file in folder:
    sample_beast(file,10)
    break



'''states =[]
for line in infile:
    data = line.split()
    if len(data) == 6:
            states.append(int(data[1].lstrip("STATE_")))
infile.close()

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

Phylo.write(tree, file="")'''


# open file stream to BEAST tree log

# iterate through stream on a by-line basis - pass through TWICE

# first pass through, count the number of trees (record states = MCMC step number)

# determine random sample of trees given user-defined sample size and burn-in

# second pass through

# locate translate block (associates integer indices with sequence labels)

# store index-label combinations in a dictionary

# locate tree block (each Newick tree string is prefixed with "tree").