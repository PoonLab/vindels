# script to retrieve all indels within the phylogenetic tree
# script will do the following:


# take in the root of the tree

# use a recursive function to iterate over every tree node 
    # if 
    # perform a comparison between each tree node sequence and its ancestor 

from ete2 import *


t = Tree("(A:1,(B:1,(E:1,D:1):0.5):0.5);" )

node = t.get_tree_root()

for leaf in node:
    print(leaf.write())

    