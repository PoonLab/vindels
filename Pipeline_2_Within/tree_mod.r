# Tree modifications 
# used for modifying different attributes of a tree and returning the result
 

require(ape)
intree <- read.tree("~/historian/data/F1.tree")
nNode <- intree$Nnode
intree$node.label <- as.character(seq(from=1,to=nNode))
intree2 <- multi2di(intree)


write.tree(intree2,'~/historian/data/F1_relabel.tree')
