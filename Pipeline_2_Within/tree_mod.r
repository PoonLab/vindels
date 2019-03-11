# Tree modifications 
# used for modifying different attributes of a tree and returning the result
 

require(ape)
intree <- read.tree("~/historian/data/F1.tree")
nNode <- intree$Nnode
intree$node.label <- as.character(seq(from=1,to=nNode))
intree <- multi2di(intree)

ntips <- Ntip(intree)
numtips <- tabulate(intree$edge[which(intree$edge[,2] <= ntips),1])

is.cherry <- which(numtips == 2)

c.nodes <- data.frame()

for (node in is.cherry){
  idx <- intree$edge[which(intree$edge[,1]==node),2]
  
  node.name <- paste0("(",intree$tip.label[idx[1]] , edgelength, ":", intree$tip.label[idx[2]], edgelength2, ")")
  
  c.nodes <- rbind(c.nodes, data.frame(name=node, ,]
  
}
  



write.tree(intree2,'~/historian/data/F1_relabel.tree')
