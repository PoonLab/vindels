# Tree modifications 
# used for modifying different attributes of a tree and returning the result
 

require(ape)
intree <- read.tree("~/PycharmProjects/hiv-withinhost/7SampleTrees/30631-a.tree.sample")

#modify all tree edge lengths using the value from the log file 
intree$edge.length <- (intree$edge.length * 2.9042e-5)

#output 
write.tree(intree,"30631-a.tree")

#bifurcating only
intree <- multi2di(intree)

#adding in numeric node labels  (1 to Nnode)
nNode <- intree$Nnode
intree$node.label <- as.character(seq(from=1,to=nNode))

#find all cherries 
ntips <- Ntip(intree)
numtips <- tabulate(intree$edge[which(intree$edge[,2] <= ntips),1])
is.cherry <- which(numtips == 2)

print(is.cherry)
c.nodes <- data.frame()

print(intree$edge)
print(intree$edge.length)
for (node in is.cherry){
	#print(node)
  idx.cherry <- which(intree$edge[,1]==node)
	tips <- intree$edge[idx.cherry,2]
	lengths <- intree$edge.length[idx.cherry]
	#print(tips)
	node.name <- paste0("(",intree$tip.label[tips[1]] , ":",lengths[1],",", intree$tip.label[tips[2]],lengths[2], ")")
	print(node.name)  
	#c.nodes <- rbind(c.nodes, data.frame(name=node, ,]
  
}
  



#write.tree(intree2,'~/historian/data/F1_relabel.tree')
